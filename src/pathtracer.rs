use std::io::{stdout, Write};
use std::sync::{Arc, Mutex};
use std::thread::{JoinHandle, spawn};
use time::now;

use physicsscene::PhysicsScene;
use pngfile::{make_file, read_frame, read_scene, write_frame, write_scene};
use renderer::Renderer;
use rendererscene::RendererScene;
use rendereroutput::RendererOutput;
use rendereroutputrow::RendererOutputRow;
use rhf::RHF;

pub struct Pathtracer {
	width: u32,
	height: u32,
	frame_number: u32,
	number_of_threads: usize,
	continue_old_simulation: bool,
	post_process_data: Vec<(f64, i32)>, // max_distance and patch_size
	search_window_radius: i32,
	spp_per_iteration: u32,
	maximum_spp: u32,
	maximum_error: f64,
	maximum_brdf_value: f64,
	perform_post_process: bool,
}

impl Pathtracer {
	pub fn new(width: u32, height: u32, frame_number: u32, number_of_threads: usize, continue_old_simulation: bool, post_process_data: Vec<(f64, i32)>, search_window_radius: i32, spp_per_iteration: u32, maximum_spp: u32, maximum_error: f64, maximum_brdf_value: f64, perform_post_process: bool) -> Self {
		Self {
			width,
			height,
			frame_number,
			number_of_threads,
			continue_old_simulation,
			post_process_data,
			search_window_radius,
			spp_per_iteration,
			maximum_spp,
			maximum_error,
			maximum_brdf_value,
			perform_post_process,
		}
	}

	pub fn only_post_process(&self) {
		let pixeldata_filename = "/Users/christian/video/pixeldata.txt".to_string();
		let read_result = read_frame(&pixeldata_filename);
		match read_result {
			Err(e) => {
				println!("{:?}", e);
			}
			Ok(mut renderer_output) => {
				self.write_to_image(&mut renderer_output, 0);
				for (i, post_process_datum) in self.post_process_data.iter().enumerate() {
					self.post_process(&mut renderer_output, *post_process_datum);
					self.write_to_image(&mut renderer_output, i+1);
				}
			}
		}
	}

	pub fn render_images_to_png(&mut self) {
		// Either start from an already created, stored scene, or create a new one.
		let mut physics_scene = PhysicsScene::new();
		if self.continue_old_simulation {
			self.read_old_physics_scene(&mut physics_scene);
		} else {
			self.frame_number = 0;
		}
		let mut renderer_scene = RendererScene::new(physics_scene.clone());
		loop {
			println!("Frame {}.", self.frame_number);
			self.write_and_read_scene(&mut physics_scene, &mut renderer_scene);
			let mut renderer_output = self.render(&mut renderer_scene);
			self.write_to_image(&mut renderer_output, 0);
			let pixeldata_filename = "/Users/christian/video/pixeldata.txt".to_string();
			let write_result = write_frame(&pixeldata_filename, &renderer_output);
			if let Err(e) = write_result {
				println!("{:?}", e);
			}
			if self.perform_post_process {
				for (i, post_process_datum) in self.post_process_data.iter().enumerate() {
					self.post_process(&mut renderer_output, *post_process_datum);
					self.write_to_image(&mut renderer_output, i+1);
				}
			}
			self.simulate_physics(&mut physics_scene);
			self.frame_number += 1;
		}
	}

	fn read_old_physics_scene(&self, physics_scene: &mut PhysicsScene) {
		let scenedata_filename = format!("/Users/christian/video/{:04}.txt", self.frame_number);
		let scene_result = read_scene(&scenedata_filename);
		match scene_result {
			Err(e) => {
				println!("{:?}", e);
			}
			Ok(scene_result) => {
				*physics_scene = scene_result;
			}
		}
	}

	fn write_and_read_scene(&self, physics_scene: &mut PhysicsScene, renderer_scene: &mut RendererScene) {
		let tm = now();
		print!("Write and read scene starts at {}:{}:{}.{}. ", tm.tm_hour, tm.tm_min, tm.tm_sec, tm.tm_nsec/10_000_000);
		let _ = stdout().flush();
		let scenedata_filename = format!("/Users/christian/video/{:04}.txt", self.frame_number);
		let write_result = write_scene(&scenedata_filename, &physics_scene);
		if let Err(e) = write_result {
			println!("{:?}", e);
		}
		let scene_result = read_scene(&scenedata_filename);
		match scene_result {
			Err(e) => {
				println!("{:?}", e);
			}
			Ok(scene_result) => {
				*physics_scene = scene_result;
			}
		}
		*renderer_scene = RendererScene::new(physics_scene.clone());
		let duration = now() - tm;
		println!("It took {}:{}:{}.{}.", duration.num_hours(), duration.num_minutes()%60, duration.num_seconds()%60, (duration.num_milliseconds()%1000)/10);
	}

	fn render(&mut self, renderer_scene: &mut RendererScene) -> RendererOutput {
		let tm = now();
		println!("Render starts at {}:{}:{}.{}. ", tm.tm_hour, tm.tm_min, tm.tm_sec, tm.tm_nsec/10_000_000);
		let rows_rendered = Arc::new(Mutex::new(vec![false; self.height as usize]));
		let mut threads: Vec<JoinHandle<_>> = Vec::new();
		for _ in 0..self.number_of_threads {
			let renderer_scene_clone = renderer_scene.clone();
			let (width, height, spp_per_iteration, maximum_spp, maximum_error, maximum_brdf_value, perform_post_process, rows_rendered_clone) = (self.width, self.height, self.spp_per_iteration, self.maximum_spp, self.maximum_error, self.maximum_brdf_value, self.perform_post_process, rows_rendered.clone());
			threads.push(spawn(move || {
				let mut renderer = Renderer::new(width, spp_per_iteration, maximum_spp, maximum_error, maximum_brdf_value, perform_post_process, renderer_scene_clone);
				for y in 0..height {
					let mut row_rendered = {
						// @TODO: Get rid of unwrap().
						let mut rows_rendered = rows_rendered_clone.lock().unwrap();
						let row_rendered = rows_rendered[y as usize];
						rows_rendered[y as usize] = true;
						row_rendered
					};
					if !row_rendered {
						renderer.render(y);
					}
				}
				renderer.get_renderer_output_rows()
			}));
		}
		let mut renderer_output_rows: Vec<RendererOutputRow> = Vec::new();
		for i in 0..self.height {
			renderer_output_rows.push(RendererOutputRow::new(i, self.width));
		}
		for (i, thread) in threads.into_iter().enumerate() {
			println!("Waiting for thread {} to finish...", i);
			let result = thread.join();
			println!("Thread {} finished.", i);
			match result {
				Ok(mut renderer_output_rows_thread) => {
					for renderer_output_row in renderer_output_rows_thread {
						let row_number = renderer_output_row.row_number as usize;
						renderer_output_rows[row_number] = renderer_output_row;
					}
				}
				Err(e) => {
					println!("{:?}", e);
				}
			}
		}
		let duration = now() - tm;
		println!("It took {}:{}:{}.{}.", duration.num_hours(), duration.num_minutes()%60, duration.num_seconds()%60, (duration.num_milliseconds()%1000)/10);
		RendererOutput::new(self.width, self.height, renderer_output_rows)
	}

	fn write_to_image(&self, renderer_output: &mut RendererOutput, prefix: usize) {
		let tm = now();
		print!("Writing to starts at {}:{}:{}.{}. ", tm.tm_hour, tm.tm_min, tm.tm_sec, tm.tm_nsec/10_000_000);
		let _ = stdout().flush();
		let image_filename = format!("/Users/christian/video/{:01}_{:04}.png", prefix, self.frame_number);
		make_file(renderer_output, &image_filename);
		let duration = now() - tm;
		println!("It took {}:{}:{}.{}.", duration.num_hours(), duration.num_minutes()%60, duration.num_seconds()%60, (duration.num_milliseconds()%1000)/10);
	}

	fn post_process(&self, renderer_output: &mut RendererOutput, post_process_datum: (f64, i32)) {
		let tm = now();
		println!("Post-process starts at {}:{}:{}.{}. ", tm.tm_hour, tm.tm_min, tm.tm_sec, tm.tm_nsec/10_000_000);
		let (colors, number_of_rays) = RHF::new(self.width as i32, self.height as i32, post_process_datum.0, post_process_datum.1, self.search_window_radius, renderer_output).rhf();
		renderer_output.colors = colors;
		renderer_output.number_of_rays = number_of_rays;
		let duration = now() - tm;
		println!("It took {}:{}:{}.{}.", duration.num_hours(), duration.num_minutes()%60, duration.num_seconds()%60, (duration.num_milliseconds()%1000)/10);
	}

	fn simulate_physics(&self, physics_scene: &mut PhysicsScene) {
		let tm = now();
		print!("Simulating starts at {}:{}:{}.{}. ", tm.tm_hour, tm.tm_min, tm.tm_sec, tm.tm_nsec/10_000_000);
		let _ = stdout().flush();
		physics_scene.simulate_physics();
		let duration = now() - tm;
		println!("It took {}:{}:{}.{}.", duration.num_hours(), duration.num_minutes()%60, duration.num_seconds()%60, (duration.num_milliseconds()%1000)/10);
	}
}
