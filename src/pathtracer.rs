use std::io::{stdout, Write};
use std::thread::{JoinHandle, spawn};
use time::now;

use math::{add, add_u16};
use pngfile::{make_file, read_scene, write_scene};
use renderer::Renderer;
use rendereroutput::RendererOutput;
use rhf::RHF;
use sceneforphysics::SceneForPhysics;
use sceneforrendering::SceneForRendering;
use tracing::Tracing;

use NUMBER_OF_BINS;

pub struct Pathtracer {
	pub width: u32,
	pub height: u32,
	pub image_scale_factor: u32,
	pub frame_number: u32,
	pub number_of_threads_render: usize,
	pub number_of_threads_post_process: usize,
	pub continue_old_simulation: bool,
	pub post_process_data: Vec<(f32, i32)>, // max_distance and patch_size
	pub search_window_radius: i32,
	pub number_of_scales: i32,
	pub number_of_rays: u64,
	pub perform_post_process: bool,
	pub tracing: Tracing,
}

impl Pathtracer {
	pub fn new(width: u32, height: u32, image_scale_factor: u32, frame_number: u32, number_of_threads_render: usize, number_of_threads_post_process: usize, continue_old_simulation: bool, post_process_data: Vec<(f32, i32)>, search_window_radius: i32, number_of_scales: i32, number_of_rays: u64, perform_post_process: bool, tracing: Tracing) -> Self {
		Self {
			width,
			height,
			image_scale_factor,
			frame_number,
			number_of_threads_render,
			number_of_threads_post_process,
			continue_old_simulation,
			post_process_data,
			search_window_radius,
			number_of_scales,
			number_of_rays,
			perform_post_process,
			tracing,
		}
	}

	pub fn render_images_to_png(&mut self) {
		// Either start from an already created, stored scene, or create a new one.
		let mut scene_for_physics = SceneForPhysics::new();
		if self.continue_old_simulation {
			self.read_old_scene_for_physics(&mut scene_for_physics);
		} else {
			self.frame_number = 0;
		}
		let mut scene_for_rendering = SceneForRendering::new(scene_for_physics.clone());
		loop {
			println!("Frame {}.", self.frame_number);
			self.write_and_read_scene(&mut scene_for_physics, &mut scene_for_rendering);
			let mut renderer_output = self.render(&mut scene_for_rendering);
			self.write_to_image(&mut renderer_output, 0);
			if self.perform_post_process {
				for (i, post_process_datum) in self.post_process_data.iter().enumerate() {
					self.post_process(&mut renderer_output, *post_process_datum);
					self.write_to_image(&mut renderer_output, i+1);
				}
			}
			self.simulate_physics(&mut scene_for_physics);
			self.frame_number += 1;
		}
	}

	fn read_old_scene_for_physics(&self, scene_for_physics: &mut SceneForPhysics) {
		let scenedata_filename = format!("/Users/christian/video/{:04}.txt", self.frame_number);
		let scene_result = read_scene(&scenedata_filename);
		match scene_result {
			Err(e) => {
				println!("{:?}", e);
			}
			Ok(scene_result) => {
				*scene_for_physics = scene_result;
			}
		}
	}

	fn write_and_read_scene(&self, scene_for_physics: &mut SceneForPhysics, scene_for_rendering: &mut SceneForRendering) {
		let tm = now();
		print!("Write and read scene starts at {}:{}:{}.{}. ", tm.tm_hour, tm.tm_min, tm.tm_sec, tm.tm_nsec/10_000_000);
		let _ = stdout().flush();
		let scenedata_filename = format!("/Users/christian/video/{:04}.txt", self.frame_number);
		let write_result = write_scene(&scenedata_filename, &scene_for_physics);
		if let Err(e) = write_result {
			println!("{:?}", e);
		}
		let scene_result = read_scene(&scenedata_filename);
		match scene_result {
			Err(e) => {
				println!("{:?}", e);
			}
			Ok(scene_result) => {
				*scene_for_physics = scene_result;
			}
		}
		*scene_for_rendering = SceneForRendering::new(scene_for_physics.clone());
		let duration = now() - tm;
		println!("It took {}:{}:{}.{}.", duration.num_hours(), duration.num_minutes()%60, duration.num_seconds()%60, (duration.num_milliseconds()%1000)/10);
	}

	fn render(&self, scene_for_rendering: &mut SceneForRendering) -> RendererOutput {
		let tm = now();
		println!("Render starts at {}:{}:{}.{}. ", tm.tm_hour, tm.tm_min, tm.tm_sec, tm.tm_nsec/10_000_000);
		let mut threads: Vec<JoinHandle<_>> = Vec::new();
		let mut renderer_output = RendererOutput::new(self.width, self.height, self.image_scale_factor);
		for i in 0_usize..self.number_of_threads_render {
			let scene_for_rendering_clone = scene_for_rendering.clone();
			let (width, height, image_scale_factor, number_of_rays, tracing, write_percentage) = (self.width, self.height, self.image_scale_factor, self.number_of_rays, self.tracing, i == 0);
			threads.push(spawn(move || {
				Renderer::new(width, height, image_scale_factor, scene_for_rendering_clone).render(number_of_rays, tracing, write_percentage)
			}));
		}
		for (i, thread) in threads.into_iter().enumerate() {
			println!("Waiting for thread {} to finish...", i);
			let result = thread.join();
			println!("Thread {} finished.", i);
			match result {
				Ok(mut renderer_output_partial) => {
					for (j, pixel) in renderer_output_partial.pixels.iter().enumerate() {
						for k in 0..NUMBER_OF_BINS {
							renderer_output.pixels[j].bins[k] = add_u16(renderer_output.pixels[j].bins[k], pixel.bins[k]);
						}
						renderer_output.pixels[j].color = add(renderer_output.pixels[j].color, pixel.color);
					}
					for (j, color) in renderer_output_partial.colors.iter().enumerate() {
						renderer_output.colors[j] = add(renderer_output.colors[j], *color);
					}
					for (j, number_of_rays) in renderer_output_partial.number_of_rays.iter().enumerate() {
						renderer_output.number_of_rays[j] += number_of_rays;
					}
				}
				Err(e) => {
					println!("{:?}", e);
				}
			}
		}
		let duration = now() - tm;
		println!("It took {}:{}:{}.{}.", duration.num_hours(), duration.num_minutes()%60, duration.num_seconds()%60, (duration.num_milliseconds()%1000)/10);
		renderer_output
	}

	fn write_to_image(&self, renderer_output: &mut RendererOutput, prefix: usize) {
		let tm = now();
		print!("Writing starts at {}:{}:{}.{}. ", tm.tm_hour, tm.tm_min, tm.tm_sec, tm.tm_nsec/10_000_000);
		let _ = stdout().flush();
		let pixeldata_filename = "/Users/christian/video/pixeldata.txt".to_string();
		let image_filename = format!("/Users/christian/video/{:01}_{:04}.png", prefix, self.frame_number);
		make_file(self.width, self.height, self.image_scale_factor, renderer_output, &pixeldata_filename, &image_filename, self.number_of_rays*(self.number_of_threads_render as u64));
		let duration = now() - tm;
		println!("It took {}:{}:{}.{}.", duration.num_hours(), duration.num_minutes()%60, duration.num_seconds()%60, (duration.num_milliseconds()%1000)/10);
	}

	fn post_process(&self, mut renderer_output: &mut RendererOutput, post_process_datum: (f32, i32)) {
		let tm = now();
		println!("Post-process starts at {}:{}:{}.{}. ", tm.tm_hour, tm.tm_min, tm.tm_sec, tm.tm_nsec/10_000_000);
		//for max_distance_and_patch_radius in self.max_distances_and_patch_radii.iter() {
			let mut threads: Vec<JoinHandle<_>> = Vec::new();
			let rows_per_thread = ((self.height*self.image_scale_factor)/(self.number_of_threads_post_process as u32)) as i32;
			for i in 0_usize..self.number_of_threads_post_process {
				let mut number_of_rays = renderer_output.number_of_rays.clone();
				let colors = renderer_output.colors.clone();
				let bins = self.create_rhf_bins(&mut number_of_rays, &mut renderer_output);
				let (width, height, image_scale_factor) = (self.width as i32, self.height as i32, self.image_scale_factor as i32);
				let (max_distance, patch_radius, search_window_radius, number_of_scales) = (post_process_datum.0, post_process_datum.1, self.search_window_radius, self.number_of_scales);
				let (vertical_start, mut vertical_end) = (rows_per_thread*(i as i32), rows_per_thread*((i as i32)+1));
				if i == self.number_of_threads_post_process - 1 {
					vertical_end = (self.height as i32)*(self.image_scale_factor as i32);
				}
				println!("v_start = {:?}, v_end = {:?}", vertical_start, vertical_end);
				threads.push(spawn(move || {
					RHF::new(width, height, image_scale_factor, max_distance, patch_radius, search_window_radius, number_of_scales, vertical_start, vertical_end, number_of_rays, colors, bins).rhf()
					}));
				}
			renderer_output.colors = Vec::new();
			for thread in threads {
				let result = thread.join();
				match result {
					Ok(mut colors_denoised_partial) => {
						renderer_output.colors.append(&mut colors_denoised_partial);
					}
					Err(e) => {
						println!("{:?}", e);
					}
				}
			}
			//}
		let duration = now() - tm;
		println!("It took {}:{}:{}.{}.", duration.num_hours(), duration.num_minutes()%60, duration.num_seconds()%60, (duration.num_milliseconds()%1000)/10);
	}

	fn create_rhf_bins(&self, number_of_rays: &mut Vec<u32>, renderer_output: &mut RendererOutput) -> Vec<[[f32; 3]; NUMBER_OF_BINS]> {
		let mut bins: Vec<[[f32; 3]; NUMBER_OF_BINS]> = Vec::new();
		for i in 0..self.width*self.height*self.image_scale_factor*self.image_scale_factor {
			let mut bin = [[0.0; 3]; NUMBER_OF_BINS];
			for (j, color) in bin.iter_mut().enumerate() {
				for (k, color_component) in color.iter_mut().enumerate() {
					*color_component = 1.0/(number_of_rays[i as usize] as f32)*f32::from(renderer_output.pixels[i as usize].bins[j as usize][k]);
				}
			}
			bins.push(bin);
		}
		bins
	}

	fn simulate_physics(&self, scene_for_physics: &mut SceneForPhysics) {
		let tm = now();
		print!("Simulating starts at {}:{}:{}.{}. ", tm.tm_hour, tm.tm_min, tm.tm_sec, tm.tm_nsec/10_000_000);
		let _ = stdout().flush();
		scene_for_physics.simulate_physics();
		let duration = now() - tm;
		println!("It took {}:{}:{}.{}.", duration.num_hours(), duration.num_minutes()%60, duration.num_seconds()%60, (duration.num_milliseconds()%1000)/10);
	}
}
