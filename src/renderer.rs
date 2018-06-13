use std::f64;

use rand::{OsRng, Rng};
use pcg_rand::Pcg32;

use hitpoint::Hitpoint;
use math::{add, brdf, dot, elementwise_mul, intensity_to_color, min, mul, norm, random_from_brdf, sub};
use ray::Ray;
use renderershape::RendererShape;
use rendereroutputpixel::RendererOutputPixel;
use rendererscene::RendererScene;

use NUMBER_OF_BINS;

//static MODULO: [usize; 5] = [0, 1, 2, 0, 1];

pub struct Renderer {
	width: u32,
	height: u32,
	spp_per_iteration: u32,
	maximum_spp: u32,
	maximum_error: f64,
	maximum_brdf_value: f64,
	perform_post_process: bool,
	scene: RendererScene,
	renderer_output_pixels: Vec<RendererOutputPixel>
}

impl Renderer {
	pub fn new(width: u32, height: u32, spp_per_iteration: u32, maximum_spp: u32, maximum_error: f64, maximum_brdf_value: f64, perform_post_process: bool, scene: RendererScene) -> Self {
		let renderer_output_pixels: Vec<RendererOutputPixel> = Vec::new();
		Self {
			width,
			height,
			spp_per_iteration,
			maximum_spp,
			maximum_error,
			maximum_brdf_value,
			perform_post_process,
			scene,
			renderer_output_pixels,
		}
	}

	pub fn get_renderer_output_pixels(&self) -> Vec<RendererOutputPixel> {
		self.renderer_output_pixels.to_vec()
	}

	pub fn render(&mut self, y: u32, x: u32) {
		if x == 0 {
			println!("rendering, y={}", y);
		}
		let number_of_light_spheres = self.scene.light_spheres.len();
		let random_exclusive_max_lights: u32 = <u32>::max_value() - <u32>::max_value()%(number_of_light_spheres as u32);
		let number_of_cameras = self.scene.cameras.len();
		self.renderer_output_pixels.push(RendererOutputPixel::new(y, x));
		let last_pos = self.renderer_output_pixels.len()-1;
		let mut pcg: Pcg32 = OsRng::new().unwrap().gen();
		let mut colors: Vec<[f64; 3]> = Vec::new();
		let mut converged = false;
		// Loop over this pixel until we estimate the error to be small enough.
		let mut iterations = 0;
		while !converged {
			iterations += 1;
			for _ in 0..self.spp_per_iteration {
				// Create a hitpoint path starting from the camera. The first hitpoint will be the hitpoint of the ray coming from the retina, through the pinhole of the camera - that is, it is a point somewhere on a renderer shape.
				let mut hitpoint_path_from_camera: Vec<Hitpoint> = Vec::new();
				self.renderer_output_pixels[last_pos].number_of_rays += 1.0;
				let index = (pcg.gen_range(0, random_exclusive_max_lights)%(number_of_cameras as u32)) as usize;
				let mut ray = self.scene.cameras[index].create_ray(x, y, self.width, self.height, &mut pcg);
				let mut color = self.scene.cameras[index].color;
				self.create_hitpoint_path(&mut ray, &mut color, &mut hitpoint_path_from_camera, &mut pcg);
				let mut total_color = [0.0, 0.0, 0.0];
				if hitpoint_path_from_camera.is_empty() {
					colors.push(total_color);
					continue;
				}
				// Connect the camera path hitpoints with points on lights.
				for hitpoint_in_camera_path in &hitpoint_path_from_camera {
					let color = self.connect_and_compute_color(hitpoint_in_camera_path, &mut pcg);
					total_color = add(total_color, color);
				}
				self.store(last_pos, total_color);
				colors.push(total_color);
			}
			if colors.is_empty() {
				break;
			}
			// Estimate the error. If it's too large, create more rays for this pixel.
			let number_of_batches: usize = 20;
			let number_of_rays = colors.len();
			let batch_size = number_of_rays/number_of_batches;
			let mut averages = self.averages(&colors, number_of_batches, batch_size);
			self.gamma_correct_averages(&mut averages);
			let use_standard_deviation = true;
			let error = if use_standard_deviation {
				self.standard_deviation(&averages)
			} else {
				self.maximum_distance(&averages)
			};
			if iterations%100 == 0 {
				println!("Currently: {} iterations at ({}, {}) with an error of {}. Iterating until it is less than {}).", iterations, x, y, error, self.maximum_error);
			}
			if error < self.maximum_error || iterations*self.spp_per_iteration >= self.maximum_spp {
				converged = true;
			}
		}
	}

	// @TODO: Fix so that it works even if colors.len()%batch_size != 0.
	fn averages(&self, colors: &[[f64; 3]], number_of_batches: usize, batch_size: usize) -> Vec<[f64; 3]> {
		let mut averages: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0]; number_of_batches];
		for i in 0..colors.len() {
			averages[i/batch_size] = add(averages[i/batch_size], colors[i]);
		}
		for average in &mut averages {
			*average = mul(1.0/(batch_size as f64), *average);
		}
		averages
	}

	fn gamma_correct_averages(&self, averages: &mut Vec<[f64; 3]>) {
		for average in averages {
			*average = intensity_to_color(*average);
		}
	}

	fn standard_deviation(&self, averages: &[[f64; 3]]) -> f64 {
		let length = averages.len();
		let mut r = 0.0;
		let mut g = 0.0;
		let mut b = 0.0;
		let mut r_squared = 0.0;
		let mut g_squared = 0.0;
		let mut b_squared = 0.0;
		for average in averages {
			r += average[0];
			g += average[1];
			b += average[2];
			r_squared += average[0]*average[0];
			g_squared += average[1]*average[1];
			b_squared += average[2]*average[2];
		}
		((r_squared + g_squared + b_squared - (r*r+g*g+b*b)/(length as f64))/(length as f64-1.0)).sqrt()
	}

	fn maximum_distance(&self, averages: &[[f64; 3]]) -> f64 {
		let mut smallest = [f64::MAX, f64::MAX, f64::MAX];
		let mut largest = [f64::MIN, f64::MIN, f64::MIN];
		for average in averages {
			for j in 0..3 {
				if average[j] < smallest[j] {
					smallest[j] = average[j];
				}
				if average[j] > largest[j] {
					largest[j] = average[j];
				}
			}
		}
		let max_r_distance = (largest[0]-smallest[0]).abs();
		let max_g_distance = (largest[1]-smallest[1]).abs();
		let max_b_distance = (largest[2]-smallest[2]).abs();
		max_r_distance+max_g_distance+max_b_distance
	}

	fn create_hitpoint_path(&mut self, mut ray: &mut Ray, color: &mut [f64; 3], hitpoint_path: &mut Vec<Hitpoint>, pcg: &mut Pcg32) {
		let bullet_probability = 0.0;
		let survival_boost_factor = 1.0/(1.0-bullet_probability);
		loop {
			let r = pcg.next_f64();
			if r < bullet_probability {
				return;
			} else {
				*color = mul(survival_boost_factor, *color);
			}
			let distance_from_retina = if hitpoint_path.is_empty() {
				0.0
			} else {
				hitpoint_path[hitpoint_path.len()-1].distance_from_retina
			};
			let hitpoint = self.closest_renderer_shape(&mut ray, distance_from_retina);
			if let Some(mut hitpoint) = hitpoint {
				let ingoing_direction = mul(-1.0, ray.direction);
				let (refractive_index_1, refractive_index_2) = if hitpoint.hit_from_outside {
					(1.0, hitpoint.material.refractive_index)
				} else {
					(hitpoint.material.refractive_index, 1.0)
				};
				let normal = if dot(ingoing_direction, hitpoint.normal) > 0.0 {
					hitpoint.normal
				} else {
					mul(-1.0, hitpoint.normal)
				};
				let (direction, _) = random_from_brdf(ingoing_direction, normal, hitpoint.material, refractive_index_1, refractive_index_2, pcg);
				ray.position = hitpoint.position;
				ray.direction = direction;
				*color = elementwise_mul(*color, hitpoint.material.color);
				hitpoint.accumulated_color = *color;
				hitpoint_path.push(hitpoint);
			} else {
				return;
			}
		}
	}

	fn connect_and_compute_color(&self, hitpoint_in_camera_path: &Hitpoint, mut pcg: &mut Pcg32) -> [f64; 3] {
		let number_of_light_spheres = self.scene.light_spheres.len();
		let number_of_cameras = self.scene.cameras.len();
		let random_exclusive_max_cameras: u32 = <u32>::max_value() - <u32>::max_value()%(number_of_cameras as u32);
		let index = (pcg.gen_range(0, random_exclusive_max_cameras)%(number_of_light_spheres as u32)) as usize;
		let light_position = self.scene.light_spheres[index].get_position(&mut pcg);
		let light_color = self.scene.light_spheres[index].color;
		let direction = sub(hitpoint_in_camera_path.position, light_position);
		let distance = norm(direction);
		let direction_normalised = mul(1.0/distance, direction);
		if dot(hitpoint_in_camera_path.normal, direction_normalised) > 0.0 {
			return [0.0, 0.0, 0.0];
		}
		let closest_hitpoint = self.closest_renderer_shape(&mut Ray::new(light_position, direction_normalised), hitpoint_in_camera_path.distance_from_retina);
		if let Some(closest_hitpoint) = closest_hitpoint {
			if distance-closest_hitpoint.distance > 1.0e-9 {
				return [0.0, 0.0, 0.0];
			}
		}
		let ingoing_direction = mul(-1.0, hitpoint_in_camera_path.incoming_direction);
		let outgoing_direction = mul(-1.0, direction_normalised);
		let normal = if dot(ingoing_direction, hitpoint_in_camera_path.normal) > 0.0 {
			hitpoint_in_camera_path.normal
		} else {
			mul(-1.0, hitpoint_in_camera_path.normal)
		};
		let (refractive_index_1, refractive_index_2) = if hitpoint_in_camera_path.hit_from_outside {
			(1.0, hitpoint_in_camera_path.material.refractive_index)
		} else {
			(hitpoint_in_camera_path.material.refractive_index, 1.0)
		};
		// @TODO: Get rid of the upper limit of the brdf.
		let brdf = min(brdf(ingoing_direction, outgoing_direction, normal, refractive_index_1, refractive_index_2, hitpoint_in_camera_path.material), self.maximum_brdf_value);
		//color = add(color, mul(brdf/(distance*distance*f64::from(number_of_light_connections)), elementwise_mul(hitpoint_in_camera_path.accumulated_color, light_color)));
		let distance_from_retina = hitpoint_in_camera_path.distance_from_retina + distance;
		mul(brdf/(distance_from_retina*distance_from_retina), elementwise_mul(hitpoint_in_camera_path.accumulated_color, light_color))
	}

	fn store(&mut self, last_pos: usize, color: [f64; 3]) {
		self.renderer_output_pixels[last_pos].pixel.color = add(self.renderer_output_pixels[last_pos].pixel.color, color);
		self.renderer_output_pixels[last_pos].color = add(self.renderer_output_pixels[last_pos].color, color);
		self.renderer_output_pixels[last_pos].number_of_bin_elements += 1;
		if self.perform_post_process {
			let color = intensity_to_color(color);
			self.renderer_output_pixels[last_pos].pixel.bins[(color[0]/255.0*(NUMBER_OF_BINS as f64)) as usize][0] += 1;
			self.renderer_output_pixels[last_pos].pixel.bins[(color[1]/255.0*(NUMBER_OF_BINS as f64)) as usize][1] += 1;
			self.renderer_output_pixels[last_pos].pixel.bins[(color[2]/255.0*(NUMBER_OF_BINS as f64)) as usize][2] += 1;
		}
	}

	// Find the closest hitpoint.
	fn closest_renderer_shape(&self, ray: &mut Ray, distance_from_retina: f64) -> Option<Hitpoint>  {
		let mut min_distance = f64::MAX;
		let mut closest_renderer_shape_index: Option<usize> = None;
		let mut closest_is_a_sphere = true;
		for (i, sphere) in self.scene.renderer_spheres.iter().enumerate() {
			let distance = sphere.distance(&ray);
			if distance < min_distance {
				min_distance = distance;
				closest_renderer_shape_index = Some(i);
			}
		}
		for (i, triangle) in self.scene.renderer_triangles.iter().enumerate() {
			let distance = triangle.distance(&ray);
			if distance < min_distance {
				min_distance = distance;
				closest_renderer_shape_index = Some(i);
				closest_is_a_sphere = false;
			}
		}
		if let Some(index) = closest_renderer_shape_index {
			let position = add(ray.position, mul(min_distance, ray.direction));
			let normal = if closest_is_a_sphere {
				self.scene.renderer_spheres[index].normal(position)
			} else {
				self.scene.renderer_triangles[index].normal(position)
			};
			let material = if closest_is_a_sphere {
				self.scene.renderer_spheres[index].material()
			} else {
				self.scene.renderer_triangles[index].material()
			};
			let hit_from_outside = dot(ray.direction, normal) < 0.0;
			Some(Hitpoint::new(position, ray.direction, min_distance, normal, material, hit_from_outside, false, [0.0, 0.0, 0.0], distance_from_retina+min_distance))
		} else {
			None
		}
	}
}
