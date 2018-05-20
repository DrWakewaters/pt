use std::f64;
use std::f64::consts::PI;

use rand::{OsRng, Rng};
use pcg_rand::Pcg32;

use hitpoint::Hitpoint;
use material::Material;
use math::{add, brdf, dot, elementwise_mul, mul, norm, normalised, random_from_brdf, sub};
use ray::Ray;
use renderershape::RendererShape;
use rendereroutput::RendererOutput;
use rendererscene::RendererScene;
use tracing::Tracing;
use tracing::Tracing::{Direct, Light, Importance, Bidirectional};

use NUMBER_OF_BINS;

//static MODULO: [usize; 5] = [0, 1, 2, 0, 1];

pub struct Renderer {
	pub width: u32,
	pub height: u32,
	pub number_of_rays: u64,
	pub perform_post_process: bool,
	pub write_percentage: bool,
	pub tracing: Tracing,
	pub scene: RendererScene,
}

impl Renderer {
	pub fn new(width: u32, height: u32, number_of_rays: u64, perform_post_process: bool, write_percentage: bool, tracing: Tracing, scene: RendererScene) -> Self {
		Self {
			width,
			height,
			number_of_rays,
			perform_post_process,
			write_percentage,
			tracing,
			scene,
		}
	}

	pub fn render(&mut self) -> RendererOutput {
		let mut pcg: Pcg32 = OsRng::new().unwrap().gen();
		let mut renderer_output = RendererOutput::new(self.width, self.height);
		match self.tracing {
			Direct => {
				self.trace(false, false, &mut pcg, &mut renderer_output);
			}
			Light => {
				//self.trace(false, true, &mut pcg, &mut renderer_output);
				self.light_tracing(&mut pcg, &mut renderer_output); // Old tracing
}
			Importance => {
				self.trace(true, false, &mut pcg, &mut renderer_output);
				//self.importance_tracing(&mut pcg, &mut renderer_output); // Old tracing
			}
			Bidirectional => {
				self.trace(true, true, &mut pcg, &mut renderer_output);
			}
		}
		renderer_output
	}

	fn trace(&mut self, trace_from_camera: bool, trace_from_light: bool, mut pcg: &mut Pcg32, renderer_output: &mut RendererOutput) {
		let number_of_lights = self.scene.lights.len();
		let random_exclusive_max_lights: u32 = <u32>::max_value() - <u32>::max_value()%(number_of_lights as u32);
		let number_of_cameras = self.scene.cameras.len();
		let random_exclusive_max_cameras: u32 = <u32>::max_value() - <u32>::max_value()%(number_of_cameras as u32);
		let spp = self.number_of_rays/((self.width*self.height) as u64);
		for x in -499..501 {
			if x%10 == 0 && self.write_percentage {
				println!("{:?} percent rendered", (x+500)/10);
			}
			for y in -499..501 {
				for _ in 0..spp {
					// Create a hitpoint path starting from the camera. The first hitpoint will be the hitpoint of the ray coming from the retina, through the pinhole of the camera - that is, it is a point somewhere on a renderer shape.
					let mut hitpoint_path_from_camera: Vec<Hitpoint> = Vec::new();
					let index = (pcg.gen_range(0, random_exclusive_max_lights)%(number_of_cameras as u32)) as usize;
					let r = pcg.next_f64()*0.5;
					let theta = pcg.next_f64()*2.0*PI;
					let dx = theta.cos()*r*0.0;
					let dy = theta.sin()*r*0.0;
					let mut point_on_retina = add(self.scene.cameras[index].retina_center, [x as f64, y as f64, 0.0]);
					let retina_position = 1000*(1000-(point_on_retina[1] as usize)) + (1000-(point_on_retina[0] as usize));
					point_on_retina = add(point_on_retina, [dx, dy, 0.0]);
					let direction = normalised(sub(self.scene.cameras[index].pinhole, point_on_retina));
					let mut ray = Ray::new(self.scene.cameras[index].pinhole, direction);
					let mut color = self.scene.cameras[index].color;
					self.create_hitpoint_path(&mut ray, &mut color, trace_from_camera, &mut hitpoint_path_from_camera, pcg);
					if hitpoint_path_from_camera.is_empty() {
						continue;
					}
					// Create a hitpoint path starting from a light.
					let mut hitpoint_path_from_light: Vec<Hitpoint> = Vec::new();
					let index = (pcg.gen_range(0, random_exclusive_max_cameras)%(number_of_lights as u32)) as usize;
					let mut direction = self.scene.lights[index].direction(&mut pcg);
					let mut ray = Ray::new(self.scene.lights[index].position, direction);
					let mut color = self.scene.lights[index].color;
					let light_hitpoint = Hitpoint::new(ray.position, [0.0, 0.0, 0.0], 0.0, [0.0, 0.0, 0.0], Material::new([1.0, 1.0, 1.0], 0.0, 0.0, 1.0, false), true, true);
					hitpoint_path_from_light.push(light_hitpoint);
					if trace_from_light {
						self.create_hitpoint_path(&mut ray, &mut color, true, &mut hitpoint_path_from_light, pcg);
					}
					if hitpoint_path_from_camera.is_empty() {
						continue;
					}
					// Connect the hitpaths.
					for hitpoint_in_light_path in &hitpoint_path_from_light {
						for hitpoint_in_camera_path in &hitpoint_path_from_camera {
							self.connect_and_store(hitpoint_in_light_path, hitpoint_in_camera_path, retina_position, renderer_output);
						}
					}
				}
			}
		}
	}

	fn create_hitpoint_path(&mut self, mut ray: &mut Ray, color: &mut [f64; 3], continue_after_first_hit: bool, hitpoint_path: &mut Vec<Hitpoint>, mut pcg: &mut Pcg32) {
		//let bullet_probability = 0.05;
		//let survival_boost_factor = 1.0/(1.0-bullet_probability);
		loop {
			/*
			let r = pcg.next_f64();
			if r < bullet_probability {
				return;
			} else {
				*color = mul(survival_boost_factor, *color);
			}
			*/
			// @TODO Store the accumulated color in the hitpoint.
			let hitpoint = self.closest_renderer_shape(&mut ray);
			if let Some(hitpoint) = hitpoint {
				let ingoing_direction = mul(-1.0, ray.direction);
				let (direction, _) = random_from_brdf(ingoing_direction, hitpoint.normal, hitpoint.material, pcg);
				ray.position = hitpoint.position;
				ray.direction = direction;
				*color = elementwise_mul(*color, hitpoint.material.color);
				hitpoint_path.push(hitpoint);
				return;
			} else {
				return;
			}
			if !continue_after_first_hit {
				return;
			}
		}
	}

	fn connect_and_store(&self, hitpoint_in_light_path: &Hitpoint, hitpoint_in_camera_path: &Hitpoint, retina_position: usize, renderer_output: &mut RendererOutput) {
		let direction = sub(hitpoint_in_camera_path.position, hitpoint_in_light_path.position);
		let distance = norm(direction);
		let direction_normalised = mul(1.0/distance, direction);
		if dot(hitpoint_in_light_path.normal, direction_normalised) < 0.0 || dot(hitpoint_in_camera_path.normal, direction_normalised) > 0.0 {
			return;
		}
		let closest_hitpoint = self.closest_renderer_shape(&mut Ray::new(hitpoint_in_light_path.position, direction_normalised));
		if let Some(closest_hitpoint) = closest_hitpoint {
			if distance-closest_hitpoint.distance > 1.0e-9 {
				return;
			}
		}
		let brdf_light = if hitpoint_in_light_path.is_end_point {
			1.0 // Should be lower? 1/pi ?
		} else {
			//println!("{}", brdf(mul(-1.0, hitpoint_in_light_path.incoming_direction), direction_normalised, hitpoint_in_light_path.normal, hitpoint_in_light_path.material));
			brdf(mul(-1.0, hitpoint_in_light_path.incoming_direction), direction_normalised, hitpoint_in_light_path.normal, hitpoint_in_light_path.material)
		};
		let brdf_camera = if hitpoint_in_camera_path.is_end_point {
			1.0 // Should be lower? 1/pi ?
		} else {
			brdf(mul(-1.0, hitpoint_in_camera_path.incoming_direction), mul(-1.0, direction_normalised), hitpoint_in_camera_path.normal, hitpoint_in_camera_path.material)
		};
		let color = mul(1_000_000.0*brdf_light*brdf_camera/(distance*distance), elementwise_mul(hitpoint_in_light_path.material.color, hitpoint_in_camera_path.material.color));
		self.store(retina_position, color, renderer_output);
	}

	fn store(&self, retina_position: usize, color: [f64; 3], renderer_output: &mut RendererOutput) {
		//let pos = y*(self.width as usize)+x;
		renderer_output.number_of_rays[retina_position] += 1;
		renderer_output.pixels[retina_position].color = add(renderer_output.pixels[retina_position].color, color);
		renderer_output.colors[retina_position] = add(renderer_output.colors[retina_position], color);
		/*
		if self.perform_post_process {
			renderer_output.pixels[pos].bins[(((color[0]).atan()*2.0/PI)*(NUMBER_OF_BINS as f64)) as usize][0] += 1;
			renderer_output.pixels[pos].bins[(((color[1]).atan()*2.0/PI)*(NUMBER_OF_BINS as f64)) as usize][1] += 1;
			renderer_output.pixels[pos].bins[(((color[2]).atan()*2.0/PI)*(NUMBER_OF_BINS as f64)) as usize][2] += 1;
		}
		*/
	}


	#[allow(dead_code)]
	fn importance_tracing(&mut self, mut pcg: &mut Pcg32, renderer_output: &mut RendererOutput) {
		let pinhole = [500.0, 500.0, -1000.0];
		let retina_center = [500.0, 500.0, -2000.0];
		let retina_normal = [0.0, 0.0, 1.0];
		/*
		for i in 0..self.number_of_rays {
			if i%(self.number_of_rays/100) == 0 && self.write_percentage {
				println!("{:?} percent rendered", (100.0*(i as f64)/(self.number_of_rays as f64)) as u32);
			}
			let direction = random_uniform_on_sphere(pcg);
			let distance = dot(sub(retina_center, pinhole), retina_normal)/dot(mul(-1.0, direction), retina_normal);
			if distance < 0.0 {
				continue;
			}
			let point_on_retina = add(mul(-1.0*distance, direction), pinhole);
			if (point_on_retina[2] - retina_center[2]).abs() > 1.0e-6 {
				println!("Error {}, {}", point_on_retina[2], retina_center[2]);
			}
			if point_on_retina[0] < 0.0 || point_on_retina[1] < 0.0 || point_on_retina[0] >= 1000.0 || point_on_retina[1] >= 1000.0 {
				continue;
			}
			let x_store = (1000.0-point_on_retina[0]) as usize;
			let y_store = (1000.0-point_on_retina[1]) as usize;
			let mut ray = Ray::new(pinhole, direction);
			let mut color = [1.0, 1.0, 1.0];
			self.importance_tracing_inner(&mut ray, &mut pcg, &mut color, renderer_output, x_store, y_store);
		}
		*/
		let spp = self.number_of_rays/((self.width*self.height) as u64);
		for x in -499..501 {
			if x%10 == 0 && self.write_percentage {
				println!("{:?} percent rendered", (x+500)/10);
			}
			for y in -499..501 {
				for _ in 0..spp {
					let r = pcg.next_f64()*0.5;
					let theta = pcg.next_f64()*2.0*PI;
					let dx = theta.cos()*r;
					let dy = theta.sin()*r;
					let mut point_on_retina = add(retina_center, [x as f64, y as f64, 0.0]);
					let x_store = 1000-(point_on_retina[0] as usize);
					let y_store = 1000-(point_on_retina[1] as usize);
					point_on_retina = add(point_on_retina, [dx, dy, 0.0]);
					let direction = normalised(sub(pinhole, point_on_retina));
					let mut ray = Ray::new(pinhole, direction);
					let mut color = [1.0, 1.0, 1.0];
					self.importance_tracing_inner(&mut ray, &mut pcg, &mut color, renderer_output, x_store, y_store);
				}
			}
		}
	}

	#[allow(dead_code)]
	fn light_tracing(&mut self, mut pcg: &mut Pcg32, mut renderer_output: &mut RendererOutput) {
		let number_of_lights = self.scene.lights.len();
		let random_exclusive_max: u32 = <u32>::max_value() - <u32>::max_value()%(number_of_lights as u32);
		for i in 0u64..self.number_of_rays {
			if i%(self.number_of_rays/100) == 0 && self.write_percentage {
				println!("{:?} percent rendered", (100.0*(i as f64)/(self.number_of_rays as f64)) as u32);
			}
			let light_index = (pcg.gen_range(0, random_exclusive_max)%(number_of_lights as u32)) as usize;
			let light_position = self.scene.lights[light_index].position;
			let light_direction = self.scene.lights[light_index].direction(pcg);
			let mut ray = Ray::new(light_position, light_direction);
			let mut color = [1.0, 1.0, 1.0];
			self.light_tracing_inner(&mut ray, &mut pcg, &mut color, &mut renderer_output);
		}
	}

	fn importance_tracing_inner(&mut self, mut ray: &mut Ray, mut pcg: &mut Pcg32, color: &mut [f64; 3], mut renderer_output: &mut RendererOutput, x_store: usize, y_store: usize) {
		let bullet_probability = 0.05;
		let survival_boost_factor = 1.0/(1.0-bullet_probability);
		loop {
			/*
			let r = pcg.next_f64();
			if r < bullet_probability {
				return;
			} else {
				*color = mul(survival_boost_factor, *color);
			}
			*/
			let hitpoint = self.closest_renderer_shape(&mut ray);
			if let Some(hitpoint) = hitpoint {
				let ingoing_direction = mul(-1.0, ray.direction);
				let (direction, _) = random_from_brdf(ingoing_direction, hitpoint.normal, hitpoint.material, pcg);
				ray.position = hitpoint.position;
				ray.direction = direction;
				*color = elementwise_mul(*color, hitpoint.material.color);
				self.direct_sample_lights(&mut ray, &mut pcg, ingoing_direction, *color, hitpoint.normal, hitpoint.material, &mut renderer_output, x_store, y_store);
			} else {
				return;
			}
		}
	}

	fn light_tracing_inner(&mut self, mut ray: &mut Ray, pcg: &mut Pcg32, color: &mut [f64; 3], mut renderer_output: &mut RendererOutput) {
		let bullet_probability = 0.05;
		let survival_boost_factor = 1.0/(1.0-bullet_probability);
		loop {
			/*
			let r = pcg.next_f64();
			if r < bullet_probability {
				return;
			} else {
				*color = mul(survival_boost_factor, *color);
			}*/
			let hitpoint = self.closest_renderer_shape(&mut ray);
			if let Some(hitpoint) = hitpoint {
				let ingoing_direction = mul(-1.0, ray.direction);
				let (direction, _) = random_from_brdf(ingoing_direction, hitpoint.normal, hitpoint.material, pcg);
				ray.position = hitpoint.position;
				ray.direction = direction;
				*color = elementwise_mul(*color, hitpoint.material.color);
				self.direct_sample_pinhole(&mut ray, ingoing_direction, *color, hitpoint.normal, hitpoint.material, &mut renderer_output);
				return;
			} else {
				return;
			}
		}
	}

	fn direct_sample_lights(&mut self, ray: &mut Ray, pcg: &mut Pcg32, ingoing_direction: [f64; 3], color: [f64; 3], normal: [f64; 3], material: Material, renderer_output: &mut RendererOutput, x_store: usize, y_store: usize) {
		let index = pcg.gen_range(0, self.scene.lights.len());
		let sample_position = self.scene.lights[index as usize].position;
		let sample_direction = sub(sample_position, ray.position);
		// @TODO Should this always be a condition, even with transmission?
		if dot(normal, sample_direction) < 0.0 {
			return;
		}
		let sample_distance = norm(sample_direction);
		let sample_direction_normalised = mul(1.0/sample_distance, sample_direction);
		let hitpoint = self.closest_renderer_shape(&mut Ray::new(ray.position, sample_direction_normalised));
		if let Some(hit) = hitpoint {
			if hit.distance < sample_distance {
				return;
			}
		}
		let brdf = brdf(ingoing_direction, sample_direction_normalised, normal, material);
		let color = mul(1_000_000.0*brdf/(sample_distance*sample_distance), color);
		self.store_old(x_store, y_store, color, renderer_output);
	}

	fn direct_sample_pinhole(&mut self, ray: &mut Ray, ingoing_direction: [f64; 3], color: [f64; 3], normal: [f64; 3], material: Material, renderer_output: &mut RendererOutput) {
		let sample_position = [500.0, 500.0, -1000.0];
		let sample_direction = sub(sample_position, ray.position);
		// @TODO Should this always be a condition, even with transmission?
		if dot(normal, sample_direction) < 0.0 {
			return;
		}
		let sample_distance = norm(sample_direction);
		let sample_direction_normalised = mul(1.0/sample_distance, sample_direction);
		let hitpoint = self.closest_renderer_shape(&mut Ray::new(ray.position, sample_direction_normalised));
		if let Some(hit) = hitpoint {
			if hit.distance < sample_distance {
				return;
			}
		}
		let brdf = brdf(ingoing_direction, sample_direction_normalised, normal, material);
		let retina_normal = [0.0, 0.0, 1.0];
		let retina_center = [500.0, 500.0, -2000.0];
		let distance_to_retina = dot(sub(retina_center, ray.position), retina_normal)/dot(sample_direction_normalised, retina_normal);
		if distance_to_retina < 0.0 {
			return;
		}
		let point_on_retina = add(mul(distance_to_retina, sample_direction_normalised), ray.position);
		if (point_on_retina[2] - retina_center[2]).abs() > 1.0e-6 {
			println!("Error {}, {}", point_on_retina[2], retina_center[2]);
		}
		if point_on_retina[0] < 0.0 || point_on_retina[1] < 0.0 || point_on_retina[0] >= 1000.0 || point_on_retina[1] >= 1000.0 {
			return;
		}
		let x_store = (1000.0-point_on_retina[0]) as usize;
		let y_store = (1000.0-point_on_retina[1]) as usize;
		let color = mul(1_000_000.0*brdf/(sample_distance*sample_distance), color);
		self.store_old(x_store, y_store, color, renderer_output);
	}

	fn store_old(&self, x_store: usize, y_store: usize, color: [f64; 3], renderer_output: &mut RendererOutput) {
		let pos = y_store*(self.width as usize)+x_store;
		renderer_output.number_of_rays[pos] += 1;
		renderer_output.pixels[pos].color = add(renderer_output.pixels[pos].color, color);
		renderer_output.colors[pos] = add(renderer_output.colors[pos], color);
		if self.perform_post_process {
			renderer_output.pixels[pos].bins[(((color[0]).atan()*2.0/PI)*(NUMBER_OF_BINS as f64)) as usize][0] += 1;
			renderer_output.pixels[pos].bins[(((color[1]).atan()*2.0/PI)*(NUMBER_OF_BINS as f64)) as usize][1] += 1;
			renderer_output.pixels[pos].bins[(((color[2]).atan()*2.0/PI)*(NUMBER_OF_BINS as f64)) as usize][2] += 1;

		}
	}

	// Find the closest hitpoint.
	fn closest_renderer_shape(&self, ray: &mut Ray) -> Option<Hitpoint>  {
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
			Some(Hitpoint::new(position, ray.direction, min_distance, normal, material, true, false))
		} else {
			None
		}
	}
}
