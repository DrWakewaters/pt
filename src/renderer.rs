use std::f64;

use rand::{OsRng, Rng};
use pcg_rand::Pcg32;

use math::{add, brdf, cross, dot, elementwise_mul, mul, norm, normalised, pick_reflection_from_brdf, sub};
use ray::Ray;
use rendereroutput::RendererOutput;
use sceneforrendering::SceneForRendering;
use tracing::Tracing;
use tracing::Tracing::{Bidirectional, Importance, Light};

use NUMBER_OF_BINS;

//static MODULO: [usize; 5] = [0, 1, 2, 0, 1];

pub struct Renderer {
	pub width: u32,
	pub height: u32,
	pub image_scale_factor: u32,
	scene: SceneForRendering,
}

impl Renderer {
	pub fn new(width: u32, height: u32, image_scale_factor: u32, scene: SceneForRendering) -> Self {
		Self {
			width,
			height,
			scene,
			image_scale_factor,
		}
	}

	pub fn render(&mut self, number_of_rays: u64, tracing: Tracing, write_percentage: bool) -> RendererOutput {
		let mut pcg: Pcg32 = OsRng::new().unwrap().gen();
		let mut renderer_output = RendererOutput::new(self.width, self.height, self.image_scale_factor);
		match tracing {
			Bidirectional => {
				//self.perform_work_bidirectional(&mut pcg, &mut renderer_output, number_of_rays, write_percentage);
			}
			Importance => {
				self.perform_work_from_eye(&mut pcg, &mut renderer_output, number_of_rays, write_percentage);
			}
			Light => {
				self.perform_work_sphere_lightsource(&mut pcg, &mut renderer_output, number_of_rays, write_percentage);
				//self.perform_work_triangle_lightsource(&mut pcg, &mut renderer_output, number_of_rays, write_percentage);
			}
		}
		renderer_output
	}

	// from eye: p -> sp1 -> sp2 -> sp3. Continue until the ray doesn't hit anything - no russian roulette. pinhole is always the same point. the retina is not part of the path, but the direction is chosen such that the same number of rays will have originate from each pixel on the retina.
	// from light: l -> sl1 -> sl2 -> sl3. Continue until the ray doesn't hit anything - no russian roulette. light is chosen randomly on the surface of the lightsource
	// connect and compute: say that we want to connect sp2 and sl2, that is make the path l -> sl1 -> sl2 => sp2 -> sp1 -> p,
	// then we just multiply the colors along the path, and multiply by a weight between sl2 and sp2. Computing this weight is a bit difficult but is done as follows:

	/*
	let f = brdf(incoming_direction, direction_to_retina, ray.normal_at_origin, ray.maximum_specular_angle_at_origin, ray.refractive_index_at_origin, ray.specular_probability_at_origin);
	let cos_retina = dot(direction_to_retina, retina_normal);
	let distance = norm(sub(ray.origin, pinhole));
	let pinhole_radius = 1.0;
	let alpha = (pinhole_radius/distance).atan();
	let proportion = 1.0 - alpha.cos();
	let color = mul(1_000_000.0*f*proportion*(1.0/(cos_retina*cos_retina*cos_retina)), ray.color);
	*/

	// l -> sl1 => p (light tracing with forcing to pinhole)
	// l => sp1 -> p (importance tracing with forcing to lightsource)
	// l => p (light or importance tracing - forcing directly between pinhole and lightsource)
	// How to trace if there is just one surface, s, between the pinhole and the lightsource? This surface can be a mirror!

	/*
	#[allow(dead_code)]
	fn perform_work_bidirectional(&mut self, mut pcg: &mut Pcg32, renderer_output: &mut RendererOutput, number_of_rays: u64, write_percentage: bool) {
		/*
		let ray_path_from_eye: Vec<> = Vec::new();
		let ray_path_from_light: Vec<> = Vec::new();
		for i in 0..number_of_rays {
			let ray_path_from_eye = build_ray_path_from_eye();
			let ray_path_from_light = build_ray_path_from_light();
			connect_ray_paths(ray_path_from_eye, ray_path_from_light);
		}
		*/
	}
	*/

	#[allow(dead_code)]
	fn perform_work_from_eye(&mut self, mut pcg: &mut Pcg32, renderer_output: &mut RendererOutput, number_of_rays: u64, write_percentage: bool) {
		let spp = number_of_rays/((self.width*self.image_scale_factor*self.height*self.image_scale_factor) as u64);
		let pinhole = [500.0, 500.0, -1000.0];
		let retina = [500.0, 500.0, 0.0];
		for x in 0..1000 {
			if x%10 == 0 && write_percentage {
				println!("{:?} percent rendered", x/10);
			}
			for y in 0..1000 {
				for _ in 0..spp {
					let rx = pcg.next_f64();
					let ry = pcg.next_f64();
					let position = add(retina, [(x as f64)+rx-500.0, (y as f64)+ry-500.0, 0.0]);
					let direction = normalised(sub(position, pinhole));
					let mut ray = Ray::new(position, direction, direction, [0.0, 0.0, 1.0], [0.0, 0.0, 0.0], 0.0, 1.0, 0.0, 1.0);
					let color = self.compute_from_eye(&mut ray, &mut pcg);
					let pos = (y*self.width*self.image_scale_factor+x) as usize;
					renderer_output.number_of_rays[pos] += 1;
					// Store information about the ray.
					/*
					renderer_output.pixels[pos].bins[(((color[0]*2.0).atan()*2.0/PI)*(NUMBER_OF_BINS as f64)) as usize][0] += 1;
					renderer_output.pixels[pos].bins[(((color[1]*2.0).atan()*2.0/PI)*(NUMBER_OF_BINS as f64)) as usize][1] += 1;
					renderer_output.pixels[pos].bins[(((color[2]*2.0).atan()*2.0/PI)*(NUMBER_OF_BINS as f64)) as usize][2] += 1;
					*/
					renderer_output.pixels[pos].color = add(renderer_output.pixels[pos].color, color);
					renderer_output.colors[pos] = add(renderer_output.colors[pos], color);
				}
			}
		}
	}

	fn compute_from_eye(&mut self, mut ray: &mut Ray, mut pcg: &mut Pcg32) -> [f64; 3] {
		let mut pixel_color: [f64; 3] = [1.0, 1.0, 1.0];
		let bullet_probability = 0.05;
		let survival_boost_factor = 1.0/(1.0-bullet_probability);
		loop {
			let r = pcg.next_f64();
			if r < bullet_probability {
				pixel_color = [0.0, 0.0, 0.0];
				break;
			}
			let (hit_object, is_lightsource, color) = self.compute_ray_object_intersection(&mut ray, &mut pcg, false);
			if !hit_object {
				pixel_color = [0.0, 0.0, 0.0];
				break;
			}
			pixel_color = mul(survival_boost_factor, elementwise_mul(pixel_color, color));
			if is_lightsource {
				break;
			}
		}
		pixel_color
	}

	#[allow(dead_code)]
	fn perform_work_sphere_lightsource(&mut self, mut pcg: &mut Pcg32, mut renderer_output: &mut RendererOutput, number_of_rays: u64, write_percentage: bool) {
		let lightsource_count = self.scene.sphere_lightsources.len();
		let random_exclusive_max: u32 = <u32>::max_value() - <u32>::max_value()%(lightsource_count as u32);
		for i in 0u64..number_of_rays {
			if i%(number_of_rays/100) == 0 && write_percentage {
				println!("{:?} percent rendered", (100.0*(i as f64)/(number_of_rays as f64)) as u32);
			}
			let lightsource_index = pcg.gen_range(0, random_exclusive_max)%(lightsource_count as u32);
			let mut ray = self.scene.sphere_lightsources[lightsource_index as usize].compute_ray_data(&mut pcg);
			self.compute(&mut ray, &mut pcg, &mut renderer_output);
		}
	}

	#[allow(dead_code)]
	fn perform_work_triangle_lightsource(&mut self, mut pcg: &mut Pcg32, mut renderer_output: &mut RendererOutput, number_of_rays: u64, write_percentage: bool) {
		let lightsource_count = self.scene.triangle_lightsources.len();
		let random_exclusive_max: u32 = <u32>::max_value() - <u32>::max_value()%(lightsource_count as u32);
		for i in 0u64..number_of_rays {
			if i%(number_of_rays/100) == 0 && write_percentage {
				println!("{:?} percent rendered", (100.0*(i as f64)/(number_of_rays as f64)) as u32);
			}
			let lightsource_index = pcg.gen_range(0, random_exclusive_max)%(lightsource_count as u32);
			let mut ray = self.scene.triangle_lightsources[lightsource_index as usize].compute_ray_data(&mut pcg);
			self.compute(&mut ray, &mut pcg, &mut renderer_output);
		}
	}

	fn compute(&mut self, mut ray: &mut Ray, mut pcg: &mut Pcg32, mut renderer_output: &mut RendererOutput) {
		loop {
			let bullet_probability = 0.05;
			let survival_boost_factor = 1.0/(1.0-bullet_probability);
			let r = pcg.next_f64();
			if r < bullet_probability {
				return;
			} else {
				ray.color = mul(survival_boost_factor, ray.color);
			}
			self.force_light_to_eye(&mut ray, &mut renderer_output);
			ray.old_direction = ray.direction;
			let (hit_object, _, color) = self.compute_ray_object_intersection(&mut ray, &mut pcg, false);
			if hit_object {
				ray.color = elementwise_mul(ray.color, color);
			} else {
				return;
			}
		}
	}

	fn force_light_to_eye(&mut self, ray: &mut Ray, renderer_output: &mut RendererOutput) {
		let pinhole = [500.0, 500.0, -1000.0];
		let direction_to_retina = normalised(sub(pinhole, ray.origin));
		if dot(ray.normal_at_origin, direction_to_retina) < 0.0 {
			return
		}
		let retina_normal = [0.0, 0.0, -1.0];
		let retina_center = [500.0, 500.0, -2000.0];
		let d = dot(sub(retina_center, pinhole), retina_normal)/dot(direction_to_retina, retina_normal);
		let retina_intersection = add(mul(d, direction_to_retina), pinhole);
		if (retina_intersection[2]-retina_center[2]).abs() > 1e-6 {
			return;
		}

		let dist = self.compute_ray_object_distance(ray.origin, direction_to_retina, true);
		let distance_to_pinhole = dot(sub(pinhole, ray.origin), retina_normal)/dot(direction_to_retina, retina_normal);
		if dist < distance_to_pinhole {
			return;
		}
		let incoming_direction = mul(-1.0, ray.old_direction);
		if d > 0.0 && retina_intersection[0] > 0.0 && retina_intersection[1] > 0.0 && retina_intersection[0] < f64::from(self.width) && retina_intersection[1] < f64::from(self.height) {
			let f = brdf(incoming_direction, direction_to_retina, ray.normal_at_origin, ray.maximum_specular_angle_at_origin, ray.refractive_index_at_origin, ray.specular_probability_at_origin);
			let cos_retina = dot(direction_to_retina, retina_normal);
			if cos_retina < 0.0 || cos_retina > 1.0 {
				println!("{:?}, {:?}", direction_to_retina, retina_normal);
				panic!();
			}
			let distance = norm(sub(ray.origin, pinhole));
			let pinhole_radius = 1.0;
			let alpha = (pinhole_radius/distance).atan();
			let proportion = 1.0 - alpha.cos();
			// Why does 1000000.0 make the pixel intensity have a maximum of seemingly exactly 1.0?
			let color = mul(1_000_000.0*f*proportion*(1.0/(cos_retina*cos_retina*cos_retina)), ray.color);
			self.store_intersection(retina_intersection, color, renderer_output);
		}
	}

	fn store_intersection(&self, retina_intersection: [f64; 3], color: [f64; 3], renderer_output: &mut RendererOutput) {
		let y = ((f64::from(self.height) - retina_intersection[1])*f64::from(self.image_scale_factor)) as u32;
		let x = ((f64::from(self.width) - retina_intersection[0])*f64::from(self.image_scale_factor)) as u32;
		if (y*self.width+x) < self.width*self.height*self.image_scale_factor*self.image_scale_factor {
			let pos = (y*self.width*self.image_scale_factor+x) as usize;
			renderer_output.number_of_rays[pos] += 1;
			// Store information about the ray.
			/*
			renderer_output.pixels[pos].bins[(((color[0]*2.0).atan()*2.0/PI)*(NUMBER_OF_BINS as f64)) as usize][0] += 1;
			renderer_output.pixels[pos].bins[(((color[1]*2.0).atan()*2.0/PI)*(NUMBER_OF_BINS as f64)) as usize][1] += 1;
			renderer_output.pixels[pos].bins[(((color[2]*2.0).atan()*2.0/PI)*(NUMBER_OF_BINS as f64)) as usize][2] += 1;
			*/
			renderer_output.pixels[pos].color = add(renderer_output.pixels[pos].color, color);
			renderer_output.colors[pos] = add(renderer_output.colors[pos], color);
		} else {
			println!("Can't write to pixels[{}] (x = {}, y = {}), color = {:?}", y*self.width*self.image_scale_factor+x, retina_intersection[0], retina_intersection[1], color);
		}
	}

	// Find the closest triangle or sphere intersection of a given ray.
	fn compute_ray_object_intersection(&self, ray: &mut Ray, pcg: &mut Pcg32, camera_ray: bool) -> (bool, bool, [f64; 3]) {
		let mut color = [0.0, 0.0, 0.0];
		let mut is_lightsource = false;
		//println!("position = {:?}, direction = {:?}", ray.origin, ray.direction);
		let (triangle_hit, triangle_index, triangle_distance, triangle_is_lightsource, triangle_refractive_index, triangle_is_opaque) = self.compute_ray_triangle_intersection(ray.origin, ray.direction, camera_ray);
		let (sphere_hit, sphere_index, sphere_distance, sphere_is_lightsource, sphere_refractive_index, sphere_is_opaque) = self.compute_ray_sphere_intersection(ray.origin, ray.direction, camera_ray);
		// The ray hit nothing.
		if !triangle_hit && !sphere_hit {
			return (false, is_lightsource, color);
		}
		let incoming_direction = mul(-1.0, ray.direction);
		// What did the ray hit? A triangle or a sphere?
		if triangle_hit && (!sphere_hit || triangle_distance < sphere_distance) {
			ray.origin = add(ray.origin, mul(triangle_distance, ray.direction));
			ray.normal_at_origin = self.scene.triangle_surfaces[triangle_index].normal;
			ray.specular_probability_at_origin = self.scene.triangle_surfaces[triangle_index].specular_probability;
			ray.maximum_specular_angle_at_origin = self.scene.triangle_surfaces[triangle_index].maximum_specular_angle;
			let (output_direction, output_refractive_index) = pick_reflection_from_brdf(incoming_direction, ray.normal_at_origin, ray.maximum_specular_angle_at_origin, ray.refractive_index_at_origin, triangle_refractive_index, ray.specular_probability_at_origin, triangle_is_opaque, pcg);
			ray.direction = output_direction;
			ray.refractive_index_at_origin_in_direction = output_refractive_index;
			is_lightsource = triangle_is_lightsource;
			color = self.scene.triangle_surfaces[triangle_index].compute_intersection_color();
		} else {
			ray.origin = add(ray.origin, mul(sphere_distance, ray.direction));
			ray.normal_at_origin = normalised(sub(ray.origin, self.scene.sphere_surfaces[sphere_index].center));
			if dot(ray.normal_at_origin, ray.direction) < 0.0 {
				// If we hit a sphere boundary from the inside.
				ray.normal_at_origin = mul(-1.0, ray.normal_at_origin);
			}
			ray.specular_probability_at_origin = self.scene.sphere_surfaces[sphere_index].specular_probability;
			ray.maximum_specular_angle_at_origin = self.scene.sphere_surfaces[sphere_index].maximum_specular_angle;
			let (output_direction, output_refractive_index) = pick_reflection_from_brdf(incoming_direction, ray.normal_at_origin, ray.maximum_specular_angle_at_origin, ray.refractive_index_at_origin, sphere_refractive_index, ray.specular_probability_at_origin, sphere_is_opaque, pcg);
			ray.direction = output_direction;
			ray.refractive_index_at_origin_in_direction = output_refractive_index;
			is_lightsource = sphere_is_lightsource;
			color = self.scene.sphere_surfaces[sphere_index].compute_intersection_color();
		}
		(true, is_lightsource, color)
	}

	// Find the distance to the closest sphere or triangle of a given ray.
	fn compute_ray_object_distance(&self, intersection: [f64; 3], direction: [f64; 3], camera_ray: bool) -> f64 {
		let (_, _, triangle_distance, _, _, _) = self.compute_ray_triangle_intersection(intersection, direction, camera_ray);
		let (_, _, sphere_distance, _, _, _) = self.compute_ray_sphere_intersection(intersection, direction, camera_ray);
		if triangle_distance < sphere_distance {
			triangle_distance
		} else {
			sphere_distance
		}
	}

	// Find the closest triangle intersection of a given ray.
	// See http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.189.5084&rep=rep1&type=pdf.
	fn compute_ray_triangle_intersection(&self, intersection: [f64; 3], direction: [f64; 3], camera_ray: bool) -> (bool, usize, f64, bool, f64, bool) {
		let mut triangle_hit = false;
		let mut triangle_index = 0;
		let mut triangle_distance = f64::MAX;
		let mut triangle_is_lightsource = false;
		let mut triangle_refractive_index = 1.0;
		let mut triangle_is_opaque = true;
		// Compute information about the closest hit.
		for (i, t) in self.scene.triangle_surfaces.iter().enumerate() {
			if camera_ray && t.invisible_for_camera_ray {
				continue;
			}
			/*
			let k = t.k as usize;
			let u = MODULO[(t.k+1) as usize];
			let v = MODULO[(t.k+2) as usize];
			let nd = 1.0/(direction[k]+t.n_u*direction[u]+t.n_v*direction[v]);
			let d = (t.n_d-intersection[k]-t.n_u*intersection[u]-t.n_v*intersection[v])*nd;
			if d < 1e-9 || d > triangle_distance {
				continue;
			}
			let hu = intersection[u]+d*direction[u];
			let hv = intersection[v]+d*direction[v];
			let lambda = hu*t.b_nu+hv*t.b_nv+t.b_d;
			if lambda < 0.0 {
				continue;
			}
			let mue = hu*t.c_nu+hv*t.c_nv+t.c_d;
			if mue < 0.0 || lambda+mue > 1.0 {
				continue;
			}*/
			let h = cross(direction, t.e2);
			let a = dot(t.e1, h);
			if a.abs() < 1.0e-6 {
				continue;
			}
			let f = 1.0/a;
			let s = sub(intersection, t.node0);
			let u = f*dot(s, h);
			if u < 0.0 || u > 1.0 {
				continue;
			}
			let q = cross(s, t.e1);
			let v = f*dot(direction, q);
			if v < 0.0 || u+v > 1.0 {
				continue;
			}
			let d = f*dot(t.e2, q);
			if d < 1.0e-6 {
				continue;
			}
			if d > triangle_distance {
				continue;
			}
			triangle_hit = true;
			triangle_index = i;
			triangle_distance = d;
			triangle_is_lightsource = t.is_lightsource;
			triangle_refractive_index = t.refractive_index;
			triangle_is_opaque = t.is_opaque;
		}
		(triangle_hit, triangle_index, triangle_distance, triangle_is_lightsource, triangle_refractive_index, triangle_is_opaque)
	}

	// Find the closest sphere intersection of a given ray.
	fn compute_ray_sphere_intersection(&self, intersection: [f64; 3], direction: [f64; 3], camera_ray: bool) -> (bool, usize, f64, bool, f64, bool) {
		let mut sphere_hit = false;
		let mut sphere_index = 0;
		let mut sphere_distance = f64::MAX;
		let mut sphere_is_lightsource = false;
		let mut sphere_refractive_index = 1.0;
		let mut sphere_is_opaque = true;
		// Compute information about the closest hit.
		for (i, s) in self.scene.sphere_surfaces.iter().enumerate() {
			if camera_ray && s.invisible_for_camera_ray {
				continue;
			}
			let b = sub(intersection, s.center);
			let a = dot(b, direction)*dot(b, direction) - dot(b, b) + s.radius*s.radius;
			if a < 0.0 {
				continue;
			}
			let d1 = -1.0*dot(sub(intersection, s.center), direction) + a.sqrt();
			let d2 = -1.0*dot(sub(intersection, s.center), direction) - a.sqrt();
			let mut d = f64::MAX;
			if d2 > 1e-6 {
				d = d2;
			} else if d1 > 1e-6 {
				d = d1;
			}
			if d < sphere_distance {
				sphere_hit = true;
				sphere_index = i;
				sphere_distance = d;
				sphere_is_lightsource = s.is_lightsource;
				sphere_refractive_index = s.refractive_index;
				sphere_is_opaque = s.is_opaque
			}
		}
		(sphere_hit, sphere_index, sphere_distance, sphere_is_lightsource, sphere_refractive_index, sphere_is_opaque)
	}
}
