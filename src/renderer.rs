use std::f64;
use std::f64::consts::PI;

use rand::{OsRng, Rng};
use pcg_rand::Pcg32;

use math::{add, brdf, dot, elementwise_mul, mul, norm, normalised, pick_reflection_from_brdf, pick_reflection_lambertian, pick_reflection_uniform, point_in_triangle, sub};
use scene::SceneForRendering;
use NUMBER_OF_BINS;

static MODULO: [usize; 5] = [0, 1, 2, 0, 1];

pub struct Renderer {
	pub width: u32,
	pub height: u32,
	pub image_scale_factor: u32,
	scene: SceneForRendering,
}

#[derive(Clone)]
pub struct RendererOutput {
	pub number_of_rays: Vec<u32>,
	pub pixels: Vec<Pixel>,
	pub colors: Vec<[f64; 3]>,
}

#[derive(Clone)]
pub struct Pixel {
	pub bins: [[u16; 3]; NUMBER_OF_BINS],
	pub color: [f64; 3],
}

impl Pixel {
	pub fn new() -> Self {
		Self {
			bins: [[0; 3]; NUMBER_OF_BINS],
			color: [0.0; 3],
		}
	}
}

impl RendererOutput {
	pub fn new(width: u32, height: u32, image_scale_factor: u32) -> Self {
		let mut number_of_rays: Vec<u32> = Vec::new();
		let mut pixels: Vec<Pixel> = Vec::new();
		let mut colors: Vec<[f64; 3]> = Vec::new();
		for _ in 0..width*height*image_scale_factor*image_scale_factor {
			number_of_rays.push(0);
			pixels.push(Pixel::new());
			colors.push([0.0, 0.0, 0.0]);
		}
		Self {
			number_of_rays,
			pixels,
			colors,
		}
	}
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

	pub fn render(&mut self, number_of_rays: u64, write_percentage: bool) -> RendererOutput {
		let mut pcg: Pcg32 = OsRng::new().unwrap().gen();
		let mut renderer_output = RendererOutput::new(self.width, self.height, self.image_scale_factor);
		self.perform_work_sphere_lightsource(&mut pcg, &mut renderer_output, number_of_rays, write_percentage);
		renderer_output
	}

	#[allow(dead_code)]
	fn perform_work_sphere_lightsource(&mut self, mut pcg: &mut Pcg32, mut renderer_output: &mut RendererOutput, number_of_rays: u64, write_percentage: bool) {
		let lightsource_count = self.scene.sphere_lightsources.len();
		let random_exclusive_max: u32 = <u32>::max_value() - <u32>::max_value()%(lightsource_count as u32);
		for i in 0u64..number_of_rays {
			if i%(number_of_rays/10) == 0 && write_percentage {
				println!("{:?} percent rendered", (100.0*(i as f64)/(number_of_rays as f64)) as u32);
			}
			let lightsource_index = pcg.gen_range(0, random_exclusive_max)%(lightsource_count as u32);
			let (point, direction, normal, mut pixel) = self.scene.sphere_lightsources[lightsource_index as usize].compute_ray_data(&mut pcg);
			self.compute(point, direction, normal, &mut pixel, &mut pcg, &mut renderer_output);
		}
	}

	#[allow(dead_code)]
	fn perform_work_triangle_lightsource(&mut self, mut pcg: &mut Pcg32, mut renderer_output: &mut RendererOutput, number_of_rays: u64, write_percentage: bool) {
		let lightsource_count = self.scene.triangle_lightsources.len();
		let random_exclusive_max: u32 = <u32>::max_value() - <u32>::max_value()%(lightsource_count as u32);
		for i in 0u64..number_of_rays {
			if i%(number_of_rays/10) == 0 && write_percentage {
				println!("{:?} percent rendered", (100.0*(i as f64)/(number_of_rays as f64)) as u32);
			}
			let lightsource_index = pcg.gen_range(0, random_exclusive_max)%(lightsource_count as u32);
			let lightsource = self.scene.triangle_lightsources[lightsource_index as usize];
			let mut point_found = false;
			while !point_found {
				let r1_u32: u32 = pcg.gen();
				let r2_u32: u32 = pcg.gen();
				let point = add(add(mul((r1_u32 as f64)/(<u32>::max_value()) as f64, lightsource.e1), mul((r2_u32 as f64)/(<u32>::max_value() as f64), lightsource.e2)), lightsource.node0);
				if point_in_triangle(point, lightsource.node0, lightsource.node1, lightsource.node2) {
					point_found = true;
					let direction = pick_reflection_lambertian(lightsource.normal, &mut pcg);
					let mut pixel = lightsource.color;
					let normal = lightsource.normal;
					self.compute(point, direction, normal, &mut pixel, &mut pcg, &mut renderer_output);
				}
			}
		}
	}

	fn compute(&mut self, intersection: [f64; 3], direction: [f64; 3], normal: [f64; 3], mut pixel: &mut[f64; 3], mut pcg: &mut Pcg32, mut renderer_output: &mut RendererOutput) {
		let mut intersection = intersection;
		let mut direction = direction;
		let mut normal = normal;
		let mut hit_object = true;
		let mut color = [0.0, 0.0, 0.0];
		let mut old_direction = direction;
		let mut maximum_specular_angle = 0.0;
		let mut refractive_index = 1.0;
		let mut specular_probability = 0.0;
		loop {
			if hit_object {
				let bullet_probability = 0.05;
				let survival_boost_factor = 1.0/(1.0-bullet_probability);
				let r_u32: u32 = pcg.gen();
				let r = (r_u32 as f64)/(<u32>::max_value() as f64);
				if r < bullet_probability {
					return;
				} else {
					*pixel = mul(survival_boost_factor, *pixel);
				}
				self.force_light_to_eye(intersection, normal, old_direction, &mut maximum_specular_angle, &mut refractive_index, &mut specular_probability, pixel, &mut renderer_output);
				old_direction = direction;
				hit_object = self.compute_ray_object_intersection(&mut intersection, &mut direction, &mut normal, &mut maximum_specular_angle, &mut refractive_index, &mut specular_probability, &mut pixel, &mut pcg, &mut color, false);
				if hit_object {
					*pixel = elementwise_mul(*pixel, color);
				}
			} else {
				return;
			}
		}
	}

	fn force_light_to_eye(&mut self, intersection: [f64; 3], normal: [f64; 3], old_direction: [f64; 3], maximum_specular_angle: &mut f64, refractive_index: &mut f64, specular_probability: &mut f64, pixel: &mut[f64; 3], renderer_output: &mut RendererOutput) {
		let pinhole = [500.0, 500.0, -1000.0];
		let direction_to_retina = normalised(sub(pinhole, intersection));
		if dot(normal, direction_to_retina) < 0.0 {
			return
		}
		let retina_normal = [0.0, 0.0, -1.0];
		let retina_center = [500.0, 500.0, -2000.0];
		let d = dot(sub(retina_center, pinhole), retina_normal)/dot(direction_to_retina, retina_normal);
		let retina_intersection = add(mul(d, direction_to_retina), pinhole);
		if (retina_intersection[2]-retina_center[2]).abs() > 1e-6 {
			return;
		}

		let dist = self.compute_ray_object_distance(intersection, direction_to_retina, true);
		let distance_to_pinhole = dot(sub(pinhole, intersection), retina_normal)/dot(direction_to_retina, retina_normal);
		if dist < distance_to_pinhole {
			return;
		}
		let incoming_direction = mul(-1.0, old_direction);
		if d > 0.0 && retina_intersection[0] > 0.0 && retina_intersection[1] > 0.0 && retina_intersection[0] < (self.width as f64) && retina_intersection[1] < (self.height as f64) {
			let mut pixel_modified = [pixel[0], pixel[1], pixel[2]];
			let f = brdf(incoming_direction, direction_to_retina, normal, *maximum_specular_angle, *refractive_index, *specular_probability);
			let cos_retina = dot(direction_to_retina, retina_normal);
			if cos_retina < 0.0 || cos_retina > 1.0 {
				println!("{:?}, {:?}", direction_to_retina, retina_normal);
				panic!();
			}
			let distance = norm(sub(intersection, pinhole));
			let pinhole_radius = 1.0;
			let alpha = (pinhole_radius/distance).atan();
			let proportion = 1.0 - alpha.cos();
			// Why does 1000000.0 make the pixel intensity have a maximum of seemingly exactly 1.0?
			pixel_modified = mul(1000000.0*f*proportion*(1.0/(cos_retina*cos_retina*cos_retina)), pixel_modified);
			self.store_intersection(retina_intersection, pixel_modified, renderer_output);
		}
	}

	fn store_intersection(&self, retina_intersection: [f64; 3], pixel: [f64; 3], renderer_output: &mut RendererOutput) {
		let y = ((self.height as f64 - retina_intersection[1])*(self.image_scale_factor as f64)) as u32;
		let x = ((self.width as f64 - retina_intersection[0])*(self.image_scale_factor as f64)) as u32;
		if (y*self.width+x) < self.width*self.height*self.image_scale_factor*self.image_scale_factor {
			let pos = (y*self.width*self.image_scale_factor+x) as usize;
			renderer_output.number_of_rays[pos] += 1;
			// Store information about the ray.
			renderer_output.pixels[pos].bins[(((pixel[0]*2.0).atan()*2.0/PI)*(NUMBER_OF_BINS as f64)) as usize][0] += 1;
			renderer_output.pixels[pos].bins[(((pixel[1]*2.0).atan()*2.0/PI)*(NUMBER_OF_BINS as f64)) as usize][1] += 1;
			renderer_output.pixels[pos].bins[(((pixel[2]*2.0).atan()*2.0/PI)*(NUMBER_OF_BINS as f64)) as usize][2] += 1;
			renderer_output.pixels[pos].color = add(renderer_output.pixels[pos].color, pixel);
			renderer_output.colors[pos] = add(renderer_output.colors[pos], pixel);
		} else {
			println!("Can't write to pixels[{}] (x = {}, y = {}), color = {:?}", y*self.width*self.image_scale_factor+x, retina_intersection[0], retina_intersection[1], pixel);
		}
	}

	fn compute_ray_object_intersection(&self, intersection: &mut[f64; 3], direction: &mut[f64; 3], normal: &mut[f64; 3], maximum_specular_angle: &mut f64, refractive_index: &mut f64, specular_probability: &mut f64, pixel: &mut[f64; 3], pcg: &mut Pcg32, color: &mut[f64; 3], camera_ray: bool) -> bool {
		let (triangle_hit, triangle_index, triangle_distance) = self.compute_ray_triangle_intersection(*intersection, *direction, camera_ray);
		let (sphere_hit, sphere_index, sphere_distance) = self.compute_ray_sphere_intersection(*intersection, *direction, camera_ray);
		// The ray hit nothing.
		if !triangle_hit && !sphere_hit {
			return false;
		}
		let incoming_direction = mul(-1.0, *direction);
		// What did the ray hit? A triangle or a sphere?
		if (triangle_hit && !sphere_hit) || (triangle_hit && sphere_hit && triangle_distance < sphere_distance) {
			*intersection = add(*intersection, mul(triangle_distance, *direction));
			*normal = self.scene.triangle_surfaces[triangle_index].normal;
			/*
			*direction = pick_reflection_uniform(*normal, pcg);
			*color = self.scene.triangle_surfaces[triangle_index].compute_intersection_color(intersection);
			*color = mul(dot(*normal, *direction)*2.0, *color);
			*/
			/*
			*direction = pick_reflection_lambertian(*normal, pcg);
			*color = self.scene.triangle_surfaces[triangle_index].compute_intersection_color(intersection);
			*/
			/*
			*specular_probability = self.scene.triangle_surfaces[triangle_index].specular_probability;
			*refractive_index = self.scene.triangle_surfaces[triangle_index].refractive_index;
			*maximum_specular_angle = self.scene.triangle_surfaces[triangle_index].maximum_specular_angle;
			*direction = pick_reflection_uniform(*normal, pcg);
			if dot(incoming_direction, *normal) < 0.0 {
				//println!("Wrong direction in compute_ray_object_intersection, triangle, {}.", dot(incoming_direction, *normal));
				return false;
			}
			let brdf = brdf(incoming_direction, *direction, *normal, *maximum_specular_angle, *refractive_index, *specular_probability);
			*color = self.scene.triangle_surfaces[triangle_index].compute_intersection_color(intersection);
			*color = mul(brdf, *color);
			*/
			*specular_probability = self.scene.triangle_surfaces[triangle_index].specular_probability;
			*refractive_index = self.scene.triangle_surfaces[triangle_index].refractive_index;
			*maximum_specular_angle = self.scene.triangle_surfaces[triangle_index].maximum_specular_angle;
			*direction = pick_reflection_from_brdf(incoming_direction, *direction, *normal, *maximum_specular_angle, *refractive_index, *specular_probability, pcg);
			if dot(incoming_direction, *normal) < 0.0 {
				return false;
			}
			*color = self.scene.triangle_surfaces[triangle_index].compute_intersection_color(intersection);
		} else {
			*intersection = add(*intersection, mul(sphere_distance, *direction));
			*normal = normalised(sub(*intersection, self.scene.sphere_surfaces[sphere_index].center));
			/*
			*direction = pick_reflection_uniform(*normal, pcg);
			*color = self.scene.sphere_surfaces[sphere_index].compute_intersection_color(intersection);
			*color = mul(dot(*normal, *direction)*2.0, *color);
			*/
			/*
			*direction = pick_reflection_lambertian(*normal, pcg);
			*color = self.scene.sphere_surfaces[sphere_index].compute_intersection_color(intersection);
			*/
			/*
			*specular_probability = self.scene.sphere_surfaces[sphere_index].specular_probability;
			*refractive_index = self.scene.sphere_surfaces[sphere_index].refractive_index;
			*maximum_specular_angle = self.scene.sphere_surfaces[sphere_index].maximum_specular_angle;
			*direction = pick_reflection_uniform(*normal, pcg);
			if dot(incoming_direction, *normal) < 0.0 {
				//println!("Wrong direction in compute_ray_object_intersection, sphere, {}.", dot(incoming_direction, *normal));
				return false;
			}
			let brdf = brdf(incoming_direction, *direction, *normal, *maximum_specular_angle, *refractive_index, *specular_probability);
			*color = self.scene.sphere_surfaces[sphere_index].compute_intersection_color(intersection);
			*color = mul(brdf, *color);
			*/
			*specular_probability = self.scene.sphere_surfaces[sphere_index].specular_probability;
			*refractive_index = self.scene.sphere_surfaces[sphere_index].refractive_index;
			*maximum_specular_angle = self.scene.sphere_surfaces[sphere_index].maximum_specular_angle;
			*direction = pick_reflection_from_brdf(incoming_direction, *direction, *normal, *maximum_specular_angle, *refractive_index, *specular_probability, pcg);
			if dot(incoming_direction, *normal) < 0.0 {
				return false;
			}
			*color = self.scene.sphere_surfaces[sphere_index].compute_intersection_color(intersection);
		}
		return true;
	}

	fn compute_ray_object_distance(&self, intersection: [f64; 3], direction: [f64; 3], camera_ray: bool) -> f64 {
		let (_, _, triangle_distance) = self.compute_ray_triangle_intersection(intersection, direction, camera_ray);
		let (_, _, sphere_distance) = self.compute_ray_sphere_intersection(intersection, direction, camera_ray);
		if triangle_distance < sphere_distance {
			triangle_distance
		} else {
			sphere_distance
		}
	}

	fn compute_ray_triangle_intersection(&self, intersection: [f64; 3], direction: [f64; 3], camera_ray: bool) -> (bool, usize, f64) {
		let mut triangle_hit = false;
		let mut triangle_index = 0;
		let mut triangle_distance = f64::MAX;
		for (i, t) in self.scene.triangle_surfaces.iter().enumerate() {
			if camera_ray && t.invisible_for_camera_ray {
				continue;
			}
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
			}
			triangle_hit = true;
			triangle_index = i;
			triangle_distance = d;
		}
		(triangle_hit, triangle_index, triangle_distance)
	}

	fn compute_ray_sphere_intersection(&self, intersection: [f64; 3], direction: [f64; 3], camera_ray: bool) -> (bool, usize, f64) {
		let mut sphere_hit = false;
		let mut sphere_index = 0;
		let mut sphere_distance = f64::MAX;
		// Compute information about the closest hit.
		// Checks for intersection of any sphere (lightsource or not).
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
			}
		}
		(sphere_hit, sphere_index, sphere_distance)
	}
}
