use std::f64;

use pcg_rand::Pcg32;
use rand::{Rng, SeedableRng};

use crate::hitpoint::Hitpoint;
use crate::material::Material;
use crate::math::{add, dot, elementwise_mul, intensity_to_color, mul, random_from_brdf};
use crate::renderershape::RendererShape;
use crate::rendereroutputpixel::RendererOutputPixel;
use crate::rendererscene::RendererScene;

use crate::NUMBER_OF_BINS;

//static MODULO: [usize; 5] = [0, 1, 2, 0, 1];

pub struct Renderer {
	width: u32,
	height: u32,
	spp_per_iteration: u32,
	maximum_spp: u32,
	maximum_error: f64,
	perform_post_process: bool,
	scene: RendererScene,
	renderer_output_pixels: Vec<RendererOutputPixel>
}

impl Renderer {
	pub fn new(width: u32, height: u32, spp_per_iteration: u32, maximum_spp: u32, maximum_error: f64, perform_post_process: bool, scene: RendererScene) -> Self {
		let renderer_output_pixels: Vec<RendererOutputPixel> = Vec::new();
		Self {
			width,
			height,
			spp_per_iteration,
			maximum_spp,
			maximum_error,
			perform_post_process,
			scene,
			renderer_output_pixels,
		}
	}

	pub fn get_renderer_output_pixels(&self) -> Vec<RendererOutputPixel> {
		self.renderer_output_pixels.to_vec()
	}

	pub fn render(&mut self, y: u32, x: u32) {
		if x == self.width/2 {
			println!("Rendering row {} of {}.", y, self.height);
		}
		self.renderer_output_pixels.push(RendererOutputPixel::new(y, x));
		let last_pos = self.renderer_output_pixels.len()-1;
		let mut pcg: Pcg32 = Pcg32::from_entropy();
		let mut colors: Vec<[f64; 3]> = Vec::new();
		let mut converged = false;
		// Loop over this pixel until we estimate the error to be small enough.
		let mut iterations = 0;
		while !converged {
			iterations += 1;
			for _ in 0..self.spp_per_iteration {
				// Create a hitpoint path starting from the camera. The first hitpoint will be the hitpoint of the ray coming from the retina, through the pinhole of the camera - that is, it is a point somewhere on a renderer shape.
				self.renderer_output_pixels[last_pos].number_of_rays += 1.0;
				let (position, direction) = self.scene.camera.create_ray(x, y, self.width, self.height, &mut pcg);
				let color = self.compute_color(position, direction, self.scene.camera.color, &mut pcg);
				self.store(last_pos, color);
				colors.push(color);
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
			if error < self.maximum_error {
				converged = true;
			} else if iterations*self.spp_per_iteration >= self.maximum_spp {
				converged = true;
				println!("Failed to converge at ({}, {}). The error is {} whereas the goal was {}.", x, y, error, self.maximum_error);
			} else if iterations%100 == 0 {
				println!("Currently: {} iterations at ({}, {}) with an error of {}. Iterating until it is less than {}.", iterations, x, y, error, self.maximum_error);
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
		let variance = (r_squared + g_squared + b_squared - (r*r+g*g+b*b)/(length as f64))/(length as f64-1.0);
		// Due to rounding erros, the computed variance could in rare cases be slightly lower than 0.0.
		if variance < 0.0 {
			0.0
		} else {
			variance.sqrt()
		}
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

	fn compute_color(&mut self, position: [f64; 3], direction: [f64; 3], color: [f64; 3], pcg: &mut Pcg32) -> [f64; 3] {
		let mut color = color;
		let bullet_probability = 0.0;
		let survival_boost_factor = 1.0/(1.0-bullet_probability);
		let mut hitpoint = Hitpoint::new(position, direction, [0.0, 0.0, 0.0], Material::none(), true, false);
		let light_sampling_probability = 0.0;
		loop {
			let r = pcg.gen::<f64>();
			if r < bullet_probability {
				return [0.0, 0.0, 0.0]
			} else {
				color = mul(survival_boost_factor, color);
			}
			let hit_renderer_shape = self.closest_renderer_shape(&mut hitpoint);
			if hit_renderer_shape {
				color = elementwise_mul(color, hitpoint.material.color);
				if hitpoint.material.is_light {
					return color
				}
				let mut light_sampling_possible = false;
				let mut brdf_modifier = 0.0;
				random_from_brdf(&mut hitpoint, &mut self.scene.light_spheres, light_sampling_probability, &mut light_sampling_possible, &mut brdf_modifier, pcg);
				if light_sampling_possible {
					// @TODO: Add support for transmission.
					color = mul(brdf_modifier, color);
				}
			} else {
				return [0.0, 0.0, 0.0]
			}
		}
	}

// rP: probability that an interation is a reflection
// sP: probability that a reflection is specular
// lBr: lambertian_brdf_reflection
// sBr: specular_brdf_reflection
// lBt: lambertian_brdf_transmission
// sBt: specular_brdf_transmission
// aP: light_sampling_probability
// aB: light_brdf

// Return correct/actual (correct probability density / the actual probability density).
// With transmission:
// correct = (1-sP)*lBr*rP + sP*sBr*rP + (1-sP)*lBt*(1-rP) + sP*sBt*(1-rP)
// actual = (1-aP)*correct + aP*aB
// Without transmission:
// correct = (1-sP)*lBr + sP*sBr
// actual = (1-aP)*correct + aP*aB

// I en process kan resultaten 1-7 ges. 1-3: 20 % sannolikhet var, 4-7: 10 % sannolikhet var.
// I själva verket dras dock 1-2 med sannolikheten 25 % var och 3-7 med 10 % sannolikhet var.
// Resultatet blir det tal som väljs. Hur får vi rätt väntevärde?
// Om vi får 1: Mult med 0.2/0.25 = 0.8
// Om vi får 2: Mult med 0.2/0.25 = 0.8
// Om vi får 3: Mult med 0.2/0.1 = 2
// Om vi får 4: Mult med 0.1/0.1 = 1
// Om vi får 5: Mult med 0.1/0.1 = 1
// Om vi får 6: Mult med 0.1/0.1 = 1
// Om vi får 7: Mult med 0.1/0.1 = 1

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
	fn closest_renderer_shape(&self, mut hitpoint: &mut Hitpoint) -> bool  {
		let mut min_distance = f64::MAX;
		let mut closest_renderer_shape_index: Option<usize> = None;
		let mut closest_is_a_sphere = true;
		let position = hitpoint.position;
		let direction = hitpoint.incoming_direction;
		for (i, sphere) in self.scene.renderer_spheres.iter().enumerate() {
			let distance = sphere.distance(position, direction);
			if distance < min_distance {
				min_distance = distance;
				closest_renderer_shape_index = Some(i);
			}
		}
		for (i, triangle) in self.scene.renderer_triangles.iter().enumerate() {
			let distance = triangle.distance(position, direction);
			if distance < min_distance {
				min_distance = distance;
				closest_renderer_shape_index = Some(i);
				closest_is_a_sphere = false;
			}
		}
		if let Some(index) = closest_renderer_shape_index {
			let position = add(position, mul(min_distance, direction));
			let mut normal = if closest_is_a_sphere {
				self.scene.renderer_spheres[index].normal(position)
			} else {
				self.scene.renderer_triangles[index].normal(position)
			};
			let material = if closest_is_a_sphere {
				self.scene.renderer_spheres[index].material()
			} else {
				self.scene.renderer_triangles[index].material()
			};
			if dot(direction, normal) > 0.0 {
				normal = mul(-1.0, normal);
			}
			hitpoint.position = position;
			hitpoint.normal = normal;
			hitpoint.material = material;
			true
		} else {
			false
		}
	}
}
