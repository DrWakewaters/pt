use std::f64::MAX;

use crate::math::{add, add_f32, mul};
use crate::rendereroutput::RendererOutput;

use crate::NUMBER_OF_BINS;

pub struct RHF {
	width: i32,
	height: i32,
	max_distance: f64,
	patch_radius: i32,
	search_window_radius: i32,
	number_of_rays: Vec<f64>,
	bins: Vec<Vec<[f32; 3]>>,
	colors: Vec<[f64; 3]>,
	colors_denoised: Vec<[f64; 3]>,
	number_of_rays_denoised: Vec<f64>,
}

impl RHF {
	pub fn new(width: i32, height: i32, max_distance: f64, patch_radius: i32, search_window_radius: i32, renderer_output: &mut RendererOutput) -> Self {
		let mut colors_denoised: Vec<[f64; 3]> = Vec::new();
		let mut number_of_rays_denoised: Vec<f64> = Vec::new();
		for _ in 0..width*height {
			colors_denoised.push([0.0, 0.0, 0.0]);
			number_of_rays_denoised.push(0.0);
		}
		Self {
			width,
			height,
			max_distance,
			patch_radius,
			search_window_radius,
			bins: renderer_output.bins.to_vec(),
			number_of_rays: renderer_output.number_of_rays.to_vec(),
			colors: renderer_output.colors.to_vec(),
			colors_denoised,
			number_of_rays_denoised,
		}
	}

	pub fn rhf(&mut self) -> (Vec<[f64; 3]>, Vec<f64>) {
		for y_0 in 0..self.height {
			for x_0 in 0..self.width {
				for y_1 in y_0-self.search_window_radius..y_0+self.search_window_radius+1 {
					for x_1 in x_0-self.search_window_radius..x_0+self.search_window_radius+1 {
						let pixel_indices = self.compute_pixel_indices(y_0, x_0, y_1, x_1);
						let distance = self.kolmogorov_smirnov_distance(&pixel_indices);
						// TODO: distance will unfortunately be lower nere the edges when we have patches. This should be taken care of.
						if distance < self.max_distance {
							for pixel_index in pixel_indices {
								let weight = 1.0/(99.0*(distance/(self.max_distance as f64)).sqrt().sqrt() + 1.0);
								self.colors_denoised[pixel_index.0 as usize] = add(self.colors_denoised[pixel_index.0 as usize], mul(weight, self.colors[pixel_index.1 as usize]));
								self.number_of_rays_denoised[pixel_index.0 as usize] += weight*self.number_of_rays[pixel_index.1 as usize];
							}
						}
					}
				}
			}
		}
		(self.colors_denoised.to_vec(), self.number_of_rays_denoised.to_vec())
	}
	fn kolmogorov_smirnov_distance(&self, pixel_indices: &[(i32, i32)]) -> f64 {
		let mut max_distances: Vec<f64> = Vec::new();
		for pixel_index in pixel_indices {
			let mut cummulative_0 = [0.0, 0.0, 0.0];
			let mut cummulative_1 = [0.0, 0.0, 0.0];
			let mut max_distance = 0.0;
			let bins_0 = &self.bins[pixel_index.0 as usize];
			let bins_1 = &self.bins[pixel_index.1 as usize];
/*
			println!("pixel_index = ({}, {})", pixel_index.0, pixel_index.1);
			for i in 0..NUMBER_OF_BINS {
				println!("bins[{}] = {:?}", i, bins_0[i as usize]);
			}
			println!("");
*/
			let number_of_rays_0 = self.number_of_rays[pixel_index.0 as usize];
			let number_of_rays_1 = self.number_of_rays[pixel_index.1 as usize];
			if number_of_rays_0 < 1.0e-6 || number_of_rays_1 < 1.0e-6 {
				max_distances.push(MAX);
				continue;
			}
			for i in 0..NUMBER_OF_BINS {
				cummulative_0 = add_f32(cummulative_0, bins_0[i]);
				cummulative_1 = add_f32(cummulative_1, bins_1[i]);
				let distance = f64::from((cummulative_1[0]-cummulative_0[0]).abs()+(cummulative_1[1]-cummulative_0[1]).abs()+(cummulative_1[2]-cummulative_0[2]).abs());
				if distance > max_distance {
					max_distance = distance;
				}
			}
			max_distances.push(max_distance);
		}
		let mut max_distance_total = 0.0;
		for max_distance in max_distances {
			max_distance_total += max_distance;
		}
		max_distance_total
	}

	fn compute_pixel_indices(&self, y_0_middle: i32, x_0_middle: i32, y_1_middle: i32, x_1_middle: i32) -> Vec<(i32, i32)> {
		let mut pixel_indices: Vec<(i32, i32)> = Vec::new();
		for vertical_shift in -self.patch_radius..self.patch_radius+1 {
			for horizontal_shift in -self.patch_radius..self.patch_radius+1 {
				let y_0_shifted = y_0_middle+vertical_shift;
				let y_1_shifted = y_1_middle+vertical_shift;
				let x_0_shifted = x_0_middle+horizontal_shift;
				let x_1_shifted = x_1_middle+horizontal_shift;
				if y_0_shifted < 0 || y_1_shifted < 0 || y_0_shifted >= self.height || y_1_shifted >= self.height {
					continue;
				}
				if x_0_shifted < 0 || x_1_shifted < 0 || x_0_shifted >= self.width || x_1_shifted >= self.width {
					continue;
				}
				pixel_indices.push((y_0_shifted*self.width+x_0_shifted, y_1_shifted*self.width+x_1_shifted));
			}
		}
		pixel_indices
	}
}
