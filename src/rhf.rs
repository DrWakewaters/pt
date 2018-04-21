use std::f32::MAX;

use math::{add, add_f32, mul};
use NUMBER_OF_BINS;

pub struct RHF {
	pub width: i32,
	pub height: i32,
	pub image_scale_factor: i32,
	pub max_distance: f32,
	pub patch_radius: i32,
	pub search_window_radius: i32,
	pub number_of_scales: i32,
	pub vertical_start: i32,
	pub vertical_end: i32,
	pub number_of_rays: Vec<u32>,
	pub bins: Vec<[[f32; 3]; NUMBER_OF_BINS]>,
	pub colors: Vec<[f64; 3]>,
	pub colors_denoised: Vec<[f64; 3]>,
	pub counter: Vec<f64>,
}

impl RHF {
	pub fn new(width: i32, height: i32, image_scale_factor: i32, max_distance: f32, patch_radius: i32, search_window_radius: i32, number_of_scales: i32, vertical_start: i32, vertical_end: i32, number_of_rays: Vec<u32>, colors: Vec<[f64; 3]>, bins: Vec<[[f32; 3]; NUMBER_OF_BINS]>) -> Self {	
		let mut colors_denoised: Vec<[f64; 3]> = Vec::new();
		let mut counter: Vec<f64> = Vec::new();
		for _ in 0..width*height*image_scale_factor*image_scale_factor {
			colors_denoised.push([0.0, 0.0, 0.0]);
			counter.push(0.0);
		}
		Self {
			width,
			height,
			image_scale_factor,
			max_distance,
			patch_radius,
			search_window_radius,
			number_of_scales,
			vertical_start,
			vertical_end,
			bins,
			number_of_rays,
			colors,
			colors_denoised,
			counter,
		}
	}
	
	pub fn rhf(&mut self) -> Vec<[f64; 3]> {
		/*
		let indices_to_study: Vec<usize> = vec![(768*self.width*self.image_scale_factor + 288) as usize, (769*self.width*self.image_scale_factor + 300) as usize, (770*self.width*self.image_scale_factor + 312) as usize, (771*self.width*self.image_scale_factor + 324) as usize, (411*self.width*self.image_scale_factor + 576) as usize, (247*self.width*self.image_scale_factor + 433) as usize];
		let mut influencer_indices: Vec<Vec<usize>> = Vec::new();
		for _ in 0..indices_to_study.len() {
			influencer_indices.push(Vec::new());
		}
		*/
		for y_0 in self.vertical_start..self.vertical_end {
			for x_0 in 0..self.width*self.image_scale_factor  {
				for y_1 in y_0-self.search_window_radius..y_0+self.search_window_radius+1 {
					for x_1 in x_0-self.search_window_radius..x_0+self.search_window_radius+1 {
						let pixel_indices = self.compute_pixel_indices(y_0, x_0, y_1, x_1);
						let distance = self.kolmogorov_smirnov_distance(&pixel_indices);
						// TODO: distance will unfortunately be lower nere the edges when we have patches. This should be taken care of.
						if distance < self.max_distance {
							for pixel_index in pixel_indices {
								self.colors_denoised[pixel_index.0 as usize] = add(self.colors_denoised[pixel_index.0 as usize], self.colors[pixel_index.1 as usize]);
								self.counter[pixel_index.0 as usize] += 1.0;
							}
						}
					}
				}
			}
		}
		let first_index = (self.vertical_start*self.width*self.image_scale_factor) as usize;
		let last_index = (self.vertical_end*self.width*self.image_scale_factor) as usize;
		let mut color_result: Vec<[f64; 3]> = Vec::new();
		for i in first_index..last_index {
			if self.counter[i].abs() < 1e-6 {
				color_result.push(self.colors[i]);
			} else {
				color_result.push(mul(1.0/self.counter[i], self.colors_denoised[i]));
			}
		}
		color_result
	}	
	pub fn kolmogorov_smirnov_distance(&self, pixel_indices: &Vec<(i32, i32)>) -> f32 {
		let mut max_distances: Vec<f32> = Vec::new();
		for pixel_index in pixel_indices {
			let mut cummulative_0 = [0.0, 0.0, 0.0];
			let mut cummulative_1 = [0.0, 0.0, 0.0];
			let mut max_distance = 0.0;
			let bins_0 = &self.bins[pixel_index.0 as usize];
			let bins_1 = &self.bins[pixel_index.1 as usize];
			let number_of_rays_0 = self.number_of_rays[pixel_index.0 as usize];
			let number_of_rays_1 = self.number_of_rays[pixel_index.1 as usize];
			if number_of_rays_0 == 0 || number_of_rays_1 == 0 {
				max_distances.push(MAX);
				continue;
			}
			for i in 0..NUMBER_OF_BINS {
				cummulative_0 = add_f32(cummulative_0, bins_0[i]);
				cummulative_1 = add_f32(cummulative_1, bins_1[i]);
				let distance = (cummulative_1[0]-cummulative_0[0]).abs()+(cummulative_1[1]-cummulative_0[1]).abs()+(cummulative_1[2]-cummulative_0[2]).abs();
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

	pub fn compute_pixel_indices(&self, y_0_middle: i32, x_0_middle: i32, y_1_middle: i32, x_1_middle: i32) -> Vec<(i32, i32)> {
		let mut pixel_indices: Vec<(i32, i32)> = Vec::new();
		for vertical_shift in -self.patch_radius..self.patch_radius+1 {
			for horizontal_shift in -self.patch_radius..self.patch_radius+1 {
				let y_0_shifted = y_0_middle+vertical_shift;
				let y_1_shifted = y_1_middle+vertical_shift;
				let x_0_shifted = x_0_middle+horizontal_shift;
				let x_1_shifted = x_1_middle+horizontal_shift;
				if y_0_shifted < 0 || y_1_shifted < 0 || y_0_shifted >= self.height*self.image_scale_factor || y_1_shifted >= self.height*self.image_scale_factor {
					continue;
				}
				if x_0_shifted < 0 || x_1_shifted < 0 || x_0_shifted >= self.width*self.image_scale_factor || x_1_shifted >= self.width*self.image_scale_factor {
					continue;
				}
				pixel_indices.push((y_0_shifted*self.width*self.image_scale_factor+x_0_shifted, y_1_shifted*self.width*self.image_scale_factor+x_1_shifted));
			}
		}
		pixel_indices
	}
}
