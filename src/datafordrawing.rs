use std::f64::consts::PI;

use dataforstoring::DataForStoring;

#[derive(Serialize, Deserialize, Debug)]
pub struct DataForDrawing {
	pub width: u32,
	pub height: u32,
	number_of_rays_total: u64,
	pub pixels: Vec<u8>
}

impl DataForDrawing {
	pub fn new(data_for_storing: DataForStoring) -> Self {
		let mut pixels: Vec<u8> = Vec::new();
		for pixel in data_for_storing.pixels {
			let mut max_intensity = pixel[0];
			if pixel[1] > max_intensity {
				max_intensity = pixel[1];
			}
			if pixel[2] > max_intensity {
				max_intensity = pixel[2];
			}
			let factor = if max_intensity > 1.0e-9 {
				let modified_max_intensity = (max_intensity/((data_for_storing.number_of_rays_total as f64)/2_000_000.0)).atan()*255.0/(PI/2.0);
				modified_max_intensity/max_intensity
			} else {
				1.0
			};
			pixels.push((pixel[0]*factor) as u8);
			pixels.push((pixel[1]*factor) as u8);
			pixels.push((pixel[2]*factor) as u8);
		}
		DataForDrawing {
			width: data_for_storing.width,
			height: data_for_storing.height,
			number_of_rays_total: data_for_storing.number_of_rays_total,
			pixels,
		}
	}
}
