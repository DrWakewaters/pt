use serde_derive::{Serialize, Deserialize};

use crate::math::{intensity_to_color, mul};
use crate::rendereroutput::RendererOutput;

#[derive(Serialize, Deserialize, Debug)]
pub struct DataForDrawing {
	pub width: u32,
	pub height: u32,
	pub colors: Vec<u8>
}

impl DataForDrawing {
	pub fn new(renderer_output: &mut RendererOutput) -> Self {
		let mut colors: Vec<u8> = Vec::new();
		for (i, color) in renderer_output.colors.iter().enumerate() {
			let color_per_ray = intensity_to_color(mul(1.0/(renderer_output.number_of_rays[i] as f64), *color));
			colors.push(color_per_ray[0] as u8);
			colors.push(color_per_ray[1] as u8);
			colors.push(color_per_ray[2] as u8);
		}
		DataForDrawing {
			width: renderer_output.width,
			height: renderer_output.height,
			colors,
		}
	}
}
