use serde_derive::{Serialize, Deserialize};

use rendereroutput::RendererOutput;

#[derive(Serialize, Deserialize, Debug)]
pub struct DataForStoring {
	width: u32,
	height: u32,
	number_of_rays: Vec<u32>,
	colors: Vec<[f64; 3]>,
}

impl DataForStoring {
	pub fn new(width: u32, height: u32, renderer_output: &mut RendererOutput) -> Self {
		let mut number_of_rays: Vec<u32> = Vec::new();
		let mut colors: Vec<[f64; 3]> = Vec::new();
		for renderer_output_row in &mut renderer_output.renderer_output_rows {
			number_of_rays.append(&mut renderer_output_row.number_of_rays);
			//@TODO: Check if this is correct.
			colors.append(&mut renderer_output_row.colors);
		}
		DataForStoring {
			width,
			height,
			number_of_rays,
			colors,
		}
	}
}
