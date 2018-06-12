use rendereroutputrow::RendererOutputRow;

use NUMBER_OF_BINS;

#[derive(Clone, Serialize, Deserialize)]
pub struct RendererOutput {
	pub width: u32,
	pub height: u32,
	pub bins: Vec<Vec<[f32; 3]>>,
	pub colors: Vec<[f64; 3]>,
	pub number_of_rays: Vec<f64>,
}

impl RendererOutput {
	pub fn new(width: u32, height: u32, mut renderer_output_rows: Vec<RendererOutputRow>) -> Self {
		let mut bins: Vec<Vec<[f32; 3]>> = Vec::new();
		let mut colors: Vec<[f64; 3]> = Vec::new();
		let mut number_of_rays: Vec<f64> = Vec::new();
		for i in 0..height {
			for j in 0..width {
				colors.push(renderer_output_rows[i as usize].colors[j as usize]);
				let mut bin = vec![[0.0; 3]; NUMBER_OF_BINS];
				let mut number_of_bin_elements = 0; // Not identical to number_of_rays. Maybe we should compute this value in renderer instead.
				for k in 0..bin.len() {
					number_of_bin_elements += renderer_output_rows[i as usize].pixels[j as usize].bins[k as usize][0];
				}
				for (k, color) in bin.iter_mut().enumerate() {
					for (l, color_component) in color.iter_mut().enumerate() {
						*color_component = (1.0/f64::from(number_of_bin_elements)*f64::from(renderer_output_rows[i as usize].pixels[j as usize].bins[k as usize][l])) as f32;
					}
				}
				bins.push(bin);
			}
		}
		for i in 0..height {
			number_of_rays.append(&mut renderer_output_rows[i as usize].number_of_rays);
		}
		Self {
			width,
			height,
			bins,
			colors,
			number_of_rays,
		}
	}
}
