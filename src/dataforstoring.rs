use rendereroutput::RendererOutput;

#[derive(Serialize, Deserialize, Debug)]
pub struct DataForStoring {
	pub width: u32,
	pub height: u32,
	pub number_of_rays_total: u64,
	pub pixels: Vec<[f64; 3]>,
}

impl DataForStoring {
	pub fn new(width: u32, height: u32, number_of_rays_total: u64, renderer_output_aggregate: &mut RendererOutput) -> Self {
		DataForStoring {
			width,
			height,
			number_of_rays_total,
			pixels: renderer_output_aggregate.colors.clone(),
		}
	}
}
