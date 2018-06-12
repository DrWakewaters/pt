use pixel::Pixel;

#[derive(Clone)]
pub struct RendererOutputRow {
	pub number_of_rays: Vec<f64>,
	pub pixels: Vec<Pixel>,
	pub colors: Vec<[f64; 3]>,
	pub row_number: u32,
}

impl RendererOutputRow {
	pub fn new(row_number: u32, row_width: u32) -> Self {
		let mut number_of_rays: Vec<f64> = Vec::new();
		let mut pixels: Vec<Pixel> = Vec::new();
		let mut colors: Vec<[f64; 3]> = Vec::new();
		for _ in 0..row_width {
			number_of_rays.push(0.0);
			pixels.push(Pixel::new());
			colors.push([0.0, 0.0, 0.0]);
		}
		Self {
			number_of_rays,
			pixels,
			colors,
			row_number,
		}
	}
}
