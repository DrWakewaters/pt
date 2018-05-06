use pixel::Pixel;

#[derive(Clone)]
pub struct RendererOutput {
	pub number_of_rays: Vec<u32>,
	pub pixels: Vec<Pixel>,
	pub colors: Vec<[f64; 3]>,
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
