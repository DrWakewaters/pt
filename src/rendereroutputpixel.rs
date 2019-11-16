use crate::pixel::Pixel;

#[derive(Clone)]
pub struct RendererOutputPixel {
	pub pixel: Pixel,
	pub color: [f64; 3],
	pub number_of_rays: f64,
	pub number_of_bin_elements: u32,
	pub y: u32,
	pub x: u32,
}

impl RendererOutputPixel {
	pub fn new(y: u32, x: u32) -> Self {
		Self {
			pixel: Pixel::new(),
			color: [0.0, 0.0, 0.0],
			number_of_rays: 0.0,
			number_of_bin_elements: 0,
			y,
			x,
		}
	}
}
