use NUMBER_OF_BINS;

#[derive(Clone)]
pub struct Pixel {
	pub bins: [[u16; 3]; NUMBER_OF_BINS],
	pub color: [f64; 3],
}

impl Pixel {
	pub fn new() -> Self {
		Self {
			bins: [[0; 3]; NUMBER_OF_BINS],
			color: [0.0; 3],
		}
	}
}
