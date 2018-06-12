use NUMBER_OF_BINS;

#[derive(Clone)]
pub struct Pixel {
	pub bins: Vec<[u16; 3]>,
	pub color: [f64; 3],
}

impl Pixel {
	pub fn new() -> Self {
		Self {
			bins: vec![[0; 3]; NUMBER_OF_BINS],
			color: [0.0; 3],
		}
	}
}
