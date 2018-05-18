pub struct Ray {
	pub position: [f64; 3],
	pub direction: [f64; 3],
}

impl Ray {
	pub fn new(position: [f64; 3], direction: [f64; 3]) -> Self {
		Self {
			position,
			direction,
		}
	}
}
