pub struct Ray {
	pub origin: [f64; 3],
	pub direction: [f64; 3],
	pub old_direction: [f64; 3],
	pub normal_at_origin: [f64; 3],
	pub color: [f64; 3],
	pub maximum_specular_angle_at_origin: f64,
	pub refractive_index_at_origin: f64,
	pub specular_probability_at_origin: f64,
}

impl Ray {
	pub fn new(origin: [f64; 3], direction: [f64; 3], old_direction: [f64; 3], normal_at_origin: [f64; 3], color: [f64; 3], maximum_specular_angle_at_origin: f64, refractive_index_at_origin: f64, specular_probability_at_origin: f64) -> Self {
		Self {
			origin,
			direction,
			old_direction,
			normal_at_origin,
			color,
			maximum_specular_angle_at_origin,
			refractive_index_at_origin,
			specular_probability_at_origin,
		}
	}
}
