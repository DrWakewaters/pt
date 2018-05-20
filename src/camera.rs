use pcg_rand::Pcg32;

use math::{add, dot, mul, random_uniform_on_sphere, sub};
use ray::Ray;

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct Camera {
	pub pinhole: [f64; 3],
	pub retina_center: [f64; 3],
	pub retina_normal: [f64; 3],
	pub color: [f64; 3],
}

impl Camera {
	pub fn new(pinhole: [f64; 3], retina_center: [f64; 3], retina_normal: [f64; 3], color: [f64; 3]) -> Self {
		Self {
			pinhole,
			retina_center,
			retina_normal,
			color,
		}
	}

	pub fn direction(&self, pcg: &mut Pcg32) -> [f64; 3] {
		random_uniform_on_sphere(pcg)
	}

	pub fn retina_position(&self, ray: &mut Ray) -> Option<usize> {
		let distance = dot(sub(self.retina_center, self.pinhole), self.retina_normal)/dot(mul(-1.0, ray.direction), self.retina_normal);
		if distance < 0.0 {
			return None;
		}
		let point_on_retina = add(mul(-1.0*distance, ray.direction), self.pinhole);
		if (point_on_retina[2] - self.retina_center[2]).abs() > 1.0e-6 {
			println!("Error: {} != {}", point_on_retina[2], self.retina_center[2]);
			return None;
		}
		if point_on_retina[0] < 0.0 || point_on_retina[1] < 0.0 || point_on_retina[0] >= 1000.0 || point_on_retina[1] >= 1000.0 {
			return None;
		}
		let x = (1000.0-point_on_retina[0]) as usize;
		let y = (1000.0-point_on_retina[1]) as usize;
		// @TODO self.width or something of the sort, instead of 1000.
		Some(y*(1000 as usize)+x)
	}
}
