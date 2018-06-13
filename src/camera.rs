use std::f64::consts::PI;

use rand::Rng;
use pcg_rand::Pcg32;

use math::{add, dot, mul, normalised, sub};
use ray::Ray;

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct Camera {
	pinhole: [f64; 3],
	retina_center: [f64; 3],
	retina_normal: [f64; 3],
	pub color: [f64; 3],
	point_on_focal_plane: [f64; 3],
	pinhole_radius: f64,
}

impl Camera {
	pub fn new(pinhole: [f64; 3], retina_center: [f64; 3], retina_normal: [f64; 3], color: [f64; 3], focal_length: f64, pinhole_radius: f64) -> Self {
		Self {
			pinhole,
			retina_center,
			retina_normal,
			color,
			point_on_focal_plane: add(pinhole, mul(focal_length, retina_normal)),
			pinhole_radius,
		}
	}

	pub fn create_ray(&self, x: u32, y: u32, width: u32, height: u32, pcg: &mut Pcg32) -> Ray {
		let r = pcg.next_f64()*0.5;
		let theta = pcg.next_f64()*2.0*PI;
		let dx = theta.cos()*r;
		let dy = theta.sin()*r;
		let point_on_retina = add(self.retina_center, [f64::from(width)/2.0 - f64::from(x) + dx, f64::from(height)/2.0 - f64::from(y) + dy, 0.0]);
		let mut direction = normalised(sub(self.pinhole, point_on_retina));
		let distance_to_focal_plane = dot(sub(self.point_on_focal_plane, self.pinhole), self.retina_normal)/dot(direction, self.retina_normal);
		let point_on_focal_plane = add(self.pinhole, mul(distance_to_focal_plane, direction));
		let pinhole_translation = mul(self.pinhole_radius, [pcg.next_f64(), pcg.next_f64(), 0.0]);
		let point_on_lens = add(self.pinhole, pinhole_translation);
		direction = normalised(sub(point_on_focal_plane, point_on_lens));
		Ray::new(point_on_lens, direction)
	}
	/*
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
	*/
}
