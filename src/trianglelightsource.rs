use pcg_rand::Pcg32;
use rand::Rng;

use math::{add, mul, point_in_triangle};
use ray::Ray;
use trianglephysics::TrianglePhysics;

#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
pub struct TriangleLightsource {
	pub node0: [f64; 3],
	pub node1: [f64; 3],
	pub node2: [f64; 3],
	pub e1: [f64; 3],
	pub e2: [f64; 3],
	pub normal: [f64; 3],
	pub color: [f64; 3],
}

impl TriangleLightsource {
	pub fn new(triangle_physics: &TrianglePhysics) -> Self {
		Self {
			node0: triangle_physics.node0,
			node1: triangle_physics.node1,
			node2: triangle_physics.node2,
			e1: triangle_physics.e1,
			e2: triangle_physics.e2,
			normal: triangle_physics.normal,
			color: triangle_physics.color,
		}
	}

	// @TODO Implement non-normal emission.
	pub fn compute_ray_data(&self, mut pcg: &mut Pcg32) -> Ray {
		loop {
			let point = add(add(mul(pcg.next_f64(), self.e1), mul(pcg.next_f64(), self.e2)), self.node0);
			if point_in_triangle(point, self.node0, self.node1, self.node2) {
				return Ray::new(point, self.normal, self.normal, self.normal, self.color, 0.0, 1.0, 0.0);
			}
		}
	}
}
