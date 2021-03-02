use std::f64;

use serde_derive::{Serialize, Deserialize};

use crate::material::Material;
use crate::math::{dot, normalised, sub};
use crate::physicssphere::PhysicsSphere;
use crate::renderershape::RendererShape;

#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
pub struct RendererSphere {
	position: [f64; 3],
	radius: f64,
	material: Material,
}

impl RendererSphere {
	pub fn new(physics_sphere: &PhysicsSphere) -> Self {
		Self {
			position: physics_sphere.position,
			radius: physics_sphere.radius,
			material: physics_sphere.material,
		}
	}
}

impl RendererShape for RendererSphere {
	fn distance(&self, position: [f64; 3], direction: [f64; 3]) -> f64 {
		let b = sub(position, self.position);
		let a = dot(b, direction)*dot(b, direction) - dot(b, b) + self.radius*self.radius;
		if a < 1.0e-6 {
			return f64::MAX;
		}
		let d1 = -1.0*dot(sub(position, self.position), direction) + a.sqrt();
		let d2 = -1.0*dot(sub(position, self.position), direction) - a.sqrt();
		if d2 > 1.0e-6 {
			return d2;
		} else if d1 > 1.0e-6 {
			return d1;
		}
		f64::MAX
	}

	fn material(&self) -> Material {
		self.material
	}

	fn normal(&self, point: [f64; 3]) -> [f64; 3] {
		normalised(sub(point, self.position))
	}
}

/*
#[allow(dead_code)]
pub fn compute_color_striped(&self, intersection: &mut[f64; 3]) -> [f64; 3] {
	// Compute the intersection relative to the center of the sphere.
	let intersection_translated = sub(*intersection, self.center);
	// Compute the intersection in the local coordinates of the sphere.
	let local_intersection = solve_linear([self.local_x, self.local_y, self.local_z], intersection_translated);
	self.get_color(local_intersection)
}

#[allow(dead_code)]
pub fn get_color(&self, local_intersection: [f64; 3]) -> [f64; 3] {
	let region_y = ((local_intersection[1]/self.radius)*5.0).floor() as i64;
	match region_y%2 {
		0 => {
			self.material.color
		}
		_ => {
			self.material.color_alt
		}
	}
}
*/
