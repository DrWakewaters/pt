use serde_derive::{Serialize, Deserialize};

use crate::material::Material;
use crate::physics::Physics;

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct PhysicsSphere {
	pub position: [f64; 3],
	pub radius: f64,
	pub material: Material,
	pub physics: Physics,
}

impl PhysicsSphere {
	pub fn new(position: [f64; 3], radius: f64, material: Material, physics: Physics) -> Self {
		Self {
			position,
			radius,
			material,
			physics,
		}
	}
}
