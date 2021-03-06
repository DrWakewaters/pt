use serde_derive::{Serialize, Deserialize};

use crate::material::Material;
use crate::physics::Physics;

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct PhysicsSphere {
	pub position: [f64; 3],
	pub radius: f64,
	pub material: Material,
	pub physics: Physics,
	pub id: i64,
	pub active: bool
}

impl PhysicsSphere {
	pub fn new(position: [f64; 3], radius: f64, material: Material, physics: Physics, id: i64, active: bool) -> Self {
		Self {
			position,
			radius,
			material,
			physics,
			id,
			active
		}
	}
}
