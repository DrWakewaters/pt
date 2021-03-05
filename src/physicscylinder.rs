use serde_derive::{Serialize, Deserialize};

use crate::material::Material;
use crate::math::normalised;
use crate::physics::Physics;

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct PhysicsCylinder {
	pub position: [f64; 3],
    pub direction: [f64; 3],
    pub length: f64,
	pub radius: f64,
	pub material: Material,
	pub physics: Physics,
	pub id: i64,
	pub active: bool
}

impl PhysicsCylinder {
	pub fn new(position: [f64; 3], direction: [f64; 3], length: f64, radius: f64, material: Material, physics: Physics, id: i64, active: bool) -> Self {
		Self {
			position,
            direction: normalised(direction),
            length,
			radius,
			material,
			physics,
			id,
			active
		}
	}
}
