use pcg_rand::Pcg32;
use serde_derive::{Serialize, Deserialize};

use crate::math::{add, mul, random_uniform_on_sphere};

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct LightSphere {
    pub position: [f64; 3],
    pub color: [f64; 3],
    pub radius: f64,
}

impl LightSphere {
    pub fn new(position: [f64; 3], color: [f64; 3], radius: f64) -> Self {
        Self {
            position,
            color,
            radius,
        }
    }

    pub fn get_position(&self, pcg: &mut Pcg32) -> [f64; 3] {
        add(self.position, mul(self.radius, random_uniform_on_sphere(pcg)))
    }
}
