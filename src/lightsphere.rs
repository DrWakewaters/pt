use serde_derive::{Serialize, Deserialize};

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
}
