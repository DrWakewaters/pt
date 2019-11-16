use serde_derive::{Serialize, Deserialize};

#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
pub struct Material {
    pub color: [f64; 3],
    pub emission: [f64; 3],
	pub specular_probability: f64,
	pub maximum_specular_angle: f64,
	pub refractive_index: f64,
	pub is_opaque: bool,
}

impl Material {
    pub fn new(color: [f64; 3], emission: [f64; 3], specular_probability: f64, maximum_specular_angle: f64, refractive_index: f64, is_opaque: bool) -> Self {
        Self {
            color,
            emission,
            specular_probability,
            maximum_specular_angle,
            refractive_index,
            is_opaque,
        }
    }
}
