use serde_derive::{Serialize, Deserialize};

#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
pub struct Material {
    pub color: [f64; 3],
	pub specular_probability: f64,
	pub maximum_specular_angle: f64,
	pub refractive_index: f64,
    pub is_opaque: bool,
    pub is_light: bool
}

impl Material {
    pub fn new(color: [f64; 3], specular_probability: f64, maximum_specular_angle: f64, refractive_index: f64, is_opaque: bool, is_light: bool) -> Self {
        Self {
            color,
            specular_probability,
            maximum_specular_angle,
            refractive_index,
            is_opaque,
            is_light,
        }
    }

    pub fn none() -> Self {
        Self{
            color: [0.0, 0.0, 0.0],
            specular_probability: 0.0,
            maximum_specular_angle: 0.0,
            refractive_index: 0.0,
            is_opaque: false,
            is_light: false,
        }
    }
}
