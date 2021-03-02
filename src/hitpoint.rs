use crate::material::Material;

pub struct Hitpoint {
    pub position: [f64; 3],
    pub incoming_direction: [f64; 3],
    pub normal: [f64; 3],
    pub material: Material,
    pub hit_from_outside: bool,
    pub transmitted: bool,
}

impl Hitpoint {
    pub fn new(position: [f64; 3], incoming_direction: [f64; 3], normal: [f64; 3], material: Material, hit_from_outside: bool, transmitted: bool) -> Self {
        Self {
            position,
            incoming_direction,
            normal,
            material,
            hit_from_outside,
            transmitted,
        }
    }
}
