use material::Material;

pub struct Hitpoint {
    pub position: [f64; 3],
    pub incoming_direction: [f64; 3],
    pub distance: f64,
    pub normal: [f64; 3],
    pub material: Material,
    pub hit_from_outside: bool,
    pub is_end_point: bool,
    pub accumulated_color: [f64; 3],
    pub distance_from_retina: f64,
}

impl Hitpoint {
    pub fn new(position: [f64; 3], incoming_direction: [f64; 3], distance: f64, normal: [f64; 3], material: Material, hit_from_outside: bool, is_end_point: bool, accumulated_color: [f64; 3], distance_from_retina: f64) -> Self {
        Self {
            position,
            incoming_direction,
            distance,
            normal,
            material,
            hit_from_outside,
            is_end_point,
            accumulated_color,
            distance_from_retina,
        }
    }
}
