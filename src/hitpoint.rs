use material::Material;

pub struct Hitpoint {
    pub position: [f64; 3],
    pub distance: f64,
    pub normal: [f64; 3],
    pub material: Material,
    pub hit_from_outside: bool
}

impl Hitpoint {
    pub fn new(position: [f64; 3], distance: f64, normal: [f64; 3], material: Material, hit_from_outside: bool) -> Self {
        Self {
            position,
            distance,
            normal,
            material,
            hit_from_outside,
        }
    }
}
