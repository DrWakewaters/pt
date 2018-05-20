use pcg_rand::Pcg32;

use math::random_uniform_on_sphere;

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct Light {
    pub position: [f64; 3],
    pub color: [f64; 3],
}

impl Light {
    pub fn new(position: [f64; 3], color: [f64; 3]) -> Self {
        Self {
            position,
            color,
        }
    }

    pub fn direction(&self, pcg: &mut Pcg32) -> [f64; 3] {
        random_uniform_on_sphere(pcg)
    }
}
