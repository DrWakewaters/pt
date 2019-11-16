use serde_derive::{Serialize, Deserialize};

#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
pub struct Physics {
    pub local_x: [f64; 3],
    pub local_y: [f64; 3],
    pub local_z: [f64; 3],
    pub rotational_axis: [f64; 3],
    pub rotational_angular_speed: f64,
    pub velocity: [f64; 3],
    pub inertia: [[f64; 3]; 3],
    pub mass: f64,
    pub is_movable: bool
}

impl Physics {
    pub fn new(local_x: [f64; 3], local_y: [f64; 3], local_z: [f64; 3], rotational_axis: [f64; 3], rotational_angular_speed: f64, velocity: [f64; 3], inertia: [[f64; 3]; 3], mass: f64, is_movable: bool) -> Self {
        Self {
            local_x,
            local_y,
            local_z,
            rotational_axis,
            rotational_angular_speed,
            velocity,
            inertia,
            mass,
            is_movable,
        }
    }
}
