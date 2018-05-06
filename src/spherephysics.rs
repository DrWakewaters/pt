#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct SpherePhysics {
	pub center: [f64; 3],
	pub radius: f64,
	pub color: [f64; 3],
	pub color_alt: [f64; 3],
	pub local_x: [f64; 3],
	pub local_y: [f64; 3],
	pub local_z: [f64; 3],
	pub rotational_axis: [f64; 3],
	pub rotational_angular_speed: f64,
	pub specular_probability: f64,
	pub maximum_specular_angle: f64,
	pub refractive_index: f64,
	pub is_lightsource: bool,
	pub is_opaque: bool,
	pub velocity: [f64; 3],
	pub inertia: [[f64; 3]; 3],
	pub mass: f64,
	pub invisible_for_camera_ray: bool,
	pub is_surface: bool
}

impl SpherePhysics {
	pub fn new(center: [f64; 3], radius: f64, color: [f64; 3], color_alt: [f64; 3], local_x: [f64; 3], local_y: [f64; 3], local_z: [f64; 3], rotational_axis: [f64; 3], rotational_angular_speed: f64, specular_probability: f64, maximum_specular_angle: f64, refractive_index: f64, is_lightsource: bool, is_opaque: bool, velocity: [f64; 3], inertia: [[f64; 3]; 3], mass: f64, invisible_for_camera_ray: bool, is_surface: bool) -> Self {
		Self {
			center,
			radius,
			color,
			color_alt,
			local_x,
			local_y,
			local_z,
			rotational_axis,
			rotational_angular_speed,
			specular_probability,
			maximum_specular_angle,
			refractive_index,
			is_lightsource,
			is_opaque,
			velocity,
			inertia,
			mass,
			invisible_for_camera_ray,
			is_surface,
		}
	}
}
