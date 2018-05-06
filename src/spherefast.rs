use spherephysics::SpherePhysics;
use math::{solve_linear, sub};

#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
pub struct SphereFast {
	pub center: [f64; 3],
	pub radius: f64,
	pub color: [f64; 3],
	pub color_alt: [f64; 3],
	pub local_x: [f64; 3],
	pub local_y: [f64; 3],
	pub local_z: [f64; 3],
	pub specular_probability: f64,
	pub maximum_specular_angle: f64,
	pub refractive_index: f64,
	pub is_lightsource: bool,
	pub invisible_for_camera_ray: bool,
	pub is_opaque: bool,
}

impl SphereFast {
	pub fn new(sphere_physics: &SpherePhysics) -> Self {
		Self {
			center: sphere_physics.center,
			radius: sphere_physics.radius,
			color: sphere_physics.color,
			color_alt: sphere_physics.color_alt,
			local_x: sphere_physics.local_x,
			local_y: sphere_physics.local_y,
			local_z: sphere_physics.local_z,
			specular_probability: sphere_physics.specular_probability,
			maximum_specular_angle: sphere_physics.maximum_specular_angle,
			refractive_index: sphere_physics.refractive_index,
			is_lightsource: sphere_physics.is_lightsource,
			invisible_for_camera_ray: sphere_physics.invisible_for_camera_ray,
			is_opaque: sphere_physics.is_opaque,
		}
	}

	#[allow(dead_code)]
	pub fn compute_intersection_color(&self) -> [f64; 3] {
		self.color
	}

	#[allow(dead_code)]
	pub fn compute_intersection_color_striped(&self, intersection: &mut[f64; 3]) -> [f64; 3] {
		if self.is_lightsource {
			return self.color;
		}
		// Compute the intersection relative to the center of the sphere.
		let intersection_translated = sub(*intersection, self.center);
		// Compute the intersection in the local coordinates of the sphere.
		let local_intersection = solve_linear([self.local_x, self.local_y, self.local_z], intersection_translated);
		self.get_color(local_intersection)
	}

	#[allow(dead_code)]
	pub fn get_color(&self, local_intersection: [f64; 3]) -> [f64; 3] {
		let region_y = ((local_intersection[1]/self.radius)*5.0).floor() as i64;
		match region_y%2 {
			0 => {
				self.color
			}
			_ => {
				self.color_alt
			}
		}
	}
}
