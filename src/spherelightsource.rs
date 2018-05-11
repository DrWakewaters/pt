use pcg_rand::Pcg32;
use rand::Rng;

use math::{add, mul, normalised, pick_reflection_lambertian, pick_reflection_uniform, sub};
use ray::Ray;
use spherephysics::SpherePhysics;

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct SphereLightsource {
	pub center: [f64; 3],
	pub radius: f64,
	pub color: [f64; 3],
	pub color_alt: [f64; 3],
	pub local_x: [f64; 3],
	pub local_y: [f64; 3],
	pub local_z: [f64; 3],
	pub rotational_axis: [f64; 3],
	pub rotational_angular_speed: f64,
}

impl SphereLightsource {
	pub fn new(sphere_physics: &SpherePhysics) -> Self {
		Self {
			center: sphere_physics.center,
			radius: sphere_physics.radius,
			color: sphere_physics.color,
			color_alt: sphere_physics.color_alt,
			local_x: sphere_physics.local_x,
			local_y: sphere_physics.local_y,
			local_z: sphere_physics.local_z,
			rotational_axis: sphere_physics.rotational_axis,
			rotational_angular_speed: sphere_physics.rotational_angular_speed,
		}
	}

	// @TODO Implement non-lambertian emission.
	pub fn compute_ray_data(&self, mut pcg: &mut Pcg32) -> Ray {
		let mut mock_normal = [0.0, 1.0, 0.0];
		//@TODO: What does this do?
		if pcg.next_f64() > 0.5 {
			mock_normal = mul(-1.0, mock_normal);
		}
		//let mut point = pick_reflection_uniform_on_stripe(mock_normal, &mut pcg, PI/6.0, 2.0*PI/6.0);
		let mut point = pick_reflection_uniform(mock_normal, &mut pcg);
		point = add(mul(self.radius, point), self.center);
		let normal = normalised(sub(point, self.center));
		let direction = pick_reflection_lambertian(normal, &mut pcg);
		let old_direction = direction;
		let color = self.color;
		Ray::new(point, direction, old_direction, normal, color, 0.0, 1.0, 0.0, 1.0)
	}
}
