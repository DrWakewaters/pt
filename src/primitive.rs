use rand::Rng;
use pcg_rand::Pcg32;

use math::{add, cross, dot, mul, normalised, pick_reflection_lambertian, pick_reflection_uniform, solve_linear, sub};

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct TrianglePhysics {
	pub node0: [f64; 3],
	pub node1: [f64; 3],
	pub node2: [f64; 3],
	pub e1: [f64; 3],
	pub e2: [f64; 3],
	pub normal: [f64; 3],
	pub color: [f64; 3],
	pub specular_probability: f64,
	pub maximum_specular_angle: f64,
	pub refractive_index: f64,
	pub is_lightsource: bool,
	pub is_opaque: bool,
	pub invisible_for_camera_ray: bool,
}

impl TrianglePhysics {
	pub fn new(nodes: &Vec<[f64; 3]>, indices: [usize; 3], color: [f64; 3], specular_probability: f64, maximum_specular_angle: f64, refractive_index: f64, is_lightsource: bool, is_opaque: bool, invisible_for_camera_ray: bool) -> Self {
		let e1 = sub(nodes[indices[1]], nodes[indices[0]]);
		let e2 = sub(nodes[indices[2]], nodes[indices[0]]);
		let normal = normalised(cross(sub(nodes[indices[1]], nodes[indices[0]]), sub(nodes[indices[2]], nodes[indices[0]])));
		Self {
			node0: nodes[indices[0]],
			node1: nodes[indices[1]],
			node2: nodes[indices[2]],
			e1,
			e2,
			normal,
			color,
			specular_probability,
			maximum_specular_angle,
			refractive_index,
			is_lightsource,
			is_opaque,
			invisible_for_camera_ray,
		}
	}
}

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
}

impl SpherePhysics {
	pub fn new(center: [f64; 3], radius: f64, color: [f64; 3], color_alt: [f64; 3], local_x: [f64; 3], local_y: [f64; 3], local_z: [f64; 3], rotational_axis: [f64; 3], rotational_angular_speed: f64, specular_probability: f64, maximum_specular_angle: f64, refractive_index: f64, is_lightsource: bool, is_opaque: bool, velocity: [f64; 3], inertia: [[f64; 3]; 3], mass: f64, invisible_for_camera_ray: bool) -> Self {
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
		}
	}
}

#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
pub struct TriangleLightsource {
	pub node0: [f64; 3],
	pub node1: [f64; 3],
	pub node2: [f64; 3],
	pub e1: [f64; 3],
	pub e2: [f64; 3],
	pub normal: [f64; 3],
	pub color: [f64; 3],
}

impl TriangleLightsource {
	pub fn new(triangle_physics: &TrianglePhysics) -> Self {
		Self {
			node0: triangle_physics.node0,
			node1: triangle_physics.node1,
			node2: triangle_physics.node2,
			e1: triangle_physics.e1,
			e2: triangle_physics.e2,
			normal: triangle_physics.normal,
			color: triangle_physics.color,
		}
	}
}

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
	pub fn compute_ray_data(&self, mut pcg: &mut Pcg32) -> ([f64; 3], [f64; 3], [f64; 3], [f64; 3]) {
		let mut mock_normal = [0.0, 1.0, 0.0];
		let rand: u32 = pcg.gen();
		if (rand as f64)/(<u32>::max_value() as f64) > 0.5 {
			mock_normal = mul(-1.0, mock_normal);
		}
		//let mut point = pick_reflection_uniform_on_stripe(mock_normal, &mut pcg, PI/6.0, 2.0*PI/6.0);
		let mut point = pick_reflection_uniform(mock_normal, &mut pcg);
		point = add(mul(self.radius, point), self.center);
		let normal = normalised(sub(point, self.center));
		let direction = pick_reflection_lambertian(normal, &mut pcg);
		let pixel = self.color;
		(point, direction, normal, pixel)
	}
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct TriangleFast {
	pub n_u: f64,
	pub n_v: f64,
	pub n_d: f64,
	pub k: u64,
	pub b_nu: f64,
	pub b_nv: f64,
	pub b_d: f64,
	pub specular_probability: f64,
	pub c_nu: f64,
	pub c_nv: f64,
	pub c_d: f64,
	pub maximum_specular_angle: f64,
	pub refractive_index: f64,
	pub color: [f64; 3],
	pub normal: [f64; 3],
	pub invisible_for_camera_ray: bool,
	pub is_lightsource: bool,
	pub is_opaque: bool,
}

impl TriangleFast {
	pub fn new(triangle_physics: &TrianglePhysics) -> Self {
		let b = triangle_physics.e2;
		let c = triangle_physics.e1;
		let node0 = triangle_physics.node0;
		let normal = triangle_physics.normal;
		let mut k = 0;
		let mut u = 1;
		let mut v = 2;
		if normal[1].abs() > normal[0].abs() && normal[1].abs() > normal[2].abs() {
			k = 1;
			u = 2;
			v = 0;
		}
		if normal[2].abs() > normal[0].abs() && normal[2].abs() > normal[1].abs() {
			k = 2;
			u = 0;
			v = 1;
		}
		let denom = b[u]*c[v]-b[v]*c[u];
		Self {
			n_u: normal[u]/normal[k],
			n_v: normal[v]/normal[k],
			n_d: dot(normal, node0)/normal[k],
			k: k as u64,
			b_nu: b[u]/denom,
			b_nv: -1.0*b[v]/denom,
			b_d: (b[v]*node0[u]-b[u]*node0[v])/denom,
			specular_probability: triangle_physics.specular_probability,
			c_nu: -1.0*c[u]/denom,
			c_nv: c[v]/denom,
			c_d: (c[u]*node0[v]-c[v]*node0[u])/denom,
			maximum_specular_angle: triangle_physics.maximum_specular_angle,
			refractive_index: triangle_physics.refractive_index,
			color: triangle_physics.color,
			normal: triangle_physics.normal,
			invisible_for_camera_ray: triangle_physics.invisible_for_camera_ray,
			is_lightsource: triangle_physics.is_lightsource,
			is_opaque: triangle_physics.is_opaque,
		}
	}

	pub fn compute_intersection_color(&self, intersection: &mut[f64; 3]) -> [f64; 3] {
		if self.is_lightsource {
			return self.color;
		}
		self.color
	}
}

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

	pub fn compute_intersection_color(&self, intersection: &mut[f64; 3]) -> [f64; 3] {
		return self.color;

		if self.is_lightsource {
			return self.color;
		}
		// Compute the intersection relative to the center of the sphere.
		let intersection_translated = sub(*intersection, self.center);
		// Compute the intersection in the local coordinates of the sphere.
		let local_intersection = solve_linear([self.local_x, self.local_y, self.local_z], intersection_translated);
		self.get_color(local_intersection)
	}

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
