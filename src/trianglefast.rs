use math::dot;
use trianglephysics::TrianglePhysics;

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

	pub fn compute_intersection_color(&self) -> [f64; 3] {
		if self.is_lightsource {
			return self.color;
		}
		self.color
	}
}
