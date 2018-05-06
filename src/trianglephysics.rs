use math::{cross, normalised, sub};

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
	pub is_surface: bool,
}

impl TrianglePhysics {
	pub fn new(nodes: &[[f64; 3]], indices: [usize; 3], color: [f64; 3], specular_probability: f64, maximum_specular_angle: f64, refractive_index: f64, is_lightsource: bool, is_opaque: bool, invisible_for_camera_ray: bool, is_surface: bool) -> Self {
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
			is_surface,
		}
	}
}
