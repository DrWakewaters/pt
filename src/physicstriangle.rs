use material::Material;
use math::{cross, normalised, sub};
use physics::Physics;

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct PhysicsTriangle {
	pub node0: [f64; 3],
	pub node1: [f64; 3],
	pub node2: [f64; 3],
	pub e1: [f64; 3],
	pub e2: [f64; 3],
	pub normal: [f64; 3],
	pub material: Material,
	pub physics: Physics,
}

impl PhysicsTriangle {
	pub fn new(nodes: &[[f64; 3]], indices: [usize; 3], material: Material, physics: Physics) -> Self {
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
			material,
			physics,
		}
	}
}
