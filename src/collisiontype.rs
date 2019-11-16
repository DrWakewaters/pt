use serde_derive::{Serialize, Deserialize};

#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
pub enum CollisionType {
	Sphere,
	Triangle,
}
