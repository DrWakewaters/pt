use serde_derive::{Serialize, Deserialize};

use crate::camera::Camera;
use crate::lightsphere::LightSphere;
use crate::physicsscene::PhysicsScene;
use crate::renderersphere::RendererSphere;
use crate::renderertriangle::RendererTriangle;

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct RendererScene {
	pub light_spheres: Vec<LightSphere>,
	pub renderer_spheres: Vec<RendererSphere>,
	pub renderer_triangles: Vec<RendererTriangle>,
	pub cameras: Vec<Camera>
}

impl RendererScene {
	pub fn new(physics_scene: PhysicsScene) -> Self {
		let mut light_spheres: Vec<LightSphere> = Vec::new();
		for light_sphere in physics_scene.light_spheres {
			light_spheres.push(light_sphere);
		}
		let mut renderer_spheres: Vec<RendererSphere> = Vec::new();
		for physics_sphere in physics_scene.physics_spheres {
			renderer_spheres.push(RendererSphere::new(&physics_sphere));
		}
		let mut renderer_triangles: Vec<RendererTriangle> = Vec::new();
		for physics_triangle in physics_scene.physics_triangles {
			renderer_triangles.push(RendererTriangle::new(&physics_triangle));
		}
		let mut cameras: Vec<Camera> = Vec::new();
		for camera in physics_scene.cameras {
			cameras.push(camera);
		}
		Self {
			light_spheres,
			renderer_spheres,
			renderer_triangles,
			cameras,
		}
	}
}
