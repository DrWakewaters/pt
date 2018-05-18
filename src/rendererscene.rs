use light::Light;
use physicsscene::PhysicsScene;
//use renderershape::RendererShape;
use renderersphere::RendererSphere;
use renderertriangle::RendererTriangle;

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct RendererScene {
	pub lights: Vec<Light>,
	pub renderer_spheres: Vec<RendererSphere>,
	pub renderer_triangles: Vec<RendererTriangle>,
}

impl RendererScene {
	pub fn new(physics_scene: PhysicsScene) -> Self {
		let mut lights: Vec<Light> = Vec::new();
		for light in physics_scene.lights {
			lights.push(light);
		}
		let mut renderer_spheres: Vec<RendererSphere> = Vec::new();
		for physics_sphere in physics_scene.physics_spheres {
			renderer_spheres.push(RendererSphere::new(&physics_sphere));
		}
		let mut renderer_triangles: Vec<RendererTriangle> = Vec::new();
		for physics_triangle in physics_scene.physics_triangles {
			renderer_triangles.push(RendererTriangle::new(&physics_triangle));
		}

		Self {
			lights,
			renderer_spheres,
			renderer_triangles,
		}
	}
}