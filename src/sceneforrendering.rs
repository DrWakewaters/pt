use sceneforphysics::SceneForPhysics;
use spherefast::SphereFast;
use spherelightsource::SphereLightsource;
use trianglefast::TriangleFast;
use trianglelightsource::TriangleLightsource;

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct SceneForRendering {
	pub triangle_lightsources: Vec<TriangleLightsource>,
	pub triangle_surfaces: Vec<TriangleFast>,
	pub sphere_lightsources: Vec<SphereLightsource>,
	pub sphere_surfaces: Vec<SphereFast>,
}

impl SceneForRendering {
	pub fn new(scene_for_physics: SceneForPhysics) -> Self {
		let mut triangle_lightsources: Vec<TriangleLightsource> = Vec::new();
		let mut triangle_surfaces: Vec<TriangleFast> = Vec::new();
		let mut sphere_lightsources: Vec<SphereLightsource> = Vec::new();
		let mut sphere_surfaces: Vec<SphereFast> = Vec::new();
		for triangle in scene_for_physics.triangles {
			if triangle.is_lightsource {
				let lightsource_triangle_for_rendering = TriangleLightsource::new(&triangle);
				triangle_lightsources.push(lightsource_triangle_for_rendering);
			}
			if triangle.is_surface {
				let surface_triangle_for_rendering = TriangleFast::new(&triangle);
				triangle_surfaces.push(surface_triangle_for_rendering);
			}
		}
		for sphere in scene_for_physics.spheres {
			if sphere.is_lightsource {
				let lightsource_sphere_for_rendering = SphereLightsource::new(&sphere);
				sphere_lightsources.push(lightsource_sphere_for_rendering);
			}
			if sphere.is_surface {
				let surface_sphere_for_rendering = SphereFast::new(&sphere);
				sphere_surfaces.push(surface_sphere_for_rendering);
			}
		}
		Self {
			triangle_lightsources,
			triangle_surfaces,
			sphere_lightsources,
			sphere_surfaces,
		}
	}
}
