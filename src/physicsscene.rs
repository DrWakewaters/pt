use serde_derive::{Serialize, Deserialize};

use crate::camera::Camera;
use crate::collisiontype::CollisionType;
use crate::math::{add, distance_squared, dot, mul, normalised, quarternion_inverse, quarternion_product, same_side, sub};
use crate::lightsphere::LightSphere;
use crate::material::Material;
use crate::physics::Physics;
use crate::physicscylinder::PhysicsCylinder;
use crate::physicssphere::PhysicsSphere;
use crate::physicstriangle::PhysicsTriangle;

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct PhysicsScene {
	pub time: f64,
	pub time_per_frame: f64,
	pub physics_cylinders: Vec<PhysicsCylinder>,
	pub physics_triangles: Vec<PhysicsTriangle>,
	pub physics_spheres: Vec<PhysicsSphere>,
	pub light_spheres: Vec<LightSphere>,
	pub cameras: Vec<Camera>,
}

impl PhysicsScene {
	pub fn new(time_per_frame: f64) -> Self {
		// This scene is modelled after https://en.wikipedia.org/wiki/Path_tracing#/media/File:Path_tracing_001.png.

		let mut id = 0;

		let white_diffuse_opaque_emissive = Material::new([0.8, 0.8, 0.8], mul(0.2, [0.8, 0.8, 0.8]), 0.0, 1.0, 0.0, true);
		let white_diffuse_opaque = Material::new([0.8, 0.8, 0.8], mul(0.0, [0.8, 0.8, 0.8]), 0.0, 0.03, 0.0, true);
		let white_semi_specular_opaque = Material::new([0.8, 0.8, 0.8], mul(0.0, [0.8, 0.8, 0.8]), 0.1, 0.03, 0.0, true);
		let white_specular_opaque = Material::new([0.8, 0.8, 0.8], mul(0.0, [0.8, 0.8, 0.8]), 1.0, 0.03, 0.0, true);

		let white_diffuse_transparent = Material::new([1.0, 1.0, 1.0], mul(0.0, [1.0, 1.0, 1.0]), 0.0, 0.01, 1.4, false);
		let white_semi_specular_transparent_emissive = Material::new([1.0, 0.96, 0.9], mul(0.0, [1.0, 0.96, 0.9]), 0.3, 0.03, 1.4, false);
		let white_semi_specular_transparent = Material::new([1.0, 1.0, 1.0], mul(0.0, [1.0, 1.0, 1.0]), 0.3, 0.03, 1.4, false);
		let white_specular_transparent = Material::new([1.0, 1.0, 1.0], mul(0.0, [1.0, 1.0, 1.0]), 1.0, 0.03, 1.4, false);

		let green_diffuse_opaque = Material::new([0.15, 0.65, 0.15], mul(0.0, [0.15, 0.65, 0.15]), 0.0, 0.03, 0.0, true);
		let red_diffuse_opaque = Material::new([0.65, 0.15, 0.15], mul(0.0, [0.65, 0.15, 0.15]), 0.0, 0.03, 0.0, true);
		let blue_diffuse_opaque = Material::new([0.15, 0.15, 0.65], mul(0.0, [0.15, 0.15, 0.65]), 0.0, 0.03, 0.0, true);
		let dark_red_diffuse_opaque_emissive = Material::new([0.6, 0.45, 0.45], mul(0.015, [0.6, 0.45, 0.45]), 0.0, 0.03, 0.0, true);
		let dark_blue_diffuse_opaque_emissive = Material::new([0.45, 0.45, 0.6], mul(0.015, [0.45, 0.45, 0.6]), 0.0, 0.03, 0.0, true);
		let light_red_diffuse_transparent_emissive = Material::new([0.95, 0.9, 0.9], mul(0.4, [0.95, 0.9, 0.9]), 0.0, 0.03, 0.0, false);
		let light_blue_diffuse_transparent_emissive = Material::new([0.9, 0.9, 0.95], mul(0.4, [0.9, 0.9, 0.95]), 0.0, 0.03, 0.0, false);

		let local_x = [1.0, 0.0, 0.0];
		let local_y = [0.0, 1.0, 0.0];
		let local_z = [0.0, 0.0, 1.0];

		let light_position = [500.0, 500.0, -2500.0];

		let light_radius = 1200.0;
		let large_radius = 160.0;
		let middle_radius = 75.0;
		let small_radius = 70.0;
		let mass = 1000.0;

		let inertia_diagonal_light_radius = 2.0/5.0*mass*light_radius*light_radius;
		let inertia_light_radius = [[inertia_diagonal_light_radius, 0.0, 0.0], [0.0, inertia_diagonal_light_radius, 0.0], [0.0, 0.0, inertia_diagonal_light_radius]];
		let inertia_diagonal_large_radius = 2.0/5.0*mass*large_radius*large_radius;
		let inertia_large_radius = [[inertia_diagonal_large_radius, 0.0, 0.0], [0.0, inertia_diagonal_large_radius, 0.0], [0.0, 0.0, inertia_diagonal_large_radius]];
		let inertia_diagonal_middle_radius = 2.0/5.0*mass*middle_radius*middle_radius;
		let inertia_middle_radius = [[inertia_diagonal_middle_radius, 0.0, 0.0], [0.0, inertia_diagonal_middle_radius, 0.0], [0.0, 0.0, inertia_diagonal_middle_radius]];
		let inertia_diagonal_small_radius = 2.0/5.0*mass*small_radius*small_radius;
		let inertia_small_radius = [[inertia_diagonal_small_radius, 0.0, 0.0], [0.0, inertia_diagonal_small_radius, 0.0], [0.0, 0.0, inertia_diagonal_small_radius]];

		let physics_light_radius = Physics::new(local_x, local_y, local_z, [1.0, 0.0, 0.0], 0.0, [0.0, 0.0, 0.0], inertia_light_radius, mass, true);
		let physics_large_radius = Physics::new(local_x, local_y, local_z, [1.0, 0.0, 0.0], 0.0, [0.0, 0.0, 0.0], inertia_large_radius, mass, true);
		let physics_middle_radius = Physics::new(local_x, local_y, local_z, [1.0, 0.0, 0.0], 0.0, [0.0, 0.0, 0.0], inertia_middle_radius, mass, true);
		let physics_small_radius = Physics::new(local_x, local_y, local_z, [1.0, 0.0, 0.0], 0.0, [0.0, 0.0, 0.0], inertia_small_radius, mass, true);

		let physics_triangle = Physics::new(local_x, local_y, local_z, [1.0, 0.0, 0.0], 0.0, [0.0, 0.0, 0.0], [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]], 0.0, false);

		// @TODO: Implement cylinder physics.
		let physics_cylinder = Physics::new(local_x, local_y, local_z, [1.0, 0.0, 0.0], 0.0, [0.0, 0.0, 0.0], [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]], 0.0, false);

		let handle_length: f64 = (200.0*200.0+100.0*100.0 as f64).sqrt()*0.999;
		let physics_cylinders: Vec<PhysicsCylinder> = vec![
			PhysicsCylinder::new([100.0, 1350.0, 900.0], [0.5, -1.0, 0.0], handle_length, 24.0, dark_red_diffuse_opaque_emissive, physics_cylinder, 0, false),
			PhysicsCylinder::new([900.0, 1350.0, 865.0], [-0.5, -1.0, 0.0], handle_length, 24.0, dark_blue_diffuse_opaque_emissive, physics_cylinder, 0, false),
			PhysicsCylinder::new([200.0, 1150.0, 900.0], [0.5, -1.0, 0.0], 0.0, 15.0, light_red_diffuse_transparent_emissive, physics_cylinder, 1, false),
			PhysicsCylinder::new([800.0, 1150.0, 865.0], [-0.5, -1.0, 0.0], 0.0, 15.0, light_blue_diffuse_transparent_emissive, physics_cylinder, 1, false)
		];

		// x increases to the right, y downwards and z inwards.
		let nodes = vec![
		// Left wall
		[-600.0, 1000.0, -1500.0],		// 0
		[-600.0, 500.0, -1500.0],		// 1
		[-600.0, 500.0, 500.0],			// 2
		[-600.0, 1000.0, 500.0],		// 3
		// Right wall
		[1600.0, 1000.0, 500.0],		// 4
		[1600.0, 500.0, 500.0],			// 5
		[1600.0, 500.0, -1500.0],		// 6
		[1600.0, 1000.0, -1500.0],		// 7
		// Far wall
		/*
		[-600.0, 1000.0, 1000.0],		// 8
		[-600.0, -600.0, 1000.0],		// 9
		[1600.0, -600.0, 1000.0],		// 10
		[1600.0, 1000.0, 1000.0],		// 11
		*/
		[-600.0, 1500.0, 1000.0],		// 8
		[-600.0, -600.0, 1000.0],		// 9
		[1600.0, -600.0, 1000.0],		// 10
		[1600.0, 1500.0, 1000.0],		// 11
		// Floor
		[-600.0, 1000.0, -2500.0],		// 12
		[-600.0, 1000.0, 1000.0],		// 13
		[1600.0, 1000.0, 1000.0],		// 14
		[1600.0, 1000.0, -2500.0],		// 15
		// Top
		[0.0, -1000.0, 0.0],			// 16
		[1000.0, -1000.0, 0.0],			// 17
		[1000.0, -1000.0, 1000.0],		// 18
		[0.0, -1000.0, 1000.0],			// 19
		];

		let physics_spheres: Vec<PhysicsSphere> = vec![
			/*
			PhysicsSphere::new([0.0, 1000.0-large_radius, 600.0], large_radius, white_diffuse_opaque, physics_large_radius),
			PhysicsSphere::new([500.0, 1000.0-large_radius, 600.0], large_radius, white_semi_specular_opaque, physics_large_radius),
			PhysicsSphere::new([1000.0, 1000.0-large_radius, 600.0], large_radius, white_specular_opaque, physics_large_radius),
			PhysicsSphere::new([0.0, 1000.0-small_radius, 370.0], small_radius, red_diffuse_opaque, physics_small_radius),
			PhysicsSphere::new([250.0, 1000.0-small_radius, 370.0], small_radius, white_semi_specular_transparent, physics_small_radius),
			PhysicsSphere::new([500.0, 1000.0-small_radius, 370.0], small_radius, blue_diffuse_opaque, physics_small_radius),
			PhysicsSphere::new([750.0, 1000.0-small_radius, 370.0], small_radius, white_specular_transparent, physics_small_radius),
			PhysicsSphere::new([1000.0, 1000.0-small_radius, 370.0], small_radius, green_diffuse_opaque, physics_small_radius),
			PhysicsSphere::new([-150.0, 200.0, 370.0], middle_radius, white_specular_transparent, physics_middle_radius),
			PhysicsSphere::new([500.0, 200.0, 370.0], middle_radius, white_specular_transparent, physics_middle_radius),
			PhysicsSphere::new([1150.0, 200.0, 370.0], middle_radius, white_specular_transparent, physics_middle_radius),
			*/
			//PhysicsSphere::new(light_position, light_radius, white_diffuse_opaque_emissive, physics_light_radius)
		];

		let wall_triangle_indices = vec![
			// Left wall.
			//[3, 1, 0],
			//[3, 2, 1],
			// Right wall.
			//[7, 5, 4],
			//[7, 6, 5],
			// Far wall.
			//[11, 9, 8],
			//[11, 10, 9],
			// Floor.
			//[15, 13, 12],
			//[15, 14, 13],
			// Top.
			//[19, 17, 16],
			//[19, 18, 17]
		];
		let mut physics_triangles_walls: Vec<PhysicsTriangle> = vec![];
		for wall_triangle_index in wall_triangle_indices {
			physics_triangles_walls.push(PhysicsTriangle::new(&nodes, wall_triangle_index, white_diffuse_opaque, physics_triangle, id, false));
		}
		id += 1;

		let letter_c_thickness = 75.0;
		let letter_c_dimensions = [300.0, 300.0, 100.0]; // W x H x D
		let mut letter_c_nodes = vec![
			// Upper layer.
			// Zeroth row
			[0.0, 0.0, 0.0],
			[letter_c_dimensions[0], 0.0, 0.0],
			// First row.
			[0.0, letter_c_thickness, 0.0],
			[letter_c_thickness, letter_c_thickness, 0.0],
			[letter_c_dimensions[0], letter_c_thickness, 0.0],
			// Second row.
			[0.0, letter_c_dimensions[1]-letter_c_thickness, 0.0],
			[letter_c_thickness, letter_c_dimensions[1]-letter_c_thickness, 0.0],
			[letter_c_dimensions[0], letter_c_dimensions[1]-letter_c_thickness, 0.0],
			// Third row.
			[0.0, letter_c_dimensions[1], 0.0],
			[letter_c_dimensions[0], letter_c_dimensions[1], 0.0],
			// Lower layer.
			// Zeroth row
			[0.0, 0.0, letter_c_dimensions[2]],
			[letter_c_dimensions[0], 0.0, letter_c_dimensions[2]],
			// First row.
			[0.0, letter_c_thickness, letter_c_dimensions[2]],
			[letter_c_thickness, letter_c_thickness, letter_c_dimensions[2]],
			[letter_c_dimensions[0], letter_c_thickness, letter_c_dimensions[2]],
			// Second row.
			[0.0, letter_c_dimensions[1]-letter_c_thickness, letter_c_dimensions[2]],
			[letter_c_thickness, letter_c_dimensions[1]-letter_c_thickness, letter_c_dimensions[2]],
			[letter_c_dimensions[0], letter_c_dimensions[1]-letter_c_thickness, letter_c_dimensions[2]],
			// Third row.
			[0.0, letter_c_dimensions[1], letter_c_dimensions[2]],
			[letter_c_dimensions[0], letter_c_dimensions[1], letter_c_dimensions[2]]
		];
		let letter_c_triangle_indices = vec![
			// Upper layer
			[0, 4, 1],
			[0, 2, 4],
			[2, 6, 3],
			[2, 5, 6],
			[5, 9, 7],
			[5, 8, 9],
			// Lower layer
			[10, 11, 14],
			[10, 14, 12],
			[12, 13, 16],
			[12, 16, 15],
			[15, 17, 19],
			[15, 19, 18],
			// Vertical sides, seen from left.
			[10, 8, 0],
			[10, 18, 8],
			// Vertical sides, seen from right.
			[1, 14, 11],
			[1, 4, 14],
			[3, 16, 13],
			[3, 6, 16],	
			[7, 19, 17],
			[7, 9, 19],
			// Horisontal sides, seen from the top.
			[0, 1, 11],
			[0, 11, 10],
			[6, 7, 17],
			[6, 17, 16],
			// Horisontal sides, seen from the bottom.
			[3, 13, 14],
			[3, 14, 4],
			[8, 18, 19],
			[8, 19, 9]	
		];
		let letter_c_translation = [-400.0, -200.0, 899.0];
		for node in letter_c_nodes.iter_mut() {
			*node = add(letter_c_translation, *node);
		}
		let mut physics_triangles_letter_c: Vec<PhysicsTriangle> = vec![];
		for letter_c_triangle_index in letter_c_triangle_indices {
			physics_triangles_letter_c.push(PhysicsTriangle::new(&letter_c_nodes, letter_c_triangle_index, white_semi_specular_transparent_emissive, physics_triangle, id, false));
		}
		id += 1;

		let letter_e_thickness = 75.0;
		let letter_e_dimensions = [300.0, 300.0, 100.0]; // W x H x D
		let mut letter_e_nodes = vec![
			// Upper layer.
			// Zeroth row
			[0.0, 0.0, 0.0],
			[letter_e_dimensions[0], 0.0, 0.0],
			// First row.
			[0.0, letter_e_thickness, 0.0],
			[letter_e_thickness, letter_e_thickness, 0.0],
			[letter_e_dimensions[0], letter_e_thickness, 0.0],
			// Second row.
			[0.0, letter_e_dimensions[1]/2.0-letter_e_thickness/2.0, 0.0],
			[letter_e_thickness, letter_e_dimensions[1]/2.0-letter_e_thickness/2.0, 0.0],
			[letter_e_dimensions[0], letter_e_dimensions[1]/2.0-letter_e_thickness/2.0, 0.0],
			// Thrid row.
			[0.0, letter_e_dimensions[1]/2.0+letter_e_thickness/2.0, 0.0],
			[letter_e_thickness, letter_e_dimensions[1]/2.0+letter_e_thickness/2.0, 0.0],
			[letter_e_dimensions[0], letter_e_dimensions[1]/2.0+letter_e_thickness/2.0, 0.0],
			// Fourth row.
			[0.0, letter_e_dimensions[1]-letter_e_thickness, 0.0],
			[letter_e_thickness, letter_e_dimensions[1]-letter_e_thickness, 0.0],
			[letter_e_dimensions[0], letter_e_dimensions[1]-letter_e_thickness, 0.0],
			// Fifth row.
			[0.0, letter_e_dimensions[1], 0.0],
			[letter_e_dimensions[0], letter_e_dimensions[1], 0.0],
			// Lower layer.
			// Zeroth row
			[0.0, 0.0, letter_e_dimensions[2]],
			[letter_e_dimensions[0], 0.0, letter_e_dimensions[2]],
			// First row.
			[0.0, letter_e_thickness, letter_e_dimensions[2]],
			[letter_e_thickness, letter_e_thickness, letter_e_dimensions[2]],
			[letter_e_dimensions[0], letter_e_thickness, letter_e_dimensions[2]],
			// Second row.
			[0.0, letter_e_dimensions[1]/2.0-letter_e_thickness/2.0, letter_e_dimensions[2]],
			[letter_e_thickness, letter_e_dimensions[1]/2.0-letter_e_thickness/2.0, letter_e_dimensions[2]],
			[letter_e_dimensions[0], letter_e_dimensions[1]/2.0-letter_e_thickness/2.0, letter_e_dimensions[2]],
			// Thrid row.
			[0.0, letter_e_dimensions[1]/2.0+letter_e_thickness/2.0, letter_e_dimensions[2]],
			[letter_e_thickness, letter_e_dimensions[1]/2.0+letter_e_thickness/2.0, letter_e_dimensions[2]],
			[letter_e_dimensions[0], letter_e_dimensions[1]/2.0+letter_e_thickness/2.0, letter_e_dimensions[2]],
			// Fourth row.
			[0.0, letter_e_dimensions[1]-letter_e_thickness, letter_e_dimensions[2]],
			[letter_e_thickness, letter_e_dimensions[1]-letter_e_thickness, letter_e_dimensions[2]],
			[letter_e_dimensions[0], letter_e_dimensions[1]-letter_e_thickness, letter_e_dimensions[2]],
			// Fifth row.
			[0.0, letter_e_dimensions[1], letter_e_dimensions[2]],
			[letter_e_dimensions[0], letter_e_dimensions[1], letter_e_dimensions[2]]
		];
		let letter_e_triangle_indices = vec![
			// Upper layer
			[0, 4, 1],
			[0, 2, 4],
			[2, 12, 3],
			[2, 11, 12],
			[11, 15, 13],
			[11, 14, 15],	
			[6, 10, 7],
			[6, 9, 10],
			// Lower layer
			[16, 17, 20],
			[16, 20, 18],
			[18, 19, 28],
			[18, 28, 27],
			[27, 29, 31],
			[27, 31, 30],	
			[22, 23, 26],
			[22, 26, 25],
			// Vertical sides, seen from left.
			[16, 14, 0],
			[16, 30, 14],
			// Vertical sides, seen from right.
			[1, 20, 17],
			[1, 4, 20],
			[3, 22, 19],
			[3, 6, 22],
			[7, 26, 23],
			[7, 10, 26],
			[9, 28, 25],
			[9, 12, 28],
			[13, 31, 29],
			[13, 15, 31],
			// Horisontal sides, seen from the top.
			[0, 1, 17],
			[0, 17, 16],
			[6, 7, 23],
			[6, 23, 22],
			[12, 13, 29],
			[12, 29, 28],
			// Horisontal sides, seen from the bottom.
			[3, 19, 20],
			[3, 20, 4],
			[9, 25, 26],
			[9, 26, 10],
			[14, 30, 31],
			[14, 31, 15]
		];
		let letter_e_translation = [-25.0, -200.0, 899.0];
		for node in letter_e_nodes.iter_mut() {
			*node = add(letter_e_translation, *node);
		}
		let mut physics_triangles_letter_e: Vec<PhysicsTriangle> = vec![];
		for letter_e_triangle_index in letter_e_triangle_indices {
			physics_triangles_letter_e.push(PhysicsTriangle::new(&letter_e_nodes, letter_e_triangle_index, white_semi_specular_transparent_emissive, physics_triangle, id, false));
		}
		id += 1;

		let letter_s_thickness = 75.0;
		let letter_s_dimensions = [300.0, 300.0, 100.0]; // W x H x D
		let mut letter_s_nodes = vec![
			// Upper layer.
			// Zeroth row
			[0.0, 0.0, 0.0],
			[letter_s_dimensions[0], 0.0, 0.0],
			// First row.
			[0.0, letter_s_thickness, 0.0],
			[letter_s_thickness, letter_s_thickness, 0.0],
			[letter_s_dimensions[0], letter_s_thickness, 0.0],
			// Second row.
			[0.0, letter_s_dimensions[1]/2.0-letter_s_thickness/2.0, 0.0],
			[letter_s_thickness, letter_s_dimensions[1]/2.0-letter_s_thickness/2.0, 0.0],
			[letter_s_dimensions[0], letter_s_dimensions[1]/2.0-letter_s_thickness/2.0, 0.0],
			// Thrid row.
			[0.0, letter_s_dimensions[1]/2.0+letter_s_thickness/2.0, 0.0],
			[letter_s_dimensions[0]-letter_s_thickness, letter_s_dimensions[1]/2.0+letter_s_thickness/2.0, 0.0],
			[letter_s_dimensions[0], letter_s_dimensions[1]/2.0+letter_s_thickness/2.0, 0.0],
			// Fourth row.
			[0.0, letter_s_dimensions[1]-letter_s_thickness, 0.0],
			[letter_s_dimensions[0]-letter_s_thickness, letter_s_dimensions[1]-letter_s_thickness, 0.0],
			[letter_s_dimensions[0], letter_s_dimensions[1]-letter_s_thickness, 0.0],
			// Fifth row.
			[0.0, letter_s_dimensions[1], 0.0],
			[letter_s_dimensions[0], letter_s_dimensions[1], 0.0],
			// Lower layer.
			// Zeroth row
			[0.0, 0.0, letter_s_dimensions[2]],
			[letter_s_dimensions[0], 0.0, letter_s_dimensions[2]],
			// First row.
			[0.0, letter_s_thickness, letter_s_dimensions[2]],
			[letter_s_thickness, letter_s_thickness, letter_s_dimensions[2]],
			[letter_s_dimensions[0], letter_s_thickness, letter_s_dimensions[2]],
			// Second row.
			[0.0, letter_s_dimensions[1]/2.0-letter_s_thickness/2.0, letter_s_dimensions[2]],
			[letter_s_thickness, letter_s_dimensions[1]/2.0-letter_s_thickness/2.0, letter_s_dimensions[2]],
			[letter_s_dimensions[0], letter_s_dimensions[1]/2.0-letter_s_thickness/2.0, letter_s_dimensions[2]],
			// Thrid row.
			[0.0, letter_s_dimensions[1]/2.0+letter_s_thickness/2.0, letter_s_dimensions[2]],
			[letter_s_dimensions[0]-letter_s_thickness, letter_s_dimensions[1]/2.0+letter_s_thickness/2.0, letter_s_dimensions[2]],
			[letter_s_dimensions[0], letter_s_dimensions[1]/2.0+letter_s_thickness/2.0, letter_s_dimensions[2]],
			// Fourth row.
			[0.0, letter_s_dimensions[1]-letter_s_thickness, letter_s_dimensions[2]],
			[letter_s_dimensions[0]-letter_s_thickness, letter_s_dimensions[1]-letter_s_thickness, letter_s_dimensions[2]],
			[letter_s_dimensions[0], letter_s_dimensions[1]-letter_s_thickness, letter_s_dimensions[2]],
			// Fifth row.
			[0.0, letter_s_dimensions[1], letter_s_dimensions[2]],
			[letter_s_dimensions[0], letter_s_dimensions[1], letter_s_dimensions[2]]
		];
		let letter_s_triangle_indices = vec![
			// Upper layer
			[0, 4, 1],
			[0, 2, 4],
			[2, 6, 3],
			[2, 5, 6],
			[5, 10, 7],
			[5, 8, 10],	
			[9, 13, 10],
			[9, 12, 13],
			[11, 15, 13],
			[11, 14, 15],
			// Lower layer
			[16, 17, 20],
			[16, 20, 18],
			[18, 19, 22],
			[18, 22, 21],
			[21, 23, 26],
			[21, 26, 24],	
			[25, 26, 29],
			[25, 29, 28],
			[27, 29, 31],
			[27, 31, 30],
			// Vertical sides, seen from left.
			[16, 8, 0],
			[16, 24, 8],
			[25, 12, 9],
			[25, 28, 12],
			[27, 14, 11],
			[27, 30, 14],
			// Vertical sides, seen from right.
			[1, 20, 17],
			[1, 4, 20],
			[3, 22, 19],
			[3, 6, 22],
			[7, 31, 23],
			[7, 15, 31],
			// Horisontal sides, seen from the top.
			[0, 1, 17],
			[0, 17, 16],
			[6, 7, 23],
			[6, 23, 22],
			[11, 12, 28],
			[11, 28, 27],
			// Horisontal sides, seen from the bottom.
			[3, 19, 20],
			[3, 20, 4],
			[8, 24, 25],
			[8, 25, 9],
			[14, 30, 31],
			[14, 31, 15]
		];
		let letter_s_translation = [350.0, -200.0, 899.0];
		for node in letter_s_nodes.iter_mut() {
			*node = add(letter_s_translation, *node);
		}
		let mut physics_triangles_letter_s: Vec<PhysicsTriangle> = vec![];
		for letter_s_triangle_index in letter_s_triangle_indices {
			physics_triangles_letter_s.push(PhysicsTriangle::new(&letter_s_nodes, letter_s_triangle_index, white_semi_specular_transparent_emissive, physics_triangle, id, false));
		}
		id += 1;

		let letter_a_thickness = 75.0;
		let letter_a_dimensions = [300.0, 300.0, 100.0]; // W x H x D
		let mut letter_a_nodes = vec![
			// Upper layer.
			// Zeroth row
			[0.0, 0.0, 0.0],
			[letter_a_thickness, 0.0, 0.0],
			[letter_a_dimensions[0]-letter_a_thickness, 0.0, 0.0],
			[letter_a_dimensions[0], 0.0, 0.0],
			// First row.
			[letter_a_thickness, letter_a_thickness, 0.0],
			[letter_a_dimensions[0]-letter_a_thickness, letter_a_thickness, 0.0],
			// Second row.
			[letter_a_thickness, letter_a_dimensions[1]/2.0-letter_a_thickness/2.0, 0.0],
			[letter_a_dimensions[0]-letter_a_thickness, letter_a_dimensions[1]/2.0-letter_a_thickness/2.0, 0.0],
			// Thrid row.
			[letter_a_thickness, letter_a_dimensions[1]/2.0+letter_a_thickness/2.0, 0.0],
			[letter_a_dimensions[0]-letter_a_thickness, letter_a_dimensions[1]/2.0+letter_a_thickness/2.0, 0.0],
			// Fourth row.
			[0.0, letter_a_dimensions[1], 0.0],
			[letter_a_thickness, letter_a_dimensions[1], 0.0],
			[letter_a_dimensions[0]-letter_a_thickness, letter_a_dimensions[1], 0.0],
			[letter_a_dimensions[0], letter_a_dimensions[1], 0.0],
			// Lower layer.
			// Zeroth row
			[0.0, 0.0, letter_a_dimensions[2]],
			[letter_a_thickness, 0.0, letter_a_dimensions[2]],
			[letter_a_dimensions[0]-letter_a_thickness, 0.0, letter_a_dimensions[2]],
			[letter_a_dimensions[0], 0.0, letter_a_dimensions[2]],
			// First row.
			[letter_a_thickness, letter_a_thickness, letter_a_dimensions[2]],
			[letter_a_dimensions[0]-letter_a_thickness, letter_a_thickness, letter_a_dimensions[2]],
			// Second row.
			[letter_a_thickness, letter_a_dimensions[1]/2.0-letter_a_thickness/2.0, letter_a_dimensions[2]],
			[letter_a_dimensions[0]-letter_a_thickness, letter_a_dimensions[1]/2.0-letter_a_thickness/2.0, letter_a_dimensions[2]],
			// Thrid row.
			[letter_a_thickness, letter_a_dimensions[1]/2.0+letter_a_thickness/2.0, letter_a_dimensions[2]],
			[letter_a_dimensions[0]-letter_a_thickness, letter_a_dimensions[1]/2.0+letter_a_thickness/2.0, letter_a_dimensions[2]],
			// Fourth row.
			[0.0, letter_a_dimensions[1], letter_a_dimensions[2]],
			[letter_a_thickness, letter_a_dimensions[1], letter_a_dimensions[2]],
			[letter_a_dimensions[0]-letter_a_thickness, letter_a_dimensions[1], letter_a_dimensions[2]],
			[letter_a_dimensions[0], letter_a_dimensions[1], letter_a_dimensions[2]]
		];
		let letter_a_triangle_indices = vec![
			// Upper layer
			[0, 11, 1],
			[0, 10, 11],
			[2, 13, 3],
			[2, 12, 13],
			[1, 5, 2],
			[1, 4, 5],	
			[6, 9, 7],
			[6, 8, 9],
			// Lower layer
			[14, 15, 25],
			[14, 25, 24],
			[16, 17, 27],
			[16, 27, 26],
			[15, 16, 19],
			[15, 19, 18],	
			[20, 21, 23],
			[20, 23, 22],
			// Vertical sides, seen from left.
			[14, 10, 0],
			[14, 24, 10],
			[19, 7, 5],
			[19, 21, 7],
			[23, 12, 9],
			[23, 26, 12],
			// Vertical sides, seen from right.
			[4, 20, 18],
			[4, 6, 20],
			[8, 25, 22],
			[8, 11, 25],
			[3, 27, 17],
			[3, 13, 27],
			// Horisontal sides, seen from the top.
			[0, 3, 17],
			[0, 17, 14],
			[6, 7, 21],
			[6, 21, 20],
			// Horisontal sides, seen from the bottom.
			[4, 18, 19],
			[4, 19, 5],
			[8, 22, 23],
			[8, 23, 9],
			[10, 24, 25],
			[10, 25, 11],
			[12, 26, 27],
			[12, 27, 13]
		];
		let letter_a_translation = [725.0, -200.0, 899.0];
		for node in letter_a_nodes.iter_mut() {
			*node = add(letter_a_translation, *node);
		}
		let mut physics_triangles_letter_a: Vec<PhysicsTriangle> = vec![];
		for letter_a_triangle_index in letter_a_triangle_indices {
			physics_triangles_letter_a.push(PhysicsTriangle::new(&letter_a_nodes, letter_a_triangle_index, white_semi_specular_transparent_emissive, physics_triangle, id, false));
		}
		id += 1;

		let letter_r_thickness = 75.0;
		let letter_r_dimensions = [300.0, 300.0, 100.0]; // W x H x D
		let mut letter_r_nodes = vec![
			// Upper layer.
			// Zeroth row
			[0.0, 0.0, 0.0],
			[letter_r_thickness, 0.0, 0.0],
			[letter_r_dimensions[0]-letter_r_thickness, 0.0, 0.0],
			[letter_r_dimensions[0], 0.0, 0.0],
			// First row.
			[letter_r_thickness, letter_r_thickness, 0.0],
			[letter_r_dimensions[0]-letter_r_thickness, letter_r_thickness, 0.0],
			// Second row.
			[letter_r_thickness, letter_r_dimensions[1]/2.0-letter_r_thickness/2.0, 0.0],
			[letter_r_dimensions[0]-letter_r_thickness, letter_r_dimensions[1]/2.0-letter_r_thickness/2.0, 0.0],
			// Thrid row.
			[letter_r_thickness, letter_r_dimensions[1]/2.0+letter_r_thickness/2.0, 0.0],
			[letter_r_dimensions[0]/2.0, letter_r_dimensions[1]/2.0+letter_r_thickness/2.0, 0.0],
			[letter_r_dimensions[0]/2.0+letter_r_thickness, letter_r_dimensions[1]/2.0+letter_r_thickness/2.0, 0.0],
			[letter_r_dimensions[0]-letter_r_thickness, letter_r_dimensions[1]/2.0+letter_r_thickness/2.0, 0.0],
			[letter_r_dimensions[0], letter_r_dimensions[1]/2.0+letter_r_thickness/2.0, 0.0],
			// Fourth row.
			[0.0, letter_r_dimensions[1], 0.0],
			[letter_r_thickness, letter_r_dimensions[1], 0.0],
			[letter_r_dimensions[0]-letter_r_thickness, letter_r_dimensions[1], 0.0],
			[letter_r_dimensions[0], letter_r_dimensions[1], 0.0],
			// Lower layer.
			// Zeroth row
			[0.0, 0.0, letter_r_dimensions[2]],
			[letter_r_thickness, 0.0, letter_r_dimensions[2]],
			[letter_r_dimensions[0]-letter_r_thickness, 0.0, letter_r_dimensions[2]],
			[letter_r_dimensions[0], 0.0, letter_r_dimensions[2]],
			// First row.
			[letter_r_thickness, letter_r_thickness, letter_r_dimensions[2]],
			[letter_r_dimensions[0]-letter_r_thickness, letter_r_thickness, letter_r_dimensions[2]],
			// Second row.
			[letter_r_thickness, letter_r_dimensions[1]/2.0-letter_r_thickness/2.0, letter_r_dimensions[2]],
			[letter_r_dimensions[0]-letter_r_thickness, letter_r_dimensions[1]/2.0-letter_r_thickness/2.0, letter_r_dimensions[2]],
			// Thrid row.
			[letter_r_thickness, letter_r_dimensions[1]/2.0+letter_r_thickness/2.0, letter_r_dimensions[2]],
			[letter_r_dimensions[0]/2.0, letter_r_dimensions[1]/2.0+letter_r_thickness/2.0, letter_r_dimensions[2]],
			[letter_r_dimensions[0]/2.0+letter_r_thickness, letter_r_dimensions[1]/2.0+letter_r_thickness/2.0, letter_r_dimensions[2]],
			[letter_r_dimensions[0]-letter_r_thickness, letter_r_dimensions[1]/2.0+letter_r_thickness/2.0, letter_r_dimensions[2]],
			[letter_r_dimensions[0], letter_r_dimensions[1]/2.0+letter_r_thickness/2.0, letter_r_dimensions[2]],
			// Fourth row.
			[0.0, letter_r_dimensions[1], letter_r_dimensions[2]],
			[letter_r_thickness, letter_r_dimensions[1], letter_r_dimensions[2]],
			[letter_r_dimensions[0]-letter_r_thickness, letter_r_dimensions[1], letter_r_dimensions[2]],
			[letter_r_dimensions[0], letter_r_dimensions[1], letter_r_dimensions[2]]
		];
		let letter_r_triangle_indices = vec![
			// Upper layer
			[0, 13, 14],
			[0, 14, 1],
			[2, 11, 12],
			[2, 12, 3],
			[1, 4, 5],
			[1, 5, 2],	
			[6, 8, 11],
			[6, 11, 7],
			[9, 15, 16],
			[9, 16, 10],
			// Lower layer
			[17, 31, 30],
			[17, 18, 31],
			[19, 29, 28],
			[19, 20, 29],
			[18, 22, 21],
			[18, 19, 22],
			[23, 28, 25],
			[23, 24, 28],
			[26, 33, 32],
			[26, 27, 33],
			// Vertical and almost vertical sides, seen from left.
			[17, 13, 0],
			[17, 30, 13],
			[22, 7, 5],
			[22, 24, 7],
			[26, 15, 9],
			[26, 32, 15],
			// Vertical and almost vertical sides, seen from right.
			[4, 23, 21],
			[4, 6, 23],
			[8, 31, 25],
			[8, 14, 31],
			[3, 29, 20],
			[3, 12, 29],
			[10, 33, 27],
			[10, 16, 33],
			// Horisontal sides, seen from the top.
			[0, 3, 20],
			[0, 20, 17],
			[6, 7, 24],
			[6, 24, 23],
			// Horisontal sides, seen from the bottom.
			[4, 21, 22],
			[4, 22, 5],
			[13, 30, 31],
			[13, 31, 14],
			[8, 25, 26],
			[8, 26, 9],
			[15, 32, 33],
			[15, 33, 16],
			[10, 27, 29],
			[10, 29, 12]
		];
		let letter_r_translation = [1100.0, -200.0, 899.0];
		for node in letter_r_nodes.iter_mut() {
			*node = add(letter_r_translation, *node);
		}
		let mut physics_triangles_letter_r: Vec<PhysicsTriangle> = vec![];
		for letter_r_triangle_index in letter_r_triangle_indices {
			physics_triangles_letter_r.push(PhysicsTriangle::new(&letter_r_nodes, letter_r_triangle_index, white_semi_specular_transparent_emissive, physics_triangle, id, false));
		}

		let mut physics_triangles: Vec<PhysicsTriangle> = vec![];
		physics_triangles.append(&mut physics_triangles_walls);
		physics_triangles.append(&mut physics_triangles_letter_c);
		physics_triangles.append(&mut physics_triangles_letter_e);
		physics_triangles.append(&mut physics_triangles_letter_s);
		physics_triangles.append(&mut physics_triangles_letter_a);
		physics_triangles.append(&mut physics_triangles_letter_r);

		let light_spheres = vec![
			LightSphere::new(light_position, [1.0, 1.0, 1.0], light_radius),
		];
		let cameras = vec! [
		Camera::new([500.0, 500.0, -1000.0], [500.0, 500.0, -2000.0], [0.0, 0.0, 1.0], [1.0, 1.0, 1.0], 1400.0, 10.0),
		];
		Self {
			time: 0.0,
			time_per_frame,
			physics_cylinders,
			physics_spheres,
			physics_triangles,
			light_spheres,
			cameras,
		}
	}

	#[allow(dead_code)]
	pub fn simulate_physics_beat_saber(&mut self) {
		self.time += self.time_per_frame;
		for physics_cylinder in &mut self.physics_cylinders {
			let time_since_start = self.time - 30.0;
			let active = time_since_start > 0.0;
			physics_cylinder.active = active;
			if time_since_start > 0.0 {
				physics_cylinder.active = true;
				if physics_cylinder.id == 1 {
					let max_length = 1000.0;
					let mut length = time_since_start*max_length/5.0;
					if length > max_length {
						length = max_length;
					}
					physics_cylinder.length = length;
				}
			}
		}
		for physics_triangle in &mut self.physics_triangles {
			let start_time = physics_triangle.id as f64 * 4.0;
			let time_since_start = self.time - start_time;
			if time_since_start > 0.0 {
				physics_triangle.active = true;
				let max_emission = 0.06;
				physics_triangle.material.emission = mul(time_since_start.atan()*max_emission*2.0/std::f64::consts::PI, physics_triangle.material.color);
			}	
		}
	}

	#[allow(dead_code)]
	pub fn simulate_physics(&mut self) {
		self.time += self.time_per_frame;
		let number_of_timesteps: u64 = 1000;
		let dt = self.time_per_frame/(number_of_timesteps as f64);
		for _ in 0..number_of_timesteps {
			let mut time_left = dt;
			let mut finished = false;
			self.simulate_sphere_sphere_gravitation(dt);
			self.simulate_sphere_scene_gravitation(dt);
			// dt is already a small time period, but here we take potentially even smaller ones.
			// Only simulate until an intersection. Then update the velocity of the intersecting sphere
			// or the velocities of the pair of intersecting spheres. Continue doing that until all of
			// dt is spent. We don't apply gravity during this process, only before it starts.
			while !finished {
				let collision_data = self.compute_collision_data();
				match collision_data {
					// We found a collision. It either happens bebfore all of dt is spent, or after. Never take a step that
					// spends more than what remains of dt!
					Some(data) => {
						if data.2 < time_left {
							self.move_spheres(data.2);
							match data.4 {
								CollisionType::Sphere => {
									// See https://www.gamasutra.com/view/feature/131424/pool_hall_lessons_fast_accurate_.php?page=3
									// (Alternative, better version at http://realtimecollisiondetection.net/blog/?p=103 ? )
									let n = normalised(sub(self.physics_spheres[data.0].position, self.physics_spheres[data.1].position));
									let a1 = dot(self.physics_spheres[data.0].physics.velocity, n);
									let a2 = dot(self.physics_spheres[data.1].physics.velocity, n);
									let s = 2.0*(a1-a2)/(self.physics_spheres[data.0].physics.mass+self.physics_spheres[data.1].physics.mass);
									self.physics_spheres[data.0].physics.velocity = sub(self.physics_spheres[data.0].physics.velocity, mul(s*self.physics_spheres[data.1].physics.mass, n));
									self.physics_spheres[data.1].physics.velocity = add(self.physics_spheres[data.1].physics.velocity, mul(s*self.physics_spheres[data.0].physics.mass, n));
									time_left -= data.2;
								}
								CollisionType::Triangle => {
									let n = self.physics_triangles[data.1].normal;
									self.physics_spheres[data.0].physics.velocity = sub(self.physics_spheres[data.0].physics.velocity, mul(2.0*dot(self.physics_spheres[data.0].physics.velocity, n), n));
									time_left -= data.2;
								}
							}
						} else {
							self.move_spheres(time_left);
							finished = true;
						}
					}
					None => {
						self.move_spheres(time_left);
						finished = true;
					}
				}
			}
		}
	}

	#[allow(dead_code)]
	fn move_spheres(&mut self, time: f64) {
		for sphere in &mut self.physics_spheres {
			// Translation.
			sphere.position = add(sphere.position, mul(time, sphere.physics.velocity));
			// Rotation.
			let rotational_angle = sphere.physics.rotational_angular_speed*time;
			let sin_rot = (rotational_angle/2.0).sin();
			let q = [(rotational_angle/2.0).cos(), sin_rot*sphere.physics.rotational_axis[0], sin_rot*sphere.physics.rotational_axis[1], sin_rot*sphere.physics.rotational_axis[2]];
			let q_inverse = quarternion_inverse(q);
			let mut local_x = [0.0, sphere.physics.local_x[0], sphere.physics.local_x[1], sphere.physics.local_x[2]];
			let mut local_y = [0.0, sphere.physics.local_y[0], sphere.physics.local_y[1], sphere.physics.local_y[2]];
			let mut local_z = [0.0, sphere.physics.local_z[0], sphere.physics.local_z[1], sphere.physics.local_z[2]];
			local_x = quarternion_product(quarternion_product(q, local_x), q_inverse);
			local_y = quarternion_product(quarternion_product(q, local_y), q_inverse);
			local_z = quarternion_product(quarternion_product(q, local_z), q_inverse);
			local_x[0] = 0.0;
			local_y[0] = 0.0;
			local_z[0] = 0.0;
			sphere.physics.local_x = [local_x[1], local_x[2], local_x[3]];
			sphere.physics.local_y = [local_y[1], local_y[2], local_y[3]];
			sphere.physics.local_z = [local_z[1], local_z[2], local_z[3]];
		}
	}

	#[allow(dead_code)]
	fn compute_collision_data(&mut self) -> Option<(usize, usize, f64, [f64; 3], CollisionType)> {
		let sphere_sphere_collision_data = self.compute_sphere_sphere_collision();
		let sphere_triangle_collision_data = self.compute_sphere_triangle_collision();
		match sphere_sphere_collision_data {
			Some(ss_data) => {
				match sphere_triangle_collision_data {
					Some(st_data) => {
						if ss_data.2 < st_data.2 {
							return sphere_sphere_collision_data;
						}
						sphere_triangle_collision_data
					}
					None => {
						sphere_sphere_collision_data
					}
				}
			}
			None => {
				sphere_triangle_collision_data
			}
		}
	}

	#[allow(dead_code)]
	fn compute_sphere_triangle_collision(&self) -> Option<(usize, usize, f64, [f64; 3], CollisionType)> {
		// Index of the sphere, index of the triangle, time until collision, point of contact, CollisionType.
		let mut collision_data: Option<(usize, usize, f64, [f64; 3], CollisionType)> = None;
		for (i, sphere) in self.physics_spheres.iter().enumerate() {
			for (j, triangle) in self.physics_triangles.iter().enumerate() {
				if dot(sphere.physics.velocity, triangle.normal).abs() < 0.0 {
					println!("Should not be able to happen.");
					continue;
				}
				let distance_center_triangle = dot(sub(triangle.node0, sphere.position), triangle.normal);
				// Is the sphere moving away from the triangle? If so, no collision.
				let time_until_center_hits_plane = distance_center_triangle/dot(sphere.physics.velocity, triangle.normal);
				// Wrong direction.
				if time_until_center_hits_plane < 0.0 {
					continue;
				}
				// See http://mathinsight.org/distance_point_plane.
				let distance_center_triangle_unsigned = distance_center_triangle.abs();
				// Too close (why hasn't it been reflected already?).
				if distance_center_triangle_unsigned < sphere.radius {
					println!("Should have been reflected already.");
					continue;
				}
				let side = if dot(sphere.physics.velocity, triangle.normal) > 0.0 {
					-1.0
				} else {
					1.0
				};
				// The time until the actual sphere hits the triangle.
				let t = (side*sphere.radius + dot(sub(triangle.node0, sphere.position), triangle.normal))/dot(sphere.physics.velocity, triangle.normal);
				if t < 0.0 {
					println!("Wrong direction. Should have been detected earlier.");
					continue;
				}
				let sphere_center_during_intersection = add(mul(t, sphere.physics.velocity), sphere.position);
				if !(same_side(sphere_center_during_intersection, triangle.node0, triangle.node1, triangle.node2) && same_side(sphere_center_during_intersection, triangle.node1, triangle.node0, triangle.node2) && same_side(sphere_center_during_intersection, triangle.node2, triangle.node0, triangle.node1)) {
					continue;
				}
				// Is this correct?
				let intersection = sub(sphere_center_during_intersection, mul(side*sphere.radius, triangle.normal));
				match collision_data {
					Some(data) => {
						if t >= 0.0 && t < data.2 {
							collision_data = Some((i, j, t, intersection, CollisionType::Triangle));
						}
					}
					None => {
						if t >= 0.0 {
							collision_data = Some((i, j, t, intersection, CollisionType::Triangle));
						}
					}
				}
			}
		}
		collision_data
	}

	#[allow(dead_code)]
	fn compute_sphere_sphere_collision(&self) -> Option<(usize, usize, f64, [f64; 3], CollisionType)> {
		// Index of the first sphere, index of the second sphere, time until collision, point of contact, CollisionType.
		let mut collision_data: Option<(usize, usize, f64, [f64; 3], CollisionType)> = None;
		// Loop over all spheres to find the the future collision that occurs closest in time to the current time.
		for (i, sphere) in self.physics_spheres.iter().enumerate() {
			for (j, other_sphere) in self.physics_spheres.iter().enumerate() {
				if i == j {
					continue;
				}
				// To determine at which time the two spheres collide if the continue with constant velocity, we note that at the time t in the future we have
				// sphere.center(t) + sphere_center(0) + sphere.velocity*t
				// other_sphere.center(t) + other_sphere_center(0) + other_sphere.velocity*t
				// so collision occurs when (sphere.center(t) - other_sphere.center(t))^2 == (sphere.radius + other_sphere.radius)^2.
				// We solve this for t.
				let delta_rs = sub(sphere.position, other_sphere.position);
				let delta_v = sub(sphere.physics.velocity, other_sphere.physics.velocity);
				// If they're not approaching each other, no collision will occur.
				if dot(delta_v, delta_rs) > 0.0 {
					continue;
				}
				let p = 2.0*dot(delta_rs, delta_v)/dot(delta_v, delta_v);
				let q = (dot(delta_rs, delta_rs) - (sphere.radius+other_sphere.radius)*(sphere.radius+other_sphere.radius))/dot(delta_v, delta_v);
				let inside_square_root = p*p/4.0 - q;
				// If inside_square_root < 0, no collision will occur.
				if inside_square_root >= 0.0 {
					// Pick the smallest of t1 and t2 that is >= 0. This is the time until collision.
					let t1 = -1.0*p/2.0 - inside_square_root.sqrt();
					let t2 = -1.0*p/2.0 + inside_square_root.sqrt();
					let mut t: Option<f64> = None;
					match collision_data {
						Some(data) => {
							if t1 >= 0.0 && t1 < data.2 {
								t = Some(t1)
							} else if t2 >= 0.0 && t2 < data.2 {
								t = Some(t2)
							}
						}
						None => {
							if t1 >= 0.0 {
								t = Some(t1)
							} else if t2 >= 0.0 {
								t = Some(t2)
							}
						}
					}
					if let Some(time) = t {
						let sphere_center_after_t = add(sphere.position, mul(time, sphere.physics.velocity));
						let other_sphere_center_after_t = add(other_sphere.position, mul(time, other_sphere.physics.velocity));
						let intersection = add(mul(other_sphere.radius/(sphere.radius+other_sphere.radius), sphere_center_after_t), mul(sphere.radius/(sphere.radius+other_sphere.radius), other_sphere_center_after_t));
						collision_data = Some((i, j, time, intersection, CollisionType::Sphere));
					}
				}
			}
		}
		collision_data
	}

	#[allow(dead_code)]
	fn simulate_sphere_sphere_gravitation(&mut self, time: f64) {
		let gravitation_constant = 0.0;
		let mut accelerations: Vec<[f64; 3]> = Vec::new();
		for (i, sphere) in self.physics_spheres.iter().enumerate() {
			let mut total_force = [0.0, 0.0, 0.0];
			for (j, other_sphere) in self.physics_spheres.iter().enumerate() {
				if i == j {
					continue;
				}
				let distance_squared = distance_squared(sphere.position, other_sphere.position);
				let direction = normalised(sub(other_sphere.position, sphere.position));
				let force = mul(gravitation_constant*sphere.physics.mass*other_sphere.physics.mass/distance_squared, direction);
				total_force = add(total_force, force);
			}
			accelerations.push(mul(1.0/sphere.physics.mass, total_force));
		}
		for (i, sphere) in self.physics_spheres.iter_mut().enumerate() {
			sphere.physics.velocity = add(sphere.physics.velocity, mul(time, accelerations[i]));
		}
	}

	#[allow(dead_code)]
	fn simulate_sphere_scene_gravitation(&mut self, time: f64) {
		let gravitational_acceleration = [1.0, 0.0, 0.0];
		for sphere in &mut self.physics_spheres {
			sphere.physics.velocity = add(sphere.physics.velocity, mul(time, gravitational_acceleration));
		}
	}
}
