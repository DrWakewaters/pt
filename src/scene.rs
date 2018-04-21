use rand::{OsRng, Rng};
use pcg_rand::Pcg32;

use math::{add, distance, distance_squared, dot, mul, normalised, quarternion_inverse, quarternion_product, same_side, sub};
use primitive::{SphereFast, SphereLightsource, SpherePhysics, TriangleFast, TriangleLightsource, TrianglePhysics};

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct SceneForPhysics {
	pub triangles: Vec<TrianglePhysics>,
	pub spheres: Vec<SpherePhysics>,
}

#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
enum CollisionType {
	Sphere,
	Triangle,
}

impl SceneForPhysics {
	pub fn new() -> Self {
   		let nodes = vec![
        	// Top
			[0.0, 0.0, 0.0],           // 0
			[0.0, 0.0, 1000.0],        // 1
			[1000.0, 0.0, 1000.0],     // 2
			[1000.0, 0.0, 0.0],        // 3
			// Bottom
			[0.0, 1000.0, 0.0],        // 4
			[0.0, 1000.0, 1000.0],     // 5
			[1000.0, 1000.0, 1000.0],  // 6
			[1000.0, 1000.0, 0.0],     // 7
		];
		let size = 50000.0;
		let nodes_large_floor = vec![
			[500.0-size, 1000.0, 500.0-size],
			[-500.0-size, 1000.0, 500.0+size],
			[500.0+size, 1000.0, 500.0+size],
			[500.0+size, 1000.0, -500.0-size],
		];
		let grey = [0.9, 0.9, 0.9];
		let red = [0.9, 0.1, 0.1];
		let green = [0.1, 0.9, 0.1];
		let specular_probability_diffuse = 0.0;
		let specular_constant_diffuse = 1.0;
		let specular_probability_specular = 1.0;
		let specular_constant_specular = 0.3;
		let refractive_index = 20.0;
		let mut spheres: Vec<SpherePhysics> = Vec::new();
		let mut pcg: Pcg32 = OsRng::new().unwrap().gen();
		let spheres_to_add = 4;
		let mut remaining_spheres_to_add = spheres_to_add;
		let sphere_radius = 100.0;
		let color = [1.0, 1.0, 1.0];
		let color_alt = [0.0, 0.0, 0.0];
		let local_x = [1.0, 0.0, 0.0];
		let local_y = [0.0, 1.0, 0.0];
		let local_z = [0.0, 0.0, 1.0];
		let mass = 1000.0;
		let inertia_diagonal = 2.0/5.0*mass*sphere_radius*sphere_radius;
		let inertia = [[inertia_diagonal, 0.0, 0.0], [0.0, inertia_diagonal, 0.0], [0.0, 0.0, inertia_diagonal]];
		while remaining_spheres_to_add > 0 {
			let rx_u32: u32 = pcg.gen();
			let ry_u32: u32 = pcg.gen();
			let rz_u32: u32 = pcg.gen();
			let rvx_u32: u32 = pcg.gen();
			let rvy_u32: u32 = pcg.gen();
			let rvz_u32: u32 = pcg.gen();
			let mut rx = (1000.0-2.0*sphere_radius)*((rx_u32 as f64)/(<u32>::max_value() as f64) - 0.5) + 500.0;
			let mut ry = (1000.0-2.0*sphere_radius)*((ry_u32 as f64)/(<u32>::max_value() as f64) - 0.5) + 500.0;
			let mut rz = (1000.0-2.0*sphere_radius)*((rz_u32 as f64)/(<u32>::max_value() as f64) - 0.5) + 500.0;
			let mut rvx = 0.0*((rvx_u32 as f64)/(<u32>::max_value() as f64) - 0.5);
			let mut rvy = 0.0*((rvy_u32 as f64)/(<u32>::max_value() as f64) - 0.5);
			let mut rvz = 0.0*((rvz_u32 as f64)/(<u32>::max_value() as f64) - 0.5);
			if remaining_spheres_to_add == 4 {
				rx = 0.0;
				ry = 880.0;
				rz = 700.0;
			} else if remaining_spheres_to_add == 3 {
				rx = 400.0;
				ry = 880.0;
				rz = 1300.0;
			} else if remaining_spheres_to_add == 2 {
				rx = 800.0;
				ry = 880.0;
				rz = 1900.0;
			} else {
				rx = 1200.0;
				ry = 880.0;
				rz = 2500.0;
			}
			let rotational_axisx = 1.0*((rvx_u32 as f64)/(<u32>::max_value() as f64) - 0.5);
			let rotational_axisy = 1.0*((rvx_u32 as f64)/(<u32>::max_value() as f64) - 0.5);
			let rotational_axisz = 1.0*((rvx_u32 as f64)/(<u32>::max_value() as f64) - 0.5);
			let rotational_axis = normalised([rotational_axisx, rotational_axisy, rotational_axisz]);
			let rotational_angular_speed = 0.0*((rvz_u32 as f64)/(<u32>::max_value() as f64) - 0.5);
			let mut add_sphere = true;
			for sphere in spheres.iter() {
				if distance(sphere.center, [rx, ry, rz]) <= 2.0*sphere_radius {
					add_sphere = false;
				}
			}
			if add_sphere {
				spheres.push(SpherePhysics::new([rx, ry, rz], sphere_radius, color, color_alt, local_x, local_y, local_z, rotational_axis, rotational_angular_speed, specular_probability_diffuse, specular_constant_diffuse, refractive_index, false, [rvx, rvy, rvz], inertia, mass, false));
				remaining_spheres_to_add -= 1;
			}
		}
        let triangles = vec![
			// Left wall.
            //TrianglePhysics::new(&nodes, [0, 4, 5], red, specular_probability_diffuse, specular_constant_diffuse, refractive_index, false, false),
            //TrianglePhysics::new(&nodes, [0, 5, 1], red, specular_probability_diffuse, specular_constant_diffuse, refractive_index, false, false),
			// Far wall.
            //TrianglePhysics::new(&nodes, [1, 5, 6], grey, specular_probability_diffuse, specular_constant_diffuse, refractive_index, false, false),
            //TrianglePhysics::new(&nodes, [1, 6, 2], grey, specular_probability_diffuse, specular_constant_diffuse, refractive_index, false, false),
			// Right wall.            
			//TrianglePhysics::new(&nodes, [2, 6, 7], green, specular_probability_diffuse, specular_constant_diffuse, refractive_index, false, false),
            //TrianglePhysics::new(&nodes, [2, 7, 3], green, specular_probability_diffuse, specular_constant_diffuse, refractive_index, false, false),
            // Near wall.
			//TrianglePhysics::new(&nodes, [3, 7, 4], grey, specular_probability_diffuse, specular_constant_diffuse, refractive_index, false, true),
            //TrianglePhysics::new(&nodes, [3, 4, 0], grey, specular_probability_diffuse, specular_constant_diffuse, refractive_index, false, true),
			// Floor.
            //TrianglePhysics::new(&nodes, [4, 7, 6], grey, specular_probability_specular, specular_constant_specular, refractive_index, false, false),
            //TrianglePhysics::new(&nodes, [4, 6, 5], grey, specular_probability_specular, specular_constant_specular, refractive_index, false, false),
            //TrianglePhysics::new(&nodes_large_floor, [0, 3, 2], grey, specular_probability_specular, specular_constant_specular, refractive_index, false, false),
            //TrianglePhysics::new(&nodes_large_floor, [0, 2, 1], grey, specular_probability_specular, specular_constant_specular, refractive_index, false, false),
            // Ceiling.
            //TrianglePhysics::new(&nodes, [0, 1, 2], grey, specular_probability_diffuse, specular_constant_diffuse, refractive_index, false, false),
            //TrianglePhysics::new(&nodes, [0, 2, 3], grey, specular_probability_diffuse, specular_constant_diffuse, refractive_index, false, false),
            TrianglePhysics::new(&nodes, [0, 1, 2], grey, specular_probability_diffuse, specular_constant_diffuse, refractive_index, true, false),
            TrianglePhysics::new(&nodes, [0, 2, 3], grey, specular_probability_diffuse, specular_constant_diffuse, refractive_index, true, false),
        ];
		Self {
			triangles,
			spheres,
		}
	}
	
	pub fn simulate_physics(&mut self) {
		let time: f64 = 0.5;
		let number_of_timesteps: u64 = 1000;
		let dt = time/(number_of_timesteps as f64);
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
									let n = normalised(sub(self.spheres[data.0].center, self.spheres[data.1].center));
									let a1 = dot(self.spheres[data.0].velocity, n);
									let a2 = dot(self.spheres[data.1].velocity, n);
									let s = 2.0*(a1-a2)/(self.spheres[data.0].mass+self.spheres[data.1].mass);
									self.spheres[data.0].velocity = sub(self.spheres[data.0].velocity, mul(s*self.spheres[data.1].mass, n));
									self.spheres[data.1].velocity = add(self.spheres[data.1].velocity, mul(s*self.spheres[data.0].mass, n));
									time_left -= data.2;
								}
								CollisionType::Triangle => {
									let n = self.triangles[data.1].normal;
									self.spheres[data.0].velocity = sub(self.spheres[data.0].velocity, mul(2.0*dot(self.spheres[data.0].velocity, n), n));
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
	
	fn move_spheres(&mut self, time: f64) {
		for sphere in self.spheres.iter_mut() {
			// Translation.
			sphere.center = add(sphere.center, mul(time, sphere.velocity));
			// Rotation.
			let rotational_angle = sphere.rotational_angular_speed*time;
			let sin_rot = (rotational_angle/2.0).sin();
			let q = [(rotational_angle/2.0).cos(), sin_rot*sphere.rotational_axis[0], sin_rot*sphere.rotational_axis[1], sin_rot*sphere.rotational_axis[2]];
			let q_inverse = quarternion_inverse(q);
			let mut local_x = [0.0, sphere.local_x[0], sphere.local_x[1], sphere.local_x[2]];
			let mut local_y = [0.0, sphere.local_y[0], sphere.local_y[1], sphere.local_y[2]];
			let mut local_z = [0.0, sphere.local_z[0], sphere.local_z[1], sphere.local_z[2]];
			local_x = quarternion_product(quarternion_product(q, local_x), q_inverse);
			local_y = quarternion_product(quarternion_product(q, local_y), q_inverse);
			local_z = quarternion_product(quarternion_product(q, local_z), q_inverse);
			local_x[0] = 0.0;
			local_y[0] = 0.0;
			local_z[0] = 0.0;
			sphere.local_x = [local_x[1], local_x[2], local_x[3]];
			sphere.local_y = [local_y[1], local_y[2], local_y[3]];
			sphere.local_z = [local_z[1], local_z[2], local_z[3]];
		}
	}
	
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
						return sphere_triangle_collision_data;
					}
					None => {
						return sphere_sphere_collision_data;
					}
				}
			}
			None => {
				return sphere_triangle_collision_data;
			}
		}
	}
	
	fn compute_sphere_triangle_collision(&self) -> Option<(usize, usize, f64, [f64; 3], CollisionType)> {
		// Index of the sphere, index of the triangle, time until collision, point of contact, CollisionType.
		let mut collision_data: Option<(usize, usize, f64, [f64; 3], CollisionType)> = None;
		for (i, sphere) in self.spheres.iter().enumerate() {
			for (j, triangle) in self.triangles.iter().enumerate() {
				if dot(sphere.velocity, triangle.normal).abs() < 0.0 {
					println!("Should not be able to happen.");
					continue;
				}
				let distance_center_triangle = dot(sub(triangle.node0, sphere.center), triangle.normal);
				// Is the sphere moving away from the triangle? If so, no collision.
				let time_until_center_hits_plane = distance_center_triangle/dot(sphere.velocity, triangle.normal);
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
				let mut side = 1.0;
				if dot(sphere.velocity, triangle.normal) > 0.0 {
					side = -1.0;
				}
				// The time until the actual sphere hits the triangle.
				let t = (side*sphere.radius + dot(sub(triangle.node0, sphere.center), triangle.normal))/dot(sphere.velocity, triangle.normal);
				if t < 0.0 {
					println!("Wrong direction. Should have been detected earlier.");
					continue;
				}
				let sphere_center_during_intersection = add(mul(t, sphere.velocity), sphere.center);
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
	
	fn compute_sphere_sphere_collision(&self) -> Option<(usize, usize, f64, [f64; 3], CollisionType)> {
		// Index of the first sphere, index of the second sphere, time until collision, point of contact, CollisionType.
		let mut collision_data: Option<(usize, usize, f64, [f64; 3], CollisionType)> = None;
		// Loop over all spheres to find the the future collision that occurs closest in time to the current time.
		for (i, sphere) in self.spheres.iter().enumerate() {
			for (j, other_sphere) in self.spheres.iter().enumerate() {
				if i == j {
					continue;
				}
				// To determine at which time the two spheres collide if the continue with constant velocity, we note that at the time t in the future we have
				// sphere.center(t) + sphere_center(0) + sphere.velocity*t
				// other_sphere.center(t) + other_sphere_center(0) + other_sphere.velocity*t
				// so collision occurs when (sphere.center(t) - other_sphere.center(t))^2 == (sphere.radius + other_sphere.radius)^2.
				// We solve this for t.
				let delta_rs = sub(sphere.center, other_sphere.center);
				let delta_v = sub(sphere.velocity, other_sphere.velocity);
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
					match t {
						Some(time) => {
							let sphere_center_after_t = add(sphere.center, mul(time, sphere.velocity));
							let other_sphere_center_after_t = add(other_sphere.center, mul(time, other_sphere.velocity));
							let intersection = add(mul(other_sphere.radius/(sphere.radius+other_sphere.radius), sphere_center_after_t), mul(sphere.radius/(sphere.radius+other_sphere.radius), other_sphere_center_after_t));
							collision_data = Some((i, j, time, intersection, CollisionType::Sphere));	
						}
						None => {
						}
					}
				}
			}
		}
		collision_data		
	}
	
	fn simulate_sphere_sphere_gravitation(&mut self, time: f64) {
		let gravitation_constant = 0.0;
		let mut accelerations: Vec<[f64; 3]> = Vec::new();
		for (i, sphere) in self.spheres.iter().enumerate() {
			let mut total_force = [0.0, 0.0, 0.0];
			for (j, other_sphere) in self.spheres.iter().enumerate() {
				if i == j {
					continue;
				}
				let distance_squared = distance_squared(sphere.center, other_sphere.center);
				let direction = normalised(sub(other_sphere.center, sphere.center));
				let force = mul(gravitation_constant*sphere.mass*other_sphere.mass/distance_squared, direction);
				total_force = add(total_force, force);
			}
			accelerations.push(mul(1.0/sphere.mass, total_force));
		}
		for (i, sphere) in self.spheres.iter_mut().enumerate() {
			sphere.velocity = add(sphere.velocity, mul(time, accelerations[i]));
		}
	}
	
	// For fun. Gravitation towards the closet surface.
	/*
	fn simulate_sphere_scene_gravitation(&mut self, time: f64) {
		let gravitational_acceleration_strength = -3.0;
		for sphere in self.spheres.iter_mut() {
			let mut closest_triangle_index = 0;
			let mut closest_triangle_distance = MAX;
			for (j, triangle) in self.triangles.iter().enumerate() {
				let distance_center_triangle = dot(sub(triangle.node0, sphere.center), triangle.normal);
				if distance_center_triangle < closest_triangle_distance {
					closest_triangle_distance = distance_center_triangle;
					closest_triangle_index = j;
				}
			}
			let gravitational_acceleration = mul(gravitational_acceleration_strength, self.triangles[closest_triangle_index].normal);
			sphere.velocity = add(sphere.velocity, mul(time, gravitational_acceleration));
		}
	}
	*/
	
	
	fn simulate_sphere_scene_gravitation(&mut self, time: f64) {
		let gravitational_acceleration = [1.0, 0.0, 0.0];
		for sphere in self.spheres.iter_mut() {
			sphere.velocity = add(sphere.velocity, mul(time, gravitational_acceleration));
		}
	}
}

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
			let surface_triangle_for_rendering = TriangleFast::new(&triangle);
			triangle_surfaces.push(surface_triangle_for_rendering);
		}
		for sphere in scene_for_physics.spheres {
			if sphere.is_lightsource {
				let lightsource_sphere_for_rendering = SphereLightsource::new(&sphere);
				sphere_lightsources.push(lightsource_sphere_for_rendering);				
			}
			let surface_sphere_for_rendering = SphereFast::new(&sphere);
			sphere_surfaces.push(surface_sphere_for_rendering);
		}
		Self {
			triangle_lightsources,
			triangle_surfaces,
			sphere_lightsources,
			sphere_surfaces,
		}
	}
}
