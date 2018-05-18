use collisiontype::CollisionType;
use math::{add, distance_squared, dot, mul, normalised, quarternion_inverse, quarternion_product, same_side, sub};
use light::Light;
use material::Material;
use physics::Physics;
use physicssphere::PhysicsSphere;
use physicstriangle::PhysicsTriangle;

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct PhysicsScene {
	pub physics_triangles: Vec<PhysicsTriangle>,
	pub physics_spheres: Vec<PhysicsSphere>,
	pub lights: Vec<Light>,
}

impl PhysicsScene {
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
		let grey_diffuse_opaque = Material::new([1.0, 1.0, 1.0], 0.0, 0.0, 0.0, true);
		let grey_specular_opaque = Material::new([1.0, 1.0, 1.0], 1.0, 0.1, 0.0, true);
		//let grey_diffuse_transparent = Material::new([0.9, 0.9, 0.9], 0.0, 0.0, 2.0, false);
		//let grey_specular_transparent = Material::new([0.9, 0.9, 0.9], 1.0, 0.1, 2.0, false);

		let local_x = [1.0, 0.0, 0.0];
		let local_y = [0.0, 1.0, 0.0];
		let local_z = [0.0, 0.0, 1.0];

		let large_radius = 150.0;
		let middle_radius = 100.0;
		let mass = 1000.0;

		let inertia_diagonal_large_radius = 2.0/5.0*mass*large_radius*large_radius;
		let inertia_large_radius = [[inertia_diagonal_large_radius, 0.0, 0.0], [0.0, inertia_diagonal_large_radius, 0.0], [0.0, 0.0, inertia_diagonal_large_radius]];
		//let inertia_diagonal_middle_radius = 2.0/5.0*mass*middle_radius*middle_radius;
		//let inertia_middle_radius = [[inertia_diagonal_middle_radius, 0.0, 0.0], [0.0, inertia_diagonal_middle_radius, 0.0], [0.0, 0.0, inertia_diagonal_middle_radius]];

		let physics_large_radius = Physics::new(local_x, local_y, local_z, [1.0, 0.0, 0.0], 0.0, [0.0, 0.0, 0.0], inertia_large_radius, mass, true);
		//let physics_middle_radius = Physics::new(local_x, local_y, local_z, [1.0, 0.0, 0.0], 0.0, [0.0, 0.0, 0.0], physics_middle_radius, mass, true);
		let physics_triangle = Physics::new(local_x, local_y, local_z, [1.0, 0.0, 0.0], 0.0, [0.0, 0.0, 0.0], [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]], 0.0, false);

		let sphere_red_specular_opaque = PhysicsSphere::new([700.0, 1000.0-(large_radius+1.0e-3), 700.0], large_radius, grey_specular_opaque, physics_large_radius);
		let sphere_red_diffuse_opaque = PhysicsSphere::new([300.0, 1000.0-(large_radius+1.0e-3), 700.0], large_radius, grey_diffuse_opaque, physics_large_radius);
		//let sphere_grey_specular_transparent_1 = PhysicsSphere::new([600.0, 1000.0-(middle_radius+1.0e-3), 450.0], middle_radius, grey_specular_transparent, physics_1);
		//let sphere_grey_specular_transparent_2 = PhysicsSphere::new([850.0, 1000.0-(middle_radius+1.0e-3), 300.0], middle_radius, grey_specular_transparent, physics_1);

		let physics_spheres = vec![sphere_red_specular_opaque, sphere_red_diffuse_opaque];

		let physics_triangles = vec![
			// Left wall.
            PhysicsTriangle::new(&nodes, [0, 4, 5], grey_diffuse_opaque, physics_triangle),
            PhysicsTriangle::new(&nodes, [0, 5, 1], grey_diffuse_opaque, physics_triangle),
			// Far wall.
            PhysicsTriangle::new(&nodes, [1, 5, 6], grey_diffuse_opaque, physics_triangle),
            PhysicsTriangle::new(&nodes, [1, 6, 2], grey_diffuse_opaque, physics_triangle),
			// Right wall.
			PhysicsTriangle::new(&nodes, [2, 6, 7], grey_diffuse_opaque, physics_triangle),
            PhysicsTriangle::new(&nodes, [2, 7, 3], grey_diffuse_opaque, physics_triangle),
            // Near wall.
			//PhysicsTriangle::new(&nodes, [3, 7, 4], grey_diffuse_opaque, physics_2),
            //PhysicsTriangle::new(&nodes, [3, 4, 0], grey_diffuse_opaque, physics_2),
			// Floor.
            PhysicsTriangle::new(&nodes, [4, 7, 6], grey_diffuse_opaque, physics_triangle),
            PhysicsTriangle::new(&nodes, [4, 6, 5], grey_diffuse_opaque, physics_triangle),
            // Ceiling.
            PhysicsTriangle::new(&nodes, [0, 1, 2], grey_diffuse_opaque, physics_triangle),
            PhysicsTriangle::new(&nodes, [0, 2, 3], grey_diffuse_opaque, physics_triangle),
        ];
		let lights = vec![
			Light::new([500.0, 1000.0-(middle_radius+1.0e-3), 400.0], [1.0, 1.0, 1.0]),
		];
		Self {
			physics_spheres,
			physics_triangles,
			lights,
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

	fn simulate_sphere_scene_gravitation(&mut self, time: f64) {
		let gravitational_acceleration = [1.0, 0.0, 0.0];
		for sphere in &mut self.physics_spheres {
			sphere.physics.velocity = add(sphere.physics.velocity, mul(time, gravitational_acceleration));
		}
	}
}
