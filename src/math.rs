use std::f64::consts::PI;

use pcg_rand::Pcg32;
use rand::Rng;

use crate::GAMMA;

use crate::hitpoint::Hitpoint;
use crate::lightsphere::LightSphere;
use crate::material::Material;

// BASIC MATH

#[inline(always)]
#[allow(dead_code)]
pub fn add(left: [f64; 3], right: [f64; 3]) -> [f64; 3] {
	[left[0]+right[0], left[1]+right[1], left[2]+right[2]]
}

#[inline(always)]
#[allow(dead_code)]
pub fn add_f32(left: [f32; 3], right: [f32; 3]) -> [f32; 3] {
	[left[0]+right[0], left[1]+right[1], left[2]+right[2]]
}

#[inline(always)]
#[allow(dead_code)]
pub fn add_f64(left: [f64; 3], right: [f64; 3]) -> [f64; 3] {
	[left[0]+right[0], left[1]+right[1], left[2]+right[2]]
}

#[inline(always)]
#[allow(dead_code)]
pub fn add_u16(left: [u16; 3], right: [u16; 3]) -> [u16; 3] {
	[left[0]+right[0], left[1]+right[1], left[2]+right[2]]
}

#[inline(always)]
#[allow(dead_code)]
pub fn sub(left: [f64; 3], right: [f64; 3]) -> [f64; 3] {
	[left[0]-right[0], left[1]-right[1], left[2]-right[2]]
}

#[inline(always)]
#[allow(dead_code)]
pub fn mul(scalar: f64, vector: [f64; 3]) -> [f64; 3] {
	[scalar*vector[0], scalar*vector[1], scalar*vector[2]]
}

#[inline(always)]
#[allow(dead_code)]
pub fn elementwise_mul(left: [f64; 3], right: [f64; 3]) -> [f64; 3] {
[left[0]*right[0], left[1]*right[1], left[2]*right[2]]
}

#[inline(always)]
#[allow(dead_code)]
pub fn dot(left: [f64; 3], right: [f64; 3]) -> f64 {
	left[0]*right[0]+left[1]*right[1]+left[2]*right[2]
}

#[inline(always)]
#[allow(dead_code)]
pub fn cross(left: [f64; 3], right: [f64; 3]) -> [f64; 3] {
	[left[1]*right[2]-left[2]*right[1], left[2]*right[0]-left[0]*right[2], left[0]*right[1]-left[1]*right[0]]
}

#[inline(always)]
#[allow(dead_code)]
pub fn norm(vector: [f64; 3]) -> f64 {
	(vector[0]*vector[0]+vector[1]*vector[1]+vector[2]*vector[2]).sqrt()
}

#[inline(always)]
#[allow(dead_code)]
pub fn distance(left: [f64; 3], right: [f64; 3]) -> f64 {
	let vector = sub(left, right);
	(vector[0]*vector[0]+vector[1]*vector[1]+vector[2]*vector[2]).sqrt()
}

#[inline(always)]
#[allow(dead_code)]
pub fn distance_squared(left: [f64; 3], right: [f64; 3]) -> f64 {
	let vector = sub(left, right);
	vector[0]*vector[0]+vector[1]*vector[1]+vector[2]*vector[2]
}

#[inline(always)]
#[allow(dead_code)]
pub fn normalised(vector: [f64; 3]) -> [f64; 3] {
	let norm = (vector[0]*vector[0]+vector[1]*vector[1]+vector[2]*vector[2]).sqrt();
	[vector[0]/norm, vector[1]/norm, vector[2]/norm]
}

#[inline(always)]
#[allow(dead_code)]
pub fn invert(matrix: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
	let determinant = matrix[0][0]*(matrix[1][1]*matrix[2][2]-matrix[1][2]*matrix[2][1]) - matrix[1][0]*(matrix[0][1]*matrix[2][2]-matrix[0][2]*matrix[2][1]) + matrix[2][0]*(matrix[0][1]*matrix[1][2]-matrix[0][2]*matrix[1][1]);
	let determinant_inverted = 1.0/determinant;
	[[determinant_inverted*matrix[1][1]*matrix[2][2]-determinant_inverted*matrix[2][1]*matrix[1][2],
	determinant_inverted*matrix[0][2]*matrix[2][1]-determinant_inverted*matrix[2][2]*matrix[0][1],
	determinant_inverted*matrix[0][1]*matrix[1][2]-determinant_inverted*matrix[1][1]*matrix[0][2]],
	[determinant_inverted*matrix[1][2]*matrix[2][0]-determinant_inverted*matrix[2][2]*matrix[1][0],
	determinant_inverted*matrix[0][0]*matrix[2][2]-determinant_inverted*matrix[2][0]*matrix[0][2],
	determinant_inverted*matrix[0][2]*matrix[1][0]-determinant_inverted*matrix[1][2]*matrix[0][0]],
	[determinant_inverted*matrix[1][0]*matrix[2][1]-determinant_inverted*matrix[2][0]*matrix[1][1],
	determinant_inverted*matrix[0][1]*matrix[2][0]-determinant_inverted*matrix[2][1]*matrix[0][0],
	determinant_inverted*matrix[0][0]*matrix[1][1]-determinant_inverted*matrix[1][0]*matrix[0][1]]]
}

#[inline(always)]
#[allow(dead_code)]
pub fn matrix_vector_mul(matrix: [[f64; 3]; 3] , vector: [f64; 3]) -> [f64; 3] {
	[matrix[0][0]*vector[0] + matrix[1][0]*vector[1] + matrix[2][0]*vector[2],
	matrix[0][1]*vector[0] + matrix[1][1]*vector[1] + matrix[2][1]*vector[2],
	matrix[0][2]*vector[0] + matrix[1][2]*vector[1] + matrix[2][2]*vector[2]]
}

// The solution to matrix*vector' = vector is vector' = matrix^(-1)*vector.
#[inline(always)]
#[allow(dead_code)]
pub fn solve_linear(matrix: [[f64; 3]; 3] , vector: [f64; 3]) -> [f64; 3] {
	let matrix_inverse = invert(matrix);
	matrix_vector_mul(matrix_inverse, vector)
}

#[inline(always)]
#[allow(dead_code)]
pub fn quarternion_inverse(quarternion: [f64; 4]) -> [f64; 4] {
	[quarternion[0], -1.0*quarternion[1], -1.0*quarternion[2], -1.0*quarternion[3]]
}

#[inline(always)]
#[allow(dead_code)]
pub fn quarternion_product(left: [f64; 4], right: [f64; 4]) -> [f64; 4] {
	[left[0]*right[0] - left[1]*right[1] - left[2]*right[2] - left[3]*right[3],
	 left[0]*right[1] + left[1]*right[0] + left[2]*right[3] - left[3]*right[2],
	 left[0]*right[2] - left[1]*right[3] + left[2]*right[0] + left[3]*right[1],
	 left[0]*right[3] + left[1]*right[2] - left[2]*right[1] + left[3]*right[0]]
}

// Given a normalised vector, compute an ON basis where this vector is a basis vector.
#[allow(dead_code)]
pub fn make_basis_vectors(direction: [f64; 3]) -> ([f64; 3], [f64; 3]) {
    let x = [1.0, 0.0, 0.0];
    let y = [0.0, 1.0, 0.0];
	let dot_product = dot(direction, x);
	let other_direction = if dot_product.abs() > 0.1 {
		y
	} else {
		x
	};
	let t_1 = normalised(cross(direction, other_direction));
	let t_2 = normalised(cross(direction, t_1));
	(t_1, t_2)
}

// Given a normalised vector, compute an ON basis where this vector is a basis vector.
#[allow(dead_code)]
pub fn make_basis_vectors_from_two(direction: [f64; 3], other_direction: [f64; 3]) -> ([f64; 3], [f64; 3]) {
	let t_1 = normalised(sub(other_direction, mul(dot(direction, other_direction), direction)));
	(t_1, normalised(cross(direction, t_1)))
}

// PICK RANDOM POINTS

// See http://mathworld.wolfram.com/SpherePointPicking.html.
#[allow(dead_code)]
pub fn random_uniform_on_sphere(pcg: &mut Pcg32) -> [f64; 3] {
	let r1 = pcg.gen::<f64>();
	let r2 = pcg.gen::<f64>();

    let theta = 2.0*PI*r1;
	let u = 2.0*(r2-0.5);
	let p = (1.0-u*u).sqrt();
	
	normalised([p*theta.cos(), p*theta.sin(), u])
}

#[allow(dead_code)]
pub fn random_lambertian_on_hemisphere(direction: [f64; 3], t_1: [f64; 3], t_2: [f64; 3], pcg: &mut Pcg32) -> [f64; 3] {
	let r1 = pcg.gen::<f64>();
	let r2 = pcg.gen::<f64>();

	let phi = 2.0*PI*r2;
	let sin_theta = (1.0-r1).sqrt();
	let cos_theta = r1.sqrt();
	
	normalised(add(add(mul(sin_theta*phi.cos(), t_1), mul(sin_theta*phi.sin(), t_2)), mul(cos_theta, direction)))
}

#[allow(dead_code)]
pub fn random_specular_on_hemisphere(direction: [f64; 3], t_1: [f64; 3], t_2: [f64; 3], maximum_angle: f64, pcg: &mut Pcg32) -> [f64; 3] {
	let r_1 = pcg.gen::<f64>();
	let r_2 = pcg.gen::<f64>();

    let phi = 2.0*PI*r_1;
    let z_min = maximum_angle.cos();
    let u = z_min + (1.0-z_min)*r_2;
    let theta = u.acos();
	let sin_theta = theta.sin();
	let cos_theta = theta.cos();

	normalised(add(add(mul(sin_theta*phi.cos(), t_1), mul(sin_theta*phi.sin(), t_2)), mul(cos_theta, direction)))
}

#[allow(dead_code)]
pub fn rotate(vector: [f64; 3], rotation_axis: [f64; 3], angle: f64) -> [f64; 3] {
	add(add(mul(angle.cos(), vector), mul(angle.sin(), cross(rotation_axis, vector))), mul(dot(rotation_axis, vector)*(1.0-angle.cos()), rotation_axis))
}

// PICK RANDOM DIRECTIONS FROM BRDFS

// A Ray hits a RendererShape.
// 1: Pick a micronormal. That normal determines everything that happens afterwards. The micronormal can't differ from the macronormal
// by more than 90 degrees. The micronormal is either identical to the macro normal, or uniformly (?) picked on a hemisphere around the macronormal.
// 2: Given that micronormal, determine if there is a transmission or reflection.
// 3a: If there is a reflection, reflect according to the reflection law.
// 3b: If there is a transmission, transmit according to Snell's law.
// Since the micronormal differs from the macronormal, strange things may happen.
// Reflection may happen so that the outgoing ray hits the surface again. If so: repeat the whole process until it escapes the surface (by reflection or transmission). This case is easy.
// Transmission might happen so that the ray is moving back the same way it came. If so, the process must be repeated, but with refractive_index_1 and refractive_index_2 inverted. This case is thus somewhat trickier.
#[allow(dead_code)]
pub fn random_from_brdf(mut hitpoint: &mut Hitpoint, light_spheres: &mut Vec<LightSphere>, light_sampling_probability: f64, light_sampling_possible: &mut bool,  brdf_modifier: &mut f64, pcg: &mut Pcg32) {

	// @TODO: Support for multiple lights.
	let direction_to_light_center = sub(light_spheres[0].position, hitpoint.position);
	let distance_to_light_center = norm(direction_to_light_center);
	let normalised_direction_to_light_center = mul(1.0/distance_to_light_center, direction_to_light_center);
	let light_angle = (light_spheres[0].radius/distance_to_light_center).atan();

	let inverted_incoming_direction = mul(-1.0, hitpoint.incoming_direction);

	// Let the specular reflection decrease in a lambertian fashion - the probability goes to zero as the micronormal reaches it's maximum "perturbation" (modification of the macronormal).
	let (t_1_specular, t_2_specular) = make_basis_vectors(hitpoint.normal);
	let angle_to_surface = PI/2.0-dot(inverted_incoming_direction, hitpoint.normal).acos();
	let max_angle = if hitpoint.material.maximum_specular_angle < angle_to_surface {
		hitpoint.material.maximum_specular_angle/2.0
	} else {
		angle_to_surface/2.0
	};
	let proportion = max_angle/(PI/2.0);
	let random_lambertian_specular = random_lambertian_on_hemisphere(hitpoint.normal, t_1_specular, t_2_specular, pcg);	

	let half_vector = normalised(add(hitpoint.normal, inverted_incoming_direction));
	let (t_1_lambertian, t_2_lambertian) = make_basis_vectors(half_vector);
	let random_lambertian_lambertian = random_lambertian_on_hemisphere(half_vector, t_1_lambertian, t_2_lambertian, pcg);

	let mut r = pcg.gen::<f64>();
	let micro_normal = if r < hitpoint.material.specular_probability {
		normalised(add(mul(1.0-proportion, hitpoint.normal), mul(proportion, random_lambertian_specular)))
	} else {
		normalised(add(half_vector, random_lambertian_lambertian))
	};

	// Compute the probability of a reflection.
	r = pcg.gen::<f64>();
	let (r_1, r_2) = if hitpoint.hit_from_outside {
		(1.0, hitpoint.material.refractive_index)
	} else {
		(hitpoint.material.refractive_index, 1.0)
	};
	let reflect = hitpoint.material.is_opaque || r < compute_f(inverted_incoming_direction, micro_normal, r_1, r_2);
	if reflect {
		let new_incoming_direction_if_reflected = normalised(sub(mul(2.0*dot(inverted_incoming_direction, micro_normal), micro_normal), inverted_incoming_direction));
		let normal_light_angle = dot(hitpoint.normal, normalised_direction_to_light_center).acos();
		if normal_light_angle > 0.0 && normal_light_angle < PI/2.0 && PI/2.0 - normal_light_angle > light_angle {
			*light_sampling_possible = true;
			let proportion = light_angle/(PI/2.0);
			let (t_1, t_2) = make_basis_vectors(normalised_direction_to_light_center);
			let random_lambertian = random_lambertian_on_hemisphere(normalised_direction_to_light_center, t_1, t_2, pcg);
			let random_direction_to_light = normalised(add(mul(1.0-proportion, normalised_direction_to_light_center), mul(proportion, random_lambertian)));
			// @TODO: Invert one of them?
			// Lambertian: to give a reclection along random_direction_to_light the half vector need to be perturbed a certain amount.
			let lambertian_micro_normal_to_light = normalised(add(random_direction_to_light, inverted_incoming_direction));
			let lambertian_brdf = brdf_lambertian_in_cone(micro_normal, lambertian_micro_normal_to_light, PI/4.0);
			//let specular_brdf = brdf_lambertian(, );
			//let correct = (1.0-hitpoint.material.specular_probability)*lambertian_brdf + hitpoint.material.specular_probability*specular_brdf;
			let correct = lambertian_brdf;
			let light_brdf = 1.0;
			let actual = (1.0-light_sampling_probability)*correct + light_sampling_probability*light_brdf;
			*brdf_modifier = if correct.abs() < 1.0e-6 && actual.abs() < 1.0e-6 {
				1.0
			} else {
				correct/actual
			};
			//println!("{}, {}, {}", *brdf_modifier, correct, actual);
			r = pcg.gen::<f64>();
			if r < light_sampling_probability {
				hitpoint.incoming_direction = random_direction_to_light;
			} else {
				hitpoint.incoming_direction = new_incoming_direction_if_reflected;
			}
		} else {
			hitpoint.incoming_direction = new_incoming_direction_if_reflected;
		}
	} else {
		// See https://en.wikipedia.org/wiki/Snell%27s_law.
		let r = r_1/r_2;
		let l = hitpoint.incoming_direction;
		let c = -1.0*dot(micro_normal, l);
		let radicand = 1.0-r*r*(1.0-c*c);
		let transmission_direction = normalised(add(mul(r, l), mul(r*c-radicand.sqrt(), micro_normal)));
		hitpoint.incoming_direction = transmission_direction;
		hitpoint.transmitted = true;
		hitpoint.hit_from_outside = !hitpoint.hit_from_outside;
		if dot(transmission_direction, hitpoint.normal) > 0.0 {
			hitpoint.transmitted = false;
			/*
			hitpoint.normal = mul(-1.0, hitpoint.normal);
			random_from_brdf(hitpoint, pcg);
			*/
		}
	}
}

// rP: probability that an interation is a reflection
// sP: probability that a reflection is specular
// lBr: lambertian_brdf_reflection
// sBr: specular_brdf_reflection
// lBt: lambertian_brdf_transmission
// sBt: specular_brdf_transmission
// aP: light_sampling_probability
// aB: light_brdf

// Return correct/actual (correct probability density / the actual probability density).
// With transmission:
// correct = (1-sP)*lBr*rP + sP*sBr*rP + (1-sP)*lBt*(1-rP) + sP*sBt*(1-rP)
// actual = (1-aP)*correct + aP*aB
// Without transmission:
// correct = (1-sP)*lBr + sP*sBr
// actual = (1-aP)*correct + aP*aB

/*
// A Ray hits a RendererShape.
// 1: Pick a micronormal. That normal determines everything that happens afterwards. The micronormal can't differ from the macronormal
// by more than 90 degrees. The micronormal is either identical to the macro normal, or uniformly (?) picked on a hemisphere around the macronormal.
// 2: Given that micronormal, determine if there is a transmission or reflection.
// 3a: If there is a reflection, reflect according to the reflection law.
// 3b: If there is a transmission, transmit according to Snell's law.
// Since the micronormal differs from the macronormal, strange things may happen.
// Reflection may happen so that the outgoing ray hits the surface again. If so: repeat the whole process until it escapes the surface (by reflection or transmission). This case is easy.
// Transmission might happen so that the ray is moving back the same way it came. If so, the process must be repeated, but with refractive_index_1 and refractive_index_2 inverted. This case is thus somewhat trickier.

#[allow(dead_code)]
pub fn random_from_brdf_old(incoming_direction: [f64; 3], normal: [f64; 3], material: Material, refractive_index_1: f64, refractive_index_2: f64, pcg: &mut Pcg32) -> ([f64; 3], bool) {
	let mut r = pcg.gen::<f64>();
	let (t_1, t_2) = make_basis_vectors(normal);	
	let micro_normal = if r < material.specular_probability {	
		let mut candidate_micro_normal = random_specular_on_hemisphere(normal, t_1, t_2, material.maximum_specular_angle/2.0, pcg);
		while dot(candidate_micro_normal, incoming_direction) < 0.0 {
			candidate_micro_normal = random_specular_on_hemisphere(normal, t_1, t_2, material.maximum_specular_angle/2.0, pcg);
		}
		candidate_micro_normal
	} else {
		let mut candidate_micro_normal = random_lambertian_on_hemisphere(normal, t_1, t_2, pcg);
		while dot(candidate_micro_normal, incoming_direction) < 0.0 {
			candidate_micro_normal = random_lambertian_on_hemisphere(normal, t_1, t_2, pcg);
		}
		candidate_micro_normal
	};
	// Compute the probability of a reflection.
	r = pcg.gen::<f64>();
	let reflect = material.is_opaque || r < compute_f(incoming_direction, micro_normal, refractive_index_1, refractive_index_2);
	if reflect {
		let outgoing_direction = normalised(sub(mul(2.0*dot(incoming_direction, micro_normal), micro_normal), incoming_direction));
		if dot(outgoing_direction, normal) > 0.0 {
			return (outgoing_direction, true)
		} else {
			return random_from_brdf_0_old(mul(-1.0, outgoing_direction), normal, material, refractive_index_1, refractive_index_2, pcg)
		}
	} else {
		// See https://en.wikipedia.org/wiki/Snell%27s_law.
		let r = refractive_index_1/refractive_index_2;
		let n = if dot(incoming_direction, micro_normal) > 0.0 {
			micro_normal
		} else {
			println!("The micro normal is incorrect.");
			mul(-1.0, micro_normal)
		};
		let l = mul(-1.0, incoming_direction);
		let c = -1.0*dot(n, l);
		let radicand = 1.0-r*r*(1.0-c*c);
		if radicand < 0.0 {
			println!("Negative radicand. Should not happen. The radicand is {}.", radicand);
		}
		let transmission_direction = normalised(add(mul(r, l), mul(r*c-radicand.sqrt(), n)));
		if dot(transmission_direction, normal) < 0.0 {
			return (transmission_direction, false)
		} else {
			return random_from_brdf_0_old(mul(-1.0, transmission_direction), mul(-1.0, normal), material, refractive_index_2, refractive_index_1, pcg)
		}
	}
}

// A Ray hits a RendererShape.
// First, determine if there is a transmission or reflection.
// If there is a reflection, determine if it's specular or diffuse and then reflect in a random direction
// chosen from a diffuse or specular distribution.
// If there is a transmission, we have to be careful, since not all directions are possible if n_1 < n_2.
// First determine if we get a specular transmission or a diffuse one.
// If we get a diffuse one and n_1 > n_2 then we just pick a random transmission direction from a hemisphere.
// If not, we have two cases
// 1: we have a diffuse transmission but n_1 < n_2
// 2: we have a specular transmission
// In either case, pick a "mock"/"simulated" income direction, and use that to compute the transmission direction
// (using Snell's law).
#[allow(dead_code)]
pub fn random_from_brdf_old_2(incoming_direction: [f64; 3], normal: [f64; 3], material: Material, refractive_index_1: f64, refractive_index_2: f64, pcg: &mut Pcg32) -> ([f64; 3], bool) {
	// Compute the probability of a reflection.
	let ratio_reflected = if material.is_opaque {
		1.0
	} else {
		compute_f(incoming_direction, normal, refractive_index_1, refractive_index_2)
	};
	let mut r = pcg.gen::<f64>();
	// Reflection?
	if r < ratio_reflected {
		r = pcg.gen::<f64>();
		if r < material.specular_probability {
			let mirror_direction = normalised(sub(mul(2.0*dot(incoming_direction, normal), normal), incoming_direction));
			let (t_1, t_2) = make_basis_vectors(mirror_direction);	
			return (random_specular_on_hemisphere(mirror_direction, t_1, t_2, material.maximum_specular_angle, pcg), true)
		} else {
			let (t_1, t_2) = make_basis_vectors(normal);	
			return (random_lambertian_on_hemisphere(normal, t_1, t_2, pcg), true)
		}
	} else {
		r = pcg.gen::<f64>();
		if r >= material.specular_probability && refractive_index_1 > refractive_index_2 {
			let (t_1, t_2) = make_basis_vectors(normal);	
			return (random_lambertian_on_hemisphere(mul(-1.0, normal), t_1, t_2, pcg), false)
		}
		let (simulated_incoming_direction, specular_reflection) = if r < material.specular_probability {
			(incoming_direction, true)
		} else {
			let (t_1, t_2) = make_basis_vectors(normal);	
			(random_lambertian_on_hemisphere(normal, t_1, t_2, pcg), false)
		};
		// See https://en.wikipedia.org/wiki/Snell%27s_law.
		let r = refractive_index_1/refractive_index_2;
		let n = if dot(simulated_incoming_direction, normal) > 0.0 {
			normal
		} else {
			println!("The normal is incorrect.");
			mul(-1.0, normal)
		};
		if material.maximum_specular_angle == 0.0 && norm(sub(n, normal)) > 1.0e-6 && specular_reflection {
			println!("The mirror reflection is not working.");
		}
		let l = mul(-1.0, simulated_incoming_direction);
		let c = -1.0*dot(n, l);
		let radicand = 1.0-r*r*(1.0-c*c);
		if radicand < 0.0 {
			println!("Negative radicand. Should not happen. The radicand is {}.", radicand);
		}
		return (normalised(add(mul(r, l), mul(r*c-radicand.sqrt(), n))), false)		
	}
}

// A Ray hits a RendererShape.
// First, determine if the interaction is specular or diffuse.
// Then, pick a random specular or lambertian reflection direction.
// Then, if the object is opaque, the interaction will be a reflection. Return the random reflection direction which was picked.
// If it's not opaque, compute which micro normal is needed to get the randomly chosen reflection direction.
// For that normal and incoming direction, compute the probability of a reflection, and then either reflect or transmit the ray.
#[allow(dead_code)]
pub fn random_from_brdf_old_3(incoming_direction: [f64; 3], normal: [f64; 3], material: Material, refractive_index_1: f64, refractive_index_2: f64, pcg: &mut Pcg32) -> ([f64; 3], bool) {
	let mut r = pcg.gen::<f64>();
	let (outgoing_direction, specular_reflection) = if r < material.specular_probability {
		let mirror_direction = normalised(sub(mul(2.0*dot(incoming_direction, normal), normal), incoming_direction));
		let (t_1, t_2) = make_basis_vectors(mirror_direction);	
		(random_specular_on_hemisphere(mirror_direction, t_1, t_2, material.maximum_specular_angle, pcg), true)
	} else {
		let (t_1, t_2) = make_basis_vectors(normal);	
		(random_lambertian_on_hemisphere(normal, t_1, t_2, pcg), false)
	};
	if material.is_opaque {
		return (outgoing_direction, true)
	}
	let micro_normal = normalised(add(incoming_direction, outgoing_direction));
	let ratio_reflected = compute_f(incoming_direction, micro_normal, refractive_index_1, refractive_index_2);
	r = pcg.gen::<f64>();
	if r < ratio_reflected {
		(outgoing_direction, true)
	// @TODO: Should it be micro_normal instead of n somewhere in this code?
	} else {
		// See https://en.wikipedia.org/wiki/Snell%27s_law.
		let r = refractive_index_1/refractive_index_2;
		let n = if dot(incoming_direction, micro_normal) > 0.0 {
			micro_normal
		} else {
			println!("The micro normal is incorrect.");
			mul(-1.0, micro_normal)
		};
		if material.maximum_specular_angle == 0.0 && norm(sub(n, normal)) > 1.0e-6 && specular_reflection {
			println!("The mirror reflection is not working.");
		}
		let l = mul(-1.0, incoming_direction);
		let c = -1.0*dot(n, l);
		let radicand = 1.0-r*r*(1.0-c*c);
		if radicand < 0.0 {
			println!("Negative radicand. Should not happen. The radicand is {}.", radicand);
		}
		(normalised(add(mul(r, l), mul(r*c-radicand.sqrt(), n))), false)
	}
}
*/
// COMPUTE BRDFS

// Compute the brdf both for a diffuse interaction and for a specular interaction and compute the weighted average.
#[allow(dead_code)]
pub fn brdf(incoming_direction: [f64; 3], outgoing_direction: [f64; 3], normal: [f64; 3], refractive_index_1: f64, refractive_index_2: f64, material: Material) -> f64 {
	let brdf_reflection_specular = brdf_specular(incoming_direction, outgoing_direction, normal, material);
	let brdf_reflection_lambertian = brdf_lambertian(outgoing_direction, normal);
	let brdf_reflection = material.specular_probability*brdf_reflection_specular + (1.0-material.specular_probability)*brdf_reflection_lambertian;
	if material.is_opaque {
		brdf_reflection
	} else {
		let must_reflect_to_hit_light = dot(outgoing_direction, normal) > 0.0;
		let ratio_reflected = compute_f(incoming_direction, normal, refractive_index_1, refractive_index_2);
		if must_reflect_to_hit_light {
			//println!("reflect, {}", ratio_reflected*brdf_reflection);
			ratio_reflected*brdf_reflection
		} else {
			// See https://en.wikipedia.org/wiki/Snell%27s_law.
			let r = refractive_index_1/refractive_index_2;
			let l = mul(-1.0, incoming_direction);
			let c = -1.0*dot(normal, l);
			let mut radicand = 1.0-r*r*(1.0-c*c);
			if radicand < 0.0 {
				//println!("Negative radicand. Should not happen. The radicand is {}.", radicand);
				radicand = 0.0;
			}
			let transmitted_direction = normalised(add(mul(r, l), mul(r*c-radicand.sqrt(), normal)));
			let brdf_transmission_specular = brdf_specular_transmission(transmitted_direction, outgoing_direction, mul(-1.0, normal), material);
			// @TODO: Fix the lambertian transmission.
			let brdf_transmission_lambertian = brdf_lambertian(transmitted_direction, mul(-1.0, normal));
			//let brdf_transmission_lambertian = 0.0;
			let brdf_transmission = material.specular_probability*brdf_transmission_specular + (1.0-material.specular_probability)*brdf_transmission_lambertian;
			//println!("transmit, {}", (1.0-ratio_reflected)*brd_transmission);
			(1.0-ratio_reflected)*brdf_transmission
		}
	}
}

#[allow(dead_code)]
fn brdf_specular_transmission(specular_transmitted_direction: [f64; 3], outgoing_direction: [f64; 3], normal: [f64; 3], material: Material) -> f64 {
	let phi = dot(outgoing_direction, specular_transmitted_direction).acos();
	let maximum_angle = min(PI/2.0 - dot(normal, specular_transmitted_direction), material.maximum_specular_angle);
	if phi < maximum_angle {
		let a = 1.0/(2.0*PI*(maximum_angle-maximum_angle.sin()));
		2.0*PI*a*(maximum_angle-phi)
	} else {
		0.0
	}
}

#[allow(dead_code)]
fn brdf_specular(incoming_direction: [f64; 3], outgoing_direction: [f64; 3], normal: [f64; 3], material: Material) -> f64 {
	let specular_direction = sub(mul(2.0*dot(incoming_direction, normal), normal), incoming_direction);
	let phi = dot(outgoing_direction, specular_direction).acos();
	let maximum_angle = min(PI/2.0 - dot(normal, specular_direction), material.maximum_specular_angle);
	if phi < maximum_angle {
		let a = 1.0/(2.0*PI*(maximum_angle-maximum_angle.sin()));
		2.0*PI*a*(maximum_angle-phi)
	} else {
		0.0
	}
}

#[allow(dead_code)]
fn brdf_lambertian(direction_1: [f64; 3], direction_2: [f64; 3]) -> f64 {
	2.0*dot(direction_1, direction_2)
}

#[allow(dead_code)]
fn brdf_lambertian_in_cone(direction_1: [f64; 3], direction_2: [f64; 3], max_angle: f64) -> f64 {
	let angle = dot(direction_1, direction_2).acos();
	if angle >= max_angle {
		0.0
	} else {
		let scaled_angle = angle*(PI/2.0)/max_angle;
		2.0*(PI/2.0)/max_angle*(scaled_angle.cos())
	}
}

// OTHER MATH

#[allow(dead_code)]
pub fn min(left: f64, right: f64) -> f64 {
	if left < right {
		left
	} else {
		right
	}
}

#[allow(dead_code)]
// See http://blackpawn.com/texts/pointinpoly/.
pub fn point_in_triangle(point: [f64; 3], first_node: [f64; 3], second_node: [f64; 3], third_node: [f64; 3]) -> bool {
	same_side(point, first_node, second_node, third_node) &&
	same_side(point, second_node, first_node, third_node) &&
	same_side(point, third_node, first_node, second_node)
}

#[inline(always)]
#[allow(dead_code)]
pub fn same_side(test_point: [f64; 3], point_inside: [f64; 3], first_node: [f64; 3], second_node: [f64; 3]) -> bool {
	dot(cross(sub(second_node, first_node), sub(test_point, first_node)), cross(sub(second_node, first_node), sub(point_inside, first_node))) >= 0.0
}

// See https://en.wikipedia.org/wiki/Fresnel_equations.
#[allow(dead_code)]
pub fn compute_f(wi: [f64; 3], n: [f64; 3], eta_1: f64, eta_2: f64) -> f64 {
	let cos_theta_i = dot(wi, n);
	let cos_theta_i_square = cos_theta_i*cos_theta_i;
	let sin_theta_i_square = 1.0-cos_theta_i_square;

	// Total internal reflection.
	if sin_theta_i_square.sqrt() >= eta_2/eta_1 {
		return 1.0;
	}

	let sin_theta_t_square = (eta_1/eta_2)*(eta_1/eta_2)*sin_theta_i_square;
	let cos_theta_t = (1.0-sin_theta_t_square).sqrt();

	let eta_1_cos_theta_i = eta_1*cos_theta_i;
	let eta_2_cos_theta_t = eta_2*cos_theta_t;
	let r_transverse_sqrt = (eta_1_cos_theta_i - eta_2_cos_theta_t)/(eta_1_cos_theta_i + eta_2_cos_theta_t);
	let r_transverse = r_transverse_sqrt*r_transverse_sqrt;

	let eta_2_cos_theta_i = eta_2*cos_theta_i;
	let eta_1_cos_theta_t = eta_1*cos_theta_t;
	let r_parallell_sqrt = (eta_2_cos_theta_i - eta_1_cos_theta_t)/(eta_2_cos_theta_i + eta_1_cos_theta_t);
	let r_parallell = r_parallell_sqrt*r_parallell_sqrt;

	(r_transverse+r_parallell)/2.0
}

#[allow(dead_code)]
pub fn intensity_to_color(color: [f64; 3]) -> [f64; 3] {
	let mut max_intensity = color[0];
	if color[1] > max_intensity {
		max_intensity = color[1];
	}
	if color[2] > max_intensity {
		max_intensity = color[2];
	}
	let factor = if max_intensity > 1.0e-12 {
		let modified_max_intensity = (GAMMA*max_intensity).atan()*255.0/(PI/2.0);
		modified_max_intensity/max_intensity
	} else {
		1.0
	};
	mul(factor, color)
}

/*

pub fn random_uniform_on_hemisphere(normal: [f64; 3], pcg: &mut Pcg32) -> [f64; 3] {
	let r1 = pcg.gen::<f64>();
	let r2 = pcg.gen::<f64>();
    let phi = 2.0*PI*r1;
    let p = (1.0-r2*r2).sqrt();
	let (t3, t1, t2) = make_basis_vectors(normal);
	// Might not be perfectly normalised due to rounding errors. Thus, normalise!
	normalised(add(add(mul(p*phi.cos(), t1), mul(p*phi.sin(), t2)), mul(r2, t3)))
}

#[allow(dead_code)]
pub fn random_lambertian_on_stripe(normal: [f64; 3], pcg: &mut Pcg32, min_theta: f64, max_theta: f64) -> [f64; 3] {
	let r1_sqrt_min = min_theta.cos();
	let r1_sqrt_max = max_theta.cos();
	let r1_min = r1_sqrt_min*r1_sqrt_min;
	let r1_max = r1_sqrt_max*r1_sqrt_max;
	let mut r1 = pcg.gen::<f64>();
	r1 *= r1_max-r1_min;
	r1 += r1_min;
	let r2 = pcg.gen::<f64>();
	let sintheta = (1.0-r1).sqrt();
	let costheta = r1.sqrt();
	let phi = 2.0*PI*r2;
	let (t3, t1, t2) = make_basis_vectors(normal);
	// Might not be perfectly normalised due to rounding errors. Thus, normalise!
	normalised(add(add(mul(sintheta*phi.cos(), t1), mul(sintheta*phi.sin(), t2)), mul(costheta, t3)))
}

#[allow(dead_code)]
pub fn random_uniform_on_stripe(normal: [f64; 3], pcg: &mut Pcg32, min_theta: f64, max_theta: f64) -> [f64; 3] {
	let r2_sqrt_min = max_theta.cos();
	let r2_sqrt_max = min_theta.cos();
	let r2_min = r2_sqrt_min*r2_sqrt_min;
	let r2_max = r2_sqrt_max*r2_sqrt_max;
	let r1 = pcg.gen::<f64>();
	let mut r2 = pcg.gen::<f64>();
	r2 *= r2_max-r2_min;
	r2 += r2_min;
    let phi = 2.0*PI*r1;
    let p = (1.0-r2*r2).sqrt();
	let (t3, t1, t2) = make_basis_vectors(normal);
	// Might not be perfectly normalised due to rounding errors. Thus, normalise!
	normalised(add(add(mul(p*phi.cos(), t1), mul(p*phi.sin(), t2)), mul(r2, t3)))
}

#[allow(dead_code)]
pub fn brdf_ct(wi: [f64; 3], wo: [f64; 3], n: [f64; 3], a: f64, eta_2: f64, ks: f64) -> f64 {
	let eta_1 = 1.0;
	let kd = 1.0 - ks;
	let won = dot(wo, n);
	let win = dot(wi, n);
	let h = normalised(add(wi, wo));
	let fd = 2.0;
	//let fs = compute_d(wo, n, h, a)*compute_f(wi, h, n, eta_1, eta_2)*compute_g(wi, wo, h, n, a)/(won*win);
	//let fs = compute_f(wi, h, n, eta_1, eta_2);
	//(kd*fd + ks*fs)*won
	2.0*won
}

#[allow(dead_code)]
pub fn compute_d(wo: [f64; 3], n: [f64; 3], h: [f64; 3], a: f64) -> f64 {
	let cos_theta_h = dot(h, n);
	let theta_h = cos_theta_h.acos();
	let tan_theta_h = theta_h.tan();
	let b = a*a+tan_theta_h*tan_theta_h;
	let numerator = a*a*compute_xi(theta_h);
	let denominator = PI*cos_theta_h*cos_theta_h*cos_theta_h*cos_theta_h*b*b;
	numerator/denominator
}

#[allow(dead_code)]
pub fn compute_g(wi: [f64; 3], wo: [f64; 3], h: [f64; 3], n: [f64; 3], a: f64) -> f64 {
	let gpwo = compute_gp(wo, h, n, a);
	let gpwi = compute_gp(wi, h, n, a);
	gpwo*gpwi
}

#[allow(dead_code)]
pub fn compute_gp(w: [f64; 3], n: [f64; 3], h: [f64; 3], a: f64) -> f64 {
	let wh = dot(w, h);
	let wn = dot(w, n);
	let xi = compute_xi(wh/wn);
	let theta_w = wn.acos();
	let tan_theta_w = theta_w.tan();
	xi*2.0/(1.0+(a*a*tan_theta_w*tan_theta_w).sqrt())
}

#[allow(dead_code)]
pub fn compute_xi(x: f64) -> f64 {
	if x > 0.0 {
		return 1.0;
	}
	0.0
}
*/
