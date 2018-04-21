use std::f64::consts::PI;

use rand::Rng;
use pcg_rand::Pcg32;

#[inline(always)]
pub fn add(left: [f64; 3], right: [f64; 3]) -> [f64; 3] {
	[left[0]+right[0], left[1]+right[1], left[2]+right[2]]
}

#[inline(always)]
pub fn add_f32(left: [f32; 3], right: [f32; 3]) -> [f32; 3] {
	[left[0]+right[0], left[1]+right[1], left[2]+right[2]]
}

#[inline(always)]
pub fn add_u16(left: [u16; 3], right: [u16; 3]) -> [u16; 3] {
	[left[0]+right[0], left[1]+right[1], left[2]+right[2]]
}

#[inline(always)]
pub fn sub(left: [f64; 3], right: [f64; 3]) -> [f64; 3] {
	[left[0]-right[0], left[1]-right[1], left[2]-right[2]]
}

#[inline(always)]
pub fn mul(scalar: f64, vector: [f64; 3]) -> [f64; 3] {
	[scalar*vector[0], scalar*vector[1], scalar*vector[2]]
}

#[inline(always)]
pub fn elementwise_mul(left: [f64; 3], right: [f64; 3]) -> [f64; 3] {
[left[0]*right[0], left[1]*right[1], left[2]*right[2]]
}

#[inline(always)]
pub fn dot(left: [f64; 3], right: [f64; 3]) -> f64 {
	left[0]*right[0]+left[1]*right[1]+left[2]*right[2]
}

#[inline(always)]
pub fn cross(left: [f64; 3], right: [f64; 3]) -> [f64; 3] {
	[left[1]*right[2]-left[2]*right[1], left[2]*right[0]-left[0]*right[2], left[0]*right[1]-left[1]*right[0]]
}

#[inline(always)]
pub fn norm(vector: [f64; 3]) -> f64 {
	(vector[0]*vector[0]+vector[1]*vector[1]+vector[2]*vector[2]).sqrt()
}

#[inline(always)]
pub fn distance(left: [f64; 3], right: [f64; 3]) -> f64 {
	let vector = sub(left, right);
	(vector[0]*vector[0]+vector[1]*vector[1]+vector[2]*vector[2]).sqrt()
}

#[inline(always)]
pub fn distance_squared(left: [f64; 3], right: [f64; 3]) -> f64 {
	let vector = sub(left, right);
	(vector[0]*vector[0]+vector[1]*vector[1]+vector[2]*vector[2])
}

#[inline(always)]
pub fn normalised(vector: [f64; 3]) -> [f64; 3] {
	let norm = (vector[0]*vector[0]+vector[1]*vector[1]+vector[2]*vector[2]).sqrt();
	[vector[0]/norm, vector[1]/norm, vector[2]/norm]
}

#[inline(always)]
pub fn same_side(test_point: [f64; 3], point_inside: [f64; 3], first_node: [f64; 3], second_node: [f64; 3]) -> bool {
	dot(cross(sub(second_node, first_node), sub(test_point, first_node)), cross(sub(second_node, first_node), sub(point_inside, first_node))) >= 0.0
}

#[inline(always)]
pub fn point_in_triangle(point: [f64; 3], first_node: [f64; 3], second_node: [f64; 3], third_node: [f64; 3]) -> bool {
	if same_side(point, first_node, second_node, third_node) {
		return false;
	}
	if same_side(point, second_node, first_node, third_node) {
		return false;
	}
	if same_side(point, third_node, first_node, second_node) {
		return false;
	}
	return true;
}

#[inline(always)]
pub fn make_basis_vectors(v1: [f64; 3]) -> ([f64; 3], [f64; 3], [f64; 3]) {
	let mut v2 = [0.0, 0.0, 1.0];
	let mut v3 = [0.0, 1.0, 0.0];
	let u1 = v1;
	let u1_normalised = normalised(u1);
	let mut u2 = sub(v2, mul(dot(u1_normalised, v2), u1_normalised));
	if norm(u2) < 1e-6 {
		v2 = [1.0, 0.0, 0.0];
		u2 = sub(v2, mul(dot(u1_normalised, v2), u1_normalised));
	}
	let u2_normalised = normalised(u2);
	let mut u3 = sub(sub(v3, mul(dot(u1_normalised, v3), u1_normalised)), mul(dot(u2_normalised, v3), u2_normalised));
	if norm(u3) < 1e-6 {
		v3 = [1.0, 0.0, 0.0];
		u3 = sub(sub(v3, mul(dot(u1_normalised, v3), u1_normalised)), mul(dot(u2_normalised, v3), u2_normalised));
	}
	let u3_normalised = normalised(u3);
	(u1_normalised, u2_normalised, u3_normalised)
}

#[inline(always)]
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
pub fn matrix_vector_mul(matrix: [[f64; 3]; 3] , vector: [f64; 3]) -> [f64; 3] {
	[matrix[0][0]*vector[0] + matrix[1][0]*vector[1] + matrix[2][0]*vector[2],
	matrix[0][1]*vector[0] + matrix[1][1]*vector[1] + matrix[2][1]*vector[2],
	matrix[0][2]*vector[0] + matrix[1][2]*vector[1] + matrix[2][2]*vector[2]]
}

// The solution to matrix*vector' = vector is vector' = matrix^(-1)*vector.
#[inline(always)]
pub fn solve_linear(matrix: [[f64; 3]; 3] , vector: [f64; 3]) -> [f64; 3] {
	let matrix_inverse = invert(matrix);
	matrix_vector_mul(matrix_inverse, vector)
}

#[inline(always)]
pub fn quarternion_inverse(quarternion: [f64; 4]) -> [f64; 4] {
	[quarternion[0], -1.0*quarternion[1], -1.0*quarternion[2], -1.0*quarternion[3]]
}

#[inline(always)]
pub fn quarternion_product(left: [f64; 4], right: [f64; 4]) -> [f64; 4] {
	[left[0]*right[0] - left[1]*right[1] - left[2]*right[2] - left[3]*right[3],
	 left[0]*right[1] + left[1]*right[0] + left[2]*right[3] - left[3]*right[2],
	 left[0]*right[2] - left[1]*right[3] + left[2]*right[0] + left[3]*right[1],
	 left[0]*right[3] + left[1]*right[2] - left[2]*right[1] + left[3]*right[0]]
}

#[inline(always)]
pub fn pick_reflection_uniform(normal: [f64; 3], pcg: &mut Pcg32) -> [f64; 3] {
	let r1_u32: u32 = pcg.gen();
	let r2_u32: u32 = pcg.gen();
	let r1 = (r1_u32 as f64)/(<u32>::max_value() as f64);
	let r2 = (r2_u32 as f64)/(<u32>::max_value() as f64);
    let phi = 2.0*PI*r1;
    let p = (1.0-r2*r2).sqrt();
	let (t3, t1, t2) = make_basis_vectors(normal);
	normalised(add(add(mul(p*phi.cos(), t1), mul(p*phi.sin(), t2)), mul(r2, t3)))
}

#[inline(always)]
pub fn pick_reflection_uniform_on_stripe(normal: [f64; 3], pcg: &mut Pcg32, min_theta: f64, max_theta: f64) -> [f64; 3] {
	let r2_sqrt_min = max_theta.cos();
	let r2_sqrt_max = min_theta.cos();
	let r2_min = r2_sqrt_min*r2_sqrt_min;
	let r2_max = r2_sqrt_max*r2_sqrt_max;
	let r1_u32: u32 = pcg.gen();
	let r2_u32: u32 = pcg.gen();
	let r1 = (r1_u32 as f64)/(<u32>::max_value() as f64);
	let mut r2 = (r2_u32 as f64)/(<u32>::max_value() as f64);
	r2 *= r2_max-r2_min;
	r2 += r2_min;
    let phi = 2.0*PI*r1;
    let p = (1.0-r2*r2).sqrt();
	let (t3, t1, t2) = make_basis_vectors(normal);
	normalised(add(add(mul(p*phi.cos(), t1), mul(p*phi.sin(), t2)), mul(r2, t3)))
}

#[inline(always)]
pub fn pick_reflection_lambertian(normal: [f64; 3], pcg: &mut Pcg32) -> [f64; 3] {
	let r1_u32: u32 = pcg.gen();
	let r2_u32: u32 = pcg.gen();
	let r1 = (r1_u32 as f64)/(<u32>::max_value() as f64);
	let r2 = (r2_u32 as f64)/(<u32>::max_value() as f64);
	let sintheta = (1.0-r1).sqrt();
	let costheta = r1.sqrt();
	let phi = 2.0*PI*r2;
	let (t3, t1, t2) = make_basis_vectors(normal);
	normalised(add(add(mul(sintheta*phi.cos(), t1), mul(sintheta*phi.sin(), t2)), mul(costheta, t3)))
}

#[inline(always)]
pub fn pick_reflection_lambertian_on_stripe(normal: [f64; 3], pcg: &mut Pcg32, min_theta: f64, max_theta: f64) -> [f64; 3] {
	let r1_sqrt_min = min_theta.cos();
	let r1_sqrt_max = max_theta.cos();
	let r1_min = r1_sqrt_min*r1_sqrt_min;
	let r1_max = r1_sqrt_max*r1_sqrt_max;
	let r1_u32: u32 = pcg.gen();
	let r2_u32: u32 = pcg.gen();
	let mut r1 = (r1_u32 as f64)/(<u32>::max_value() as f64);
	r1 *= r1_max-r1_min;
	r1 += r1_min;
	let r2 = (r2_u32 as f64)/(<u32>::max_value() as f64);
	let sintheta = (1.0-r1).sqrt();
	let costheta = r1.sqrt();
	let phi = 2.0*PI*r2;
	let (t3, t1, t2) = make_basis_vectors(normal);
	normalised(add(add(mul(sintheta*phi.cos(), t1), mul(sintheta*phi.sin(), t2)), mul(costheta, t3)))
}

#[inline(always)]
pub fn dfg(wi: [f64; 3], wo: [f64; 3], n: [f64; 3], a: f64, eta_2: f64, ks: f64) -> f64 {
	let eta_1 = 1.0;
	let kd = 1.0 - ks;
	let won = dot(wo, n);
	let win = dot(wi, n);
	let h = normalised(add(wi, wo));
	let fd = 2.0;
	//let fs = compute_d(wo, n, h, a)*compute_f(wi, h, n, eta_1, eta_2)*compute_g(wi, wo, h, n, a)/(won*win);
	let fs = compute_f(wi, h, n, eta_1, eta_2);
	(kd*fd + ks*fs)*won
}

#[inline(always)]
pub fn compute_d(wo: [f64; 3], n: [f64; 3], h: [f64; 3], a: f64) -> f64 {
	let cos_theta_h = dot(h, n);
	let theta_h = cos_theta_h.acos();
	let tan_theta_h = theta_h.tan();
	let b = a*a+tan_theta_h*tan_theta_h;
	let numerator = a*a*compute_xi(theta_h);
	let denominator = PI*cos_theta_h*cos_theta_h*cos_theta_h*cos_theta_h*b*b;
	numerator/denominator
}

#[inline(always)]
pub fn compute_f(wi: [f64; 3], h: [f64; 3], n: [f64; 3], eta_1: f64, eta_2: f64) -> f64 {
	let cos_theta_i = dot(wi, n);
	let cos_theta_i_square = cos_theta_i*cos_theta_i;
	let sin_theta_t_square = (eta_1/eta_2)*(eta_1/eta_2)*(1.0-cos_theta_i_square);
	let cos_theta_t = 1.0-sin_theta_t_square;

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

#[inline(always)]
pub fn compute_g(wi: [f64; 3], wo: [f64; 3], h: [f64; 3], n: [f64; 3], a: f64) -> f64 {
	return 1.0;
	let gpwo = compute_gp(wo, h, n, a);
	let gpwi = compute_gp(wi, h, n, a);
	gpwo*gpwi
}

#[inline(always)]
pub fn compute_gp(w: [f64; 3], n: [f64; 3], h: [f64; 3], a: f64) -> f64 {
	let wh = dot(w, h);
	let wn = dot(w, n);
	let xi = compute_xi(wh/wn);
	let theta_w = wn.acos();
	let tan_theta_w = theta_w.tan();
	xi*2.0/(1.0+(a*a*tan_theta_w*tan_theta_w).sqrt())
}

#[inline(always)]
pub fn compute_xi(x: f64) -> f64 {
	if x > 0.0 {
		return 1.0;
	}
	return 0.0;
}
