use std::f64;

use serde_derive::{Serialize, Deserialize};

use crate::material::Material;
use crate::math::{cross, dot, sub};
use crate::physicstriangle::PhysicsTriangle;
use crate::ray::Ray;
use crate::renderershape::RendererShape;

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct RendererTriangle {
	node0: [f64; 3],
	node1: [f64; 3],
	node2: [f64; 3],
	e1: [f64; 3],
	e2: [f64; 3],
	normal: [f64; 3],
	material: Material,
}

impl RendererTriangle {
	pub fn new(physics_triangle: &PhysicsTriangle) -> Self {
		Self {
			node0: physics_triangle.node0,
			node1: physics_triangle.node1,
			node2: physics_triangle.node2,
			e1: physics_triangle.e1,
			e2: physics_triangle.e2,
			normal: physics_triangle.normal,
			material: physics_triangle.material,
		}
	}
}

// See http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.189.5084&rep=rep1&type=pdf.
impl RendererShape for RendererTriangle {
	fn distance(&self, ray: &Ray) -> f64 {
		let h = cross(ray.direction, self.e2);
		let a = dot(self.e1, h);
		if a.abs() < 1.0e-6 {
			return f64::MAX;
		}
		let f = 1.0/a;
		let s = sub(ray.position, self.node0);
		let u = f*dot(s, h);
		if u < 1.0e-6 || u > 1.0-1.0e-6 {
			return f64::MAX;
		}
		let q = cross(s, self.e1);
		let v = f*dot(ray.direction, q);
		if v < 1.0e-6 || u+v > 1.0-1.0e-6 {
			return f64::MAX;
		}
		let d = f*dot(self.e2, q);
		if d < 1.0e-6 {
			return f64::MAX;
		}
		d
	}

	fn material(&self) -> Material {
		self.material
	}

	fn normal(&self, _point: [f64; 3]) -> [f64; 3] {
		self.normal
	}
}

/*
// http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.159.5135&rep=rep1&type=pdf
pub n_u: f64,
pub n_v: f64,
pub n_d: f64,
pub k: u64,
pub b_nu: f64,
pub b_nv: f64,
pub b_d: f64,
pub c_nu: f64,
pub c_nv: f64,
pub c_d: f64,
pub normal: [f64; 3],
pub material: Material,
*/

/*
let b = physics_triangle.e2;
let c = physics_triangle.e1;
let node0 = physics_triangle.node0;
let normal = physics_triangle.normal;
let mut k = 0;
let mut u = 1;
let mut v = 2;
if normal[1].abs() > normal[0].abs() && normal[1].abs() > normal[2].abs() {
	k = 1;
	u = 2;
	v = 0;
}
if normal[2].abs() > normal[0].abs() && normal[2].abs() > normal[1].abs() {
	k = 2;
	u = 0;
	v = 1;
}
let denom = b[u]*c[v]-b[v]*c[u];
Self {
	n_u: normal[u]/normal[k],
	n_v: normal[v]/normal[k],
	n_d: dot(normal, node0)/normal[k],
	k: k as u64,
	b_nu: b[u]/denom,
	b_nv: -1.0*b[v]/denom,
	b_d: (b[v]*node0[u]-b[u]*node0[v])/denom,
	c_nu: -1.0*c[u]/denom,
	c_nv: c[v]/denom,
	c_d: (c[u]*node0[v]-c[v]*node0[u])/denom,
	normal: physics_triangle.normal,
	material: physics_triangle.material,
}
*/

/*
let k = t.k as usize;
let u = MODULO[(t.k+1) as usize];
let v = MODULO[(t.k+2) as usize];
let nd = 1.0/(direction[k]+t.n_u*direction[u]+t.n_v*direction[v]);
let d = (t.n_d-intersection[k]-t.n_u*intersection[u]-t.n_v*intersection[v])*nd;
if d < 1e-9 || d > triangle_distance {
	continue;
}
let hu = intersection[u]+d*direction[u];
let hv = intersection[v]+d*direction[v];
let lambda = hu*t.b_nu+hv*t.b_nv+t.b_d;
if lambda < 0.0 {
	continue;
}
let mue = hu*t.c_nu+hv*t.c_nv+t.c_d;
if mue < 0.0 || lambda+mue > 1.0 {
	continue;
}*/
