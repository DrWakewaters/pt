use std::f64;

use serde_derive::{Serialize, Deserialize};

use crate::material::Material;
use crate::math::{add, dot, distance, mul, norm_squared, normalised, sub};
use crate::physicscylinder::PhysicsCylinder;
use crate::ray::Ray;
use crate::renderershape::RendererShape;

#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
pub struct RendererCylinder {
	position: [f64; 3],
    direction: [f64; 3],
    length: f64,
	radius: f64,
	material: Material,
    pub active: bool
}

impl RendererCylinder {
	pub fn new(physics_cylinder: &PhysicsCylinder) -> Self {
		Self {
			position: physics_cylinder.position,
            direction: physics_cylinder.direction,
            length: physics_cylinder.length,
			radius: physics_cylinder.radius,
			material: physics_cylinder.material,
            active: physics_cylinder.active,
        }
	}

    fn distance_to_mantle(&self, ray: &Ray) -> f64 {
        let a = self.radius*self.radius*norm_squared(self.direction);
        let b = sub(ray.position, self.position);
        let c = sub(b, self.direction);
        let d = ray.direction;
        let e = b[0];
        let f = b[1];
        let g = b[2];
        let j = c[0];
        let k = c[1];
        let l = c[2];
        let m = d[0];
        let n = d[1];
        let o = d[2];

        let t_squared_terms = f*f*o*o - 2.0*f*g*n*o - 2.0*f*k*o*o + g*g*n*n - 2.0*g*l*n*n + 2.0*g*k*n*o + k*k*o*o + l*l*n*n - 2.0*k*l*n*o + g*g*m*m - 2.0*g*l*m*m + 2.0*g*j*m*o - 2.0*e*g*m*o + j*j*o*o - 2.0*j*l*m*o - 2.0*e*j*o*o + l*l*m*m + 2.0*e*l*m*o + e*e*o*o + f*f*m*m + 2.0*f*j*m*n - 2.0*f*k*m*m + j*j*n*n - 2.0*e*f*m*n - 2.0*j*k*m*n - 2.0*e*j*n*n + k*k*m*m + 2.0*e*k*m*n + e*e*n*n + 2.0*f*l*n*o;
        let t_terms = 2.0*f*f*l*o - 2.0*f*g*k*o - 2.0*f*g*l*n - 2.0*f*k*l*o + 2.0*f*l*l*n + 2.0*g*g*k*n + 2.0*g*k*k*o - 2.0*g*k*l*n + 2.0*g*g*j*m + 2.0*g*j*j*o - 2.0*g*j*l*m - 2.0*e*g*j*o - 2.0*e*g*l*m - 2.0*e*j*l*o + 2.0*e*l*l*m + 2.0*e*e*l*o + 2.0*f*f*j*m + 2.0*f*j*j*n - 2.0*f*j*k*m - 2.0*e*f*j*n - 2.0*e*f*k*m - 2.0*e*j*k*n + 2.0*e*k*k*m + 2.0*e*e*k*n;
        let constant_terms = f*f*l*l - 2.0*f*g*k*l + g*g*k*k + g*g*j*j - 2.0*e*g*j*l + e*e*l*l + f*f*j*j - 2.0*e*f*j*k + e*e*k*k - a;
        let p = t_terms/t_squared_terms;
        let q = constant_terms/t_squared_terms;
        let discriminant = (p/2.0)*(p/2.0) - q;
        if discriminant < 1.0e-6 {
            return f64::MAX;
        }
        let t_1 = -0.5*p - discriminant.sqrt();
        let t_2 = -0.5*p + discriminant.sqrt();
        if t_2 < 1.0e-6 {
            return f64::MAX;
        }
        let distance_to_infinite_cylinder =  if t_1 < 1.0e-6 {
            t_2
        } else {
            t_1
        };
        let point_on_infinite_cylinder = add(ray.position, mul(distance_to_infinite_cylinder, ray.direction));
        let start_point = self.position;
        let end_point = add(self.position, mul(self.length, self.direction));
        if dot(sub(point_on_infinite_cylinder, start_point), sub(end_point, start_point)) > 0.0 && dot(sub(point_on_infinite_cylinder, end_point), sub(start_point, end_point)) > 0.0 {
            return distance_to_infinite_cylinder;
        }
        f64::MAX
	}

    fn distance_to_cap(&self, ray: &Ray) -> f64 {
        let start_point = self.position;
        let end_point = add(self.position, mul(self.length, self.direction));
        let distance_to_start_cap = self.distance_to_cap_inner(ray, self.direction, start_point, self.radius);
        let distance_to_end_cap = self.distance_to_cap_inner(ray, self.direction, end_point, self.radius);
        if distance_to_start_cap < distance_to_end_cap  {
            return distance_to_start_cap;
        }
        distance_to_end_cap
    }

    fn distance_to_cap_inner(&self, ray: &Ray, normal_on_cap: [f64; 3], center_point_on_cap: [f64; 3], radius: f64) -> f64 {
        let q = dot(normal_on_cap, sub(center_point_on_cap, ray.position))/dot(normal_on_cap, ray.direction);
        if q < 1.0e-6 {
            return f64::MAX;
        }
        let w = add(ray.position, mul(q, ray.direction));
        if distance(w, center_point_on_cap) < radius {
            return q;
        }
        f64::MAX
    }
}

impl RendererShape for RendererCylinder {
	fn distance(&self, ray: &Ray) -> f64 {
        let mantle_distance = self.distance_to_mantle(ray);
        return mantle_distance;
        /*
        let cap_distance = self.distance_to_cap(ray);
        if mantle_distance < cap_distance {
            return mantle_distance;
        }
        cap_distance
        */
    }

	fn material(&self) -> Material {
		self.material
	}

	fn normal(&self, point: [f64; 3]) -> [f64; 3] {
        /*
        let start_point = self.position;
        let end_point = add(self.position, mul(self.length, self.direction));
        if distance(start_point, point) < self.radius {
            return mul(-1.0, self.direction);
        } else if distance(end_point, point) < self.radius {
            return self.direction;
        }
        */
        normalised(sub(point, add(self.position, mul(dot(self.direction, sub(point, self.position)), self.direction))))
	}
}