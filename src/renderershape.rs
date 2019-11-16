use crate::material::Material;
use crate::ray::Ray;

pub trait RendererShape {
    fn distance(&self, ray: &Ray) -> f64;
    fn material(&self) -> Material;
    fn normal(&self, point: [f64; 3]) -> [f64; 3];
}
