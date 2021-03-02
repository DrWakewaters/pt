use crate::material::Material;

pub trait RendererShape {
    fn distance(&self, position: [f64; 3], direction: [f64; 3]) -> f64;
    fn material(&self) -> Material;
    fn normal(&self, point: [f64; 3]) -> [f64; 3];
}
