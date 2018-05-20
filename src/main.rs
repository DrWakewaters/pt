#![feature(placement_in_syntax)]
#![feature(collection_placement)]
#[macro_use]

extern crate serde_derive;
extern crate serde_json;
extern crate png;
extern crate pcg_rand;
extern crate rand;
extern crate time;

mod camera;
mod collisiontype;
mod datafordrawing;
mod dataforstoring;
mod hitpoint;
mod light;
mod material;
mod math;
mod pathtracer;
mod physics;
mod physicsscene;
mod physicssphere;
mod physicstriangle;
mod pixel;
mod pngfile;
mod ray;
mod renderer;
mod rendereroutput;
mod rendererscene;
mod renderershape;
mod renderersphere;
mod renderertriangle;
mod rhf;
mod tracing;

use pathtracer::Pathtracer;
use tracing::Tracing;

const NUMBER_OF_BINS: usize = 128;

fn main() {
	let mut pathtracer = Pathtracer::new(1000, 1000, 0, 8, 3, false, vec![(0.3, 1)], 5, 1, 50_000_000, true, Tracing::Light);
	pathtracer.render_images_to_png();
}
