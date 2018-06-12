#![feature(placement_in_syntax)]
#![feature(collection_placement)]

extern crate bincode;
extern crate pcg_rand;
extern crate png;
extern crate rand;
#[macro_use] extern crate serde_derive;
extern crate time;

mod camera;
mod collisiontype;
mod datafordrawing;
mod hitpoint;
mod lightsphere;
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
mod rendereroutputrow;
mod rendererscene;
mod renderershape;
mod renderersphere;
mod renderertriangle;
mod rhf;

use pathtracer::Pathtracer;

const NUMBER_OF_BINS: usize = 1;

fn main() {
	let mut pathtracer = Pathtracer::new(1000, 1000, 0, 8, false, vec![(0.2, 0)], 15, 1_000, 1_000_000, 15.0, 1.0, false);
	let only_post_process = false;
	if only_post_process {
		pathtracer.only_post_process();
	} else {
		pathtracer.render_images_to_png();
	}
}
