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
mod rendereroutputpixel;
mod rendererscene;
mod renderershape;
mod renderersphere;
mod renderertriangle;
mod rhf;

use std::env::current_dir;

use pathtracer::Pathtracer;

// Set it to 1 if we're not doing post-processing.
const NUMBER_OF_BINS: usize = 1;
const GAMMA: f64 = 100.0;

fn main() {
	let post_process_data = false;
	match current_dir() {
		Ok(directory) => {
			let mut pathtracer = Pathtracer::new(1000, 1000, 0, 8, false, vec![(0.2, 0)], 15, 3000, 500_000, 20.0, 1000.0, false, directory);
			if post_process_data {
				pathtracer.post_process_data();
			} else {
				pathtracer.render_images_to_png();
			}
		}
		Err(e) => {
			println!("{:?}", e);
		}
	}
}
