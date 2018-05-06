// To install the nightly channel of rust:
// rustup install nightly
// To set the nightly channel of rust to default:
// rustup default nightly
// To create a video file from the output images:
// ffmpeg -r 60 -f image2 -s 1000x1000 -start_number [start_frame_number] -i %04d.png -vframes [number_of_frames] -vcodec libx264 -crf 10 -pix_fmt yuv420p test.mp4
// Eg ffmpeg -r 60 -f image2 -s 1000x1000 -start_number 0 -i %04d.png -vframes 60 -vcodec libx264 -crf 10 -pix_fmt yuv420p test.mp4

// Timings using NUMBER_OF_BINS = 64, width = 1000, height = 1000, image_scale_factor = 1, max_distance = 4.0, patch_radius = 1, search_window_radius = 5, number_of_scales = 1, number_of_rays (per thread) = ???.
// Note that the renderer does constant work per thread regardless of the number of threads - N threads means N times more work in total.
// 	#threads 	render		post-process	render speedup		post-process speedup
//	1 			55s			150				1					1
//	2 			60s			82s				1.8					1.8
//	3 			62s			49s				2.7					3.1
// 	4  			62s			49s				3.5					3.1
// 	5  			65s			51s				4.2					2.9
// 	6  			65s			49s				5.1					3.1
// 	7  			70s			52s				5.5					2.9
//	8			70s			55s				6.3					2.7

#![feature(placement_in_syntax)]
#![feature(collection_placement)]
#[macro_use]

extern crate serde_derive;
extern crate serde_json;
extern crate png;
extern crate pcg_rand;
extern crate rand;
extern crate time;

mod collisiontype;
mod datafordrawing;
mod dataforstoring;
mod math;
mod pathtracer;
mod pixel;
mod pngfile;
mod ray;
mod renderer;
mod rendereroutput;
mod rhf;
mod sceneforrendering;
mod sceneforphysics;
mod spherefast;
mod spherelightsource;
mod spherephysics;
mod trianglefast;
mod trianglelightsource;
mod trianglephysics;

use pathtracer::Pathtracer;

const NUMBER_OF_BINS: usize = 2;

fn main() {
	let mut pathtracer = Pathtracer::new(1000, 1000, 1, 0, 8, 3, false, vec![(3.0, 1)], 10, 1, 100_000_000, false, false);
	pathtracer.render_images_to_png();
}
