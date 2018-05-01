use std::io::{BufReader, BufWriter, Error, Read, Write};
use std::f64::consts::PI;
use std::result::Result;
use std::fs::File;

use png::{BitDepth, ColorType, Encoder, HasParameters};
use scene::SceneForPhysics;
use serde_json::{from_str, to_string};

use renderer::RendererOutput;

#[derive(Serialize, Deserialize, Debug)]
struct DataForStoring {
	width: u32,
	height: u32,
	image_scale_factor: u32,
	pixels: Vec<[f64; 3]>,
}

impl DataForStoring {
	fn new(width: u32, height: u32, image_scale_factor: u32, renderer_output_aggregate: &mut RendererOutput) -> Self {
		DataForStoring {
			width,
			height,
			image_scale_factor,
			pixels: renderer_output_aggregate.colors.clone(),
		}
	}
}

#[derive(Serialize, Deserialize, Debug)]
struct DataForDrawing {
	width: u32,
	height: u32,
	image_scale_factor: u32,
	pixels: Vec<u8>
}

impl DataForDrawing {
	fn new(data_for_storing: DataForStoring) -> Self {
		let mut pixels: Vec<u8> = Vec::new();
		for pixel in data_for_storing.pixels {
			let mut max_intensity = pixel[0];
			if pixel[1] > max_intensity {
				max_intensity = pixel[1];
			}
			if pixel[2] > max_intensity {
				max_intensity = pixel[2];
			}
			let mut factor = 1.0;
			if max_intensity > 1e-9 {
				let modified_max_intensity = (max_intensity/20.0).atan()*255.0/(PI/2.0);
				factor = modified_max_intensity/max_intensity;
			}
			pixels.push((pixel[0]*factor) as u8);
			pixels.push((pixel[1]*factor) as u8);
			pixels.push((pixel[2]*factor) as u8);
		}
		DataForDrawing {
			width: data_for_storing.width,
			height: data_for_storing.height,
			image_scale_factor: data_for_storing.image_scale_factor,
			pixels,
		}
	}
}

pub fn make_file(width: u32, height: u32, image_scale_factor: u32, renderer_output_aggregate: &mut RendererOutput, text_filename: &String, image_filename: &String) {
	let data_for_storing = DataForStoring::new(width, height, image_scale_factor, renderer_output_aggregate);
	let write_result = write_frame(&text_filename, &data_for_storing);
	match write_result {
		Err(e) => {
			println!("{:?}", e);
		}
		_ => {
		}
	}
	let read_result = read_frame(&text_filename);
	match read_result {
		Err(e) => {
			println!("{:?}", e);
		}
		Ok(data_for_storing) => {
			let data_for_drawing = DataForDrawing::new(data_for_storing);
			let draw_result = draw_frame(&image_filename, data_for_drawing);
			match draw_result {
				Err(e) => {
					println!("{:?}", e);
				}
				_ => {
				}
			}
		}
	}
}

fn write_frame(text_filename: &String, data_for_storing: &DataForStoring) -> Result<(), Error> {
	let serialized = to_string(data_for_storing)?;
	let mut file = File::create(text_filename)?;
	file.write_all(serialized.as_bytes())?;
	Ok(())
}

fn read_frame(text_filename: &String) -> Result<DataForStoring, Error> {
	let file = File::open(text_filename)?;
	let mut serialized = String::new();
	let mut buf_reader = BufReader::new(file);
	buf_reader.read_to_string(&mut serialized)?;
	let deserialized: DataForStoring = from_str(&serialized)?;
	Ok(deserialized)
}

fn draw_frame(image_filename: &String, data_for_drawing: DataForDrawing) -> Result<(), Error> {
	let file = File::create(image_filename)?;
	let ref mut buf_writer = BufWriter::new(file);
	let mut encoder = Encoder::new(buf_writer, data_for_drawing.width*data_for_drawing.image_scale_factor, data_for_drawing.height*data_for_drawing.image_scale_factor);
	encoder.set(ColorType::RGB).set(BitDepth::Eight);
	let mut writer = encoder.write_header()?;
	writer.write_image_data(&data_for_drawing.pixels)?;
	Ok(())
}

pub fn write_scene(text_filename: &String, scene: &SceneForPhysics) -> Result<(), Error> {
	let serialized = to_string(scene)?;
	let mut file = File::create(text_filename)?;
	file.write_all(serialized.as_bytes())?;
	Ok(())
}

pub fn read_scene(text_filename: &String) -> Result<SceneForPhysics, Error> {
	let file = File::open(text_filename)?;
	let mut serialized = String::new();
	let mut buf_reader = BufReader::new(file);
	buf_reader.read_to_string(&mut serialized)?;
	let deserialized: SceneForPhysics = from_str(&serialized)?;
	Ok(deserialized)
}
