use std::io::{BufReader, BufWriter, Error, Read, Write};
use std::result::Result;
use std::fs::File;

use png::{BitDepth, ColorType, Encoder, HasParameters};
use serde_json::{from_str, to_string};

use datafordrawing::DataForDrawing;
use dataforstoring::DataForStoring;
use rendereroutput::RendererOutput;
use physicsscene::PhysicsScene;

pub fn make_file(width: u32, height: u32, renderer_output_aggregate: &mut RendererOutput, text_filename: &str, image_filename: &str, number_of_rays_total: u64) {
	let data_for_storing = DataForStoring::new(width, height, number_of_rays_total, renderer_output_aggregate);
	let write_result = write_frame(&text_filename, &data_for_storing);
	if let Err(e) = write_result {
		println!("{:?}", e);
	}
	let read_result = read_frame(&text_filename);
	match read_result {
		Err(e) => {
			println!("{:?}", e);
		}
		Ok(data_for_storing) => {
			let data_for_drawing = DataForDrawing::new(data_for_storing);
			let draw_result = draw_frame(&image_filename, &data_for_drawing);
			if let Err(e) = draw_result {
				println!("{:?}", e);
			}
		}
	}
}

fn write_frame(text_filename: &str, data_for_storing: &DataForStoring) -> Result<(), Error> {
	let serialized = to_string(data_for_storing)?;
	let mut file = File::create(text_filename)?;
	file.write_all(serialized.as_bytes())?;
	Ok(())
}

fn read_frame(text_filename: &str) -> Result<DataForStoring, Error> {
	let file = File::open(text_filename)?;
	let mut serialized = String::new();
	let mut buf_reader = BufReader::new(file);
	buf_reader.read_to_string(&mut serialized)?;
	let deserialized: DataForStoring = from_str(&serialized)?;
	Ok(deserialized)
}

fn draw_frame(image_filename: &str, data_for_drawing: &DataForDrawing) -> Result<(), Error> {
	let file = File::create(image_filename)?;
	let buf_writer = &mut BufWriter::new(file);
	let mut encoder = Encoder::new(buf_writer, data_for_drawing.width, data_for_drawing.height);
	encoder.set(ColorType::RGB).set(BitDepth::Eight);
	let mut writer = encoder.write_header()?;
	writer.write_image_data(&data_for_drawing.pixels)?;
	Ok(())
}

pub fn write_scene(text_filename: &str, scene: &PhysicsScene) -> Result<(), Error> {
	let serialized = to_string(scene)?;
	let mut file = File::create(text_filename)?;
	file.write_all(serialized.as_bytes())?;
	Ok(())
}

pub fn read_scene(text_filename: &str) -> Result<PhysicsScene, Error> {
	let file = File::open(text_filename)?;
	let mut serialized = String::new();
	let mut buf_reader = BufReader::new(file);
	buf_reader.read_to_string(&mut serialized)?;
	let deserialized: PhysicsScene = from_str(&serialized)?;
	Ok(deserialized)
}
