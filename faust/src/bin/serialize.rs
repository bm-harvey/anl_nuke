use faust::event::Event;
use faust::event_reader::Reader;
use nukers::data_set::DataSet;
use std::fs::File;
use std::io::BufWriter;
use std::io::Write;

use std::path::Path;

use rayon::prelude::*;

fn main() {
    let start = std::time::Instant::now();

    //let input_files_iter = std::fs::read_dir("K:tamu_data/exp/si28_c_35/ascii/").unwrap();
    //let input_files_iter = std::fs::read_dir("K:tamu_data/exp/o16_c_35/ascii/").unwrap();
    let input_files_iter = std::fs::read_dir("/data/sjygroup/sjy20/bmharvey/acs/c12_si_35/ascii").unwrap();
    //let input_files_iter = std::fs::read_dir("/data/sjygroup/sjy20/bmharvey/acs/si28_c_35/ascii").unwrap();
    //let input_files_iter = std::fs::read_dir("K:tamu_data/exp/c12_si_35/ascii/").unwrap();

    input_files_iter.for_each(|file| {
        let file = file.unwrap();
        let file_name: String = file.path().to_str().unwrap().into();
        let mut reader = Reader::builder()
            .with_capacity(10_000_000)
            .with_input_file(&file_name)
            .build();

        let input_path = Path::new(&file_name);

        input_path.file_name().unwrap().to_str().unwrap();

        let output_file = input_path
            .file_name()
            .unwrap()
            .to_str()
            .unwrap()
            .split('.')
            .collect::<Vec<_>>()[0]
            .to_string()
            + ".rkyv";

        //let output_file_name = format!("K:\\tamu_data\\exp\\o16_c_35\\rkyv\\{}", output_file);
        //let output_file_name = format!("/data/sjygroup/sjy20/bmharvey/acs/si28_c_35/rkyv/{}", output_file);
        let output_file_name = format!("/data/sjygroup/sjy20/bmharvey/acs/c12_si_35/rkyv/{}", output_file);
        //let output_file_name = format!("K:\\tamu_data\\exp\\c12_si_35\\rkyv\\{}", output_file);
        //let output_file_name = format!("K:\\tamu_data\\exp\\si28_c_35\\rkyv\\{}", output_file);

        let mut data_set = DataSet::<Event>::new();

        let mut idx = 0;

        loop {
            let event: Option<Event> = reader.next_event();
            if event.is_none() {
                break;
            }

            let event: Event = event.unwrap();
            if event.mult() == 0 {
                continue;
            }
            data_set.add_event(event);

            idx += 1;
            if idx % 1_000_000 == 0 {
                println!("Reading event {} from {}", idx, &file_name.clone());
            }
        }

        println!("Generating file : {}", output_file_name);
        let out_file = File::create(output_file_name).unwrap();
        let mut out_buf = BufWriter::new(out_file);
        out_buf
            .write_all(rkyv::to_bytes::<_, 256>(&data_set).unwrap().as_slice())
            .unwrap();
    });

    println!(
        "Time to read and write data: {} seconds",
        start.elapsed().as_secs_f32()
    );

    /*

    let mut idx = 0;
    let mut file_idx = 0;

    let mut data_set = DataSet::<Event>::new();

    loop {
        let event: Option<Event> = reader.next_event();
        if event.is_none() {
            break;
        }

        let event: Event = event.unwrap();
        if event.mult() == 0 {
            continue;
        }
        idx += 1;
        data_set.add_event(event);

        if idx % 100_000 == 0 {
            println!("Reading event {}", idx);
        }

        if idx % 10_000_000 == 0 {
            file_idx += 1;
            let file_name = format!("K:tamu_data/exp/si28_c_35/rkyv/file_{}.rkyv", file_idx);
            //let file_name = format!("./data/si28_c_35_rkyv/file_{}.rkyv", file_idx);
            println!("Generating file : {}", file_name);
            let out_file = File::create(file_name).unwrap();
            let mut out_buf = BufWriter::new(out_file);
            out_buf
                .write_all(rkyv::to_bytes::<_, 256>(&data_set).unwrap().as_slice())
                .unwrap();
            data_set.clear();
        }
    }

    if !data_set.is_empty() {
        file_idx += 1;
        let file_name = format!("K:tamu_data/exp/si28_c_35/rkyv/file_{}.rkyv", file_idx);
        println!("Generating file : {}", file_name);
        let out_file = File::create(file_name).unwrap();
        let mut out_buf = BufWriter::new(out_file);
        out_buf
            .write_all(rkyv::to_bytes::<_, 256>(&data_set).unwrap().as_slice())
            .unwrap();
        data_set.clear();
    }
    println!(
        "Serializing {} events into {} files too {} seconds",
        idx,
        file_idx,
        start.elapsed().as_secs_f32()
    );
    */
}
