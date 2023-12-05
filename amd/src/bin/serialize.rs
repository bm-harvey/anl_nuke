use amd::{amd_reader::EventReader, event::Event};

use std::fs::File;
use std::io::BufWriter;

use std::io::Write;

use nukers::data_set::DataSet;

fn main() {
    let timer_start = std::time::Instant::now();

    let num_time_steps = 31;
    //let path = "../amd_zn/zn70_zn70_35_g/";
    let path = "d:tamu_data/amd/zn70_zn70_35_g/";
    let mut event_reader = EventReader::new(path, num_time_steps, 140).unwrap();

    //let out_file = File::create("../amd_zn/zn70_zn70_35_g.rmp").unwrap();

    let mut idx = 0;
let mut file_idx = 0;

    let mut data_set = DataSet::<Event>::new();
    loop {
        let event: Option<Event> = event_reader.next_event();
        if event.is_none() {
            break;
        }

        idx += 1;
        if idx % 100 == 0 {
            println!("Reading event {}", idx);
        }

        let event: Event = event.unwrap();
        data_set.add_event(event);

        if idx % 1000 == 0 {
            file_idx += 1;
            let file_name = format!("../data/amd/zn70_zn70_35_g_rkyv/zn70_zn70_35_g_{}.rkyv", file_idx);
            println!("Generating file : {}", file_name);
            let out_file =
                File::create(file_name).unwrap();
            let mut out_buf = BufWriter::new(out_file);
            out_buf
                .write_all(rkyv::to_bytes::<_, 256>(&data_set).unwrap().as_slice())
                .unwrap();
            data_set.clear();
        }
        
    }

    if !data_set.is_empty() { 
    file_idx += 1;
    let out_file =
        File::create(format!("../data/amd/zn70_zn70_35_g_rkyv/zn70_zn70_35_g_{}.rkyv", file_idx)).unwrap();
    let mut out_buf = BufWriter::new(out_file);
    out_buf
        .write_all(rkyv::to_bytes::<_, 256>(&data_set).unwrap().as_slice())
        .unwrap();
    data_set.clear();
    }

    println!(
        "Time to serialize {} events  : {} s ",
        idx,
        timer_start.elapsed().as_secs_f32()
    );
}
