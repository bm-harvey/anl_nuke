// Prevents additional console window on Windows in release, DO NOT REMOVE!!
#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

mod analysis;
use nukers::anl_module;
use nukers::data_set;
use std::path::Path;

#[tauri::command]
fn random_number() -> Vec<[f64; 2]> {
    let size = 10_000;
    let mut result = Vec::with_capacity(size);
    for _idx in 0..size {
        result.push([
            rand::random::<f64>() * 20. - 10.,
            rand::random::<f64>() * 20. - 10.,
        ]);
    }
    result
}

fn main() {
    let files = anl_module::Anl::map_files(Path::new("D:\\tamu_data\\exp\\si28_c_35"));
    let data = data_set::DataCollection::new(&files);
    let anl = analysis::Anl::new(files);

    tauri::Builder::default()
        .invoke_handler(tauri::generate_handler![random_number])
        .run(tauri::generate_context!())
        .expect("error while running tauri application");
}
