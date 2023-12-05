use std::path::Path;

fn main() {
    let path_to_file = Path::new("param\\detnum_corner_z_x_y_r.txt");

    faust_filter::faust_filter::FaustFilter::new(&path_to_file);
}
