use cxx::SharedPtr;
use nukers::anl_module::Anl;
use nukers::anl_module::AnlModule;
use roost::hist;
use roost::hist::RtH1D;

struct CountAlphas {
    alpha_mult: SharedPtr<RtH1D>,
}

impl CountAlphas {
    pub fn new() -> CountAlphas {
        let alpha_mult =
            roost::hist::new_h1d("alpha_mult".into(), "alpha_mult".into(), 13, -0.5, 12.5);
        CountAlphas { alpha_mult }
    }
}

impl AnlModule<faust::event::Event> for CountAlphas {
    fn name(&self) -> String {
        "count_alphas_with_z10".into()
    }
    fn analyze_event(&mut self, event: &faust::event::Event, _idx: usize) {
        if event.z_mult(10) == 1 {
            self.alpha_mult.fill(event.pid_mult(2,4) as f64);
        }
    }

    fn generate_output(&mut self, _output_directory: &std::path::Path) {
        let file_name = _output_directory.join(format!("{}.root", self.name()));
        let file = roost::file::create(file_name.to_str().unwrap().to_string());
        self.alpha_mult.write();
        file.close();
    }
}

fn main() {
    let mut count_alphas = CountAlphas::new();

    let anl = Anl::<faust::event::Event>::new()
        .with_input_directory("K:\\tamu_data\\exp\\si28_c_35\\rkyv")
        .with_output_directory("K:\\tamu_data\\exp\\si28_c_35\\anl")
        .with_real_module(count_alphas) 
        .run();
}
