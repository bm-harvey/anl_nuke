use faust::event::Event;
use nukers::anl_module::Anl;
use nukers::anl_module::AnlModule;
use roost::hist::RtH1D;
use cxx::SharedPtr;

fn main() {
    Anl::<Event>::new()
        .with_input_directory("K:\\tamu_data\\exp\\c12_si_35\\rkyv")
        .with_output_directory("K:\\tamu_data\\exp\\c12_si_35\\anl")
        .with_real_module(Angle::new())
        .run();
}

struct Angle {
    theta: SharedPtr<RtH1D>,
}

impl Angle {
    fn new() -> Angle {
        Angle {
            theta: roost::hist::new_h1d(String::from("theta"), String::from("theta"), 1_024, 0.0, 30.)
        }
    }
}
impl AnlModule<Event> for Angle {
    fn name(&self) -> String {
        "angle".to_string()
    }

    fn analyze_event(&mut self, _event: &Event, _idx: usize) {
        if _event.mult() != 1 {
            return;
        }
        if _event.pid_mult(6, 12) != 1 {
            return;
        }

        let theta = _event.particles()[0].momentum_MeV_per_c().theta_deg();
        self.theta.fill(theta);
    }

    fn generate_output(&mut self, output_directory: &std::path::Path) {
        let file_name: String = output_directory
            .join(format!("{}.root", self.name()))
            .to_str()
            .unwrap()
            .into();

        println!("Creating output file: {}", file_name);
        let file = roost::file::create(file_name);
        self.theta.write();
    }
}
