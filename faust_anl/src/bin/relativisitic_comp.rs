
use faust::event::Event;
use faust_anl::general_particle_selection::{ParticleSelectionRule, GeneralParticleFilter};
use faust_anl::relative_energy::AllCombinationsIter;
use faust_anl::source::RelativisticSource;
use nukers::anl_module::{Anl, AnlModule};

use cxx::SharedPtr;
use roost::hist::RtH1D;

fn main() {
    let filter =
        GeneralParticleFilter::new(faust_anl::MatchingPattern::Minimum).with_particles(2, 2, 4);

    Anl::new()
        //.with_input_directory("K:\\tamu_data\\exp\\si28_c_35\\rkyv")
        //.with_output_directory("K:\\tamu_data\\exp\\si28_c_35\\anl") 
        .with_input_directory("K:\\tamu_data\\exp\\c12_si_35\\rkyv")
        .with_output_directory("K:\\tamu_data\\exp\\c12_si_35\\anl")
        .with_filter(filter)
        .with_real_module(RelativisticTest::new())
        //.with_max_real_events(10_000_000)
        .run();
}

struct RelativisticTest {
    rules: [ParticleSelectionRule; 1],
    e_ex: SharedPtr<RtH1D>,
    e_ex_relativistic: SharedPtr<RtH1D>,
}

impl RelativisticTest {
    fn new() -> Self {
        let rules = [
            ParticleSelectionRule::new(2, 2, 4),
        ];

        Self {
            rules,
            e_ex: roost::hist::new_h1d(
                String::from("e_ex"),
                String::from(";E_{x} [MeV];Yield"),
                4_096,
                -2.,
                200.,
            ),
            e_ex_relativistic: roost::hist::new_h1d(
                String::from("e_ex_relativistic"),
                String::from(";E_{x} [MeV];Yield"),
                4_096,
                -2.,
                200.,
            ),
        }
    }
}

impl AnlModule<Event> for RelativisticTest {
    fn name(&self) -> String {
        "rel_test".to_string()
    }

    fn analyze_event(&mut self, event: &Event, _idx: usize) {

        let combos_iter = AllCombinationsIter::new(&self.rules, event);

        for source in combos_iter {

            let rel_source = RelativisticSource::new(&source);
            
            let e_ex = source.excitation_energy_MeV();
            let e_ex_relativistic = rel_source.excitation_energy_MeV();

            self.e_ex.fill(e_ex);
            self.e_ex_relativistic.fill(e_ex_relativistic);
        }

    }


    fn generate_output(&mut self, output_directory: &std::path::Path) {
        let file_name: String = output_directory
            .join(format!("{}.root", self.name()))
            .to_str()
            .unwrap()
            .into();
        println!("Creating output file: {}", file_name);
        let file = roost::file::create(file_name);
        self.e_ex.write();
        self.e_ex_relativistic.write();
        file.close();
    }
}
