use faust::event::Event;
use faust_anl::general_particle_selection::ParticleSelectionRule;
use faust_anl::relative_energy::AllCombinationsIter;
use faust_anl::source::Source;
use nukers::anl_module::{Anl, AnlModule, EventFilter};

use cxx::SharedPtr;
use faust::nuclear_masses::NUCLEAR_DB;
use itertools::Itertools;
use roost::hist::{RtH1D, RtH2D};
use roost::tree::{RtBranch, RtIntBranch, RtTree};

fn main() {
    Anl::new()
        .with_input_directory("K:\\tamu_data\\exp\\si28_c_35\\rkyv")
        .with_output_directory("K:\\tamu_data\\exp\\si28_c_35\\anl")
        .with_real_module(AngularCorrelations::new())
        //.with_max_real_events(1_000_000)
        .run();
}

struct AngularCorrelations {
    rules: [ParticleSelectionRule; 1],
    e_ex: SharedPtr<RtH1D>,
    e_ex_thru_hoyle: SharedPtr<RtH1D>,
    open_angle_thru_hoyle: SharedPtr<RtH1D>,
}

impl AngularCorrelations {
    fn new() -> Self {
        let rules = [ParticleSelectionRule::new(4, 2, 4)];

        Self {
            rules,
            e_ex: roost::hist::new_h1d(
                String::from("e_rel"),
                String::from(";E^{*} [MeV];Yield"),
                4_096,
                -2.,
                120.,
            ),
            e_ex_thru_hoyle: roost::hist::new_h1d(
                String::from("e_rel_thru_hoyle"),
                String::from("{}^{12}C(0^+);E^{*} [MeV];Yield"),
                4_096,
                -2.,
                120.,
            ),
            open_angle_thru_hoyle: roost::hist::new_h1d(
                String::from("open_angle_thru_hoyle"),
                String::from("{}^{12}C(0^+);#theta [deg];Yield"),
                1_024,
                0.,
                180.,
            ),
        }
    }

    fn has_be8_gs(source: &Source) -> bool {
        source.iter().combinations(2).any(|combo| {
            let mut maybe_be8_gs = Source::new();
            for particle in combo {
                maybe_be8_gs.add_particle(particle);
            }
            maybe_be8_gs.excitation_energy_MeV() < 0.2
        })
    }
    fn has_c12_3m(source: &Source) -> bool {
        source.iter().combinations(3).any(|combo| {
            let mut maybe_c12_3m = Source::new();
            for particle in combo {
                maybe_c12_3m.add_particle(particle);
            }
            let e_star = maybe_c12_3m.excitation_energy_MeV();

            9. < e_star && e_star < 11.
        })
    }
    fn has_hoyle(source: &Source) -> bool {
        source.iter().combinations(3).any(|combo| {
            let mut maybe_hoyle = Source::new();
            for particle in combo {
                maybe_hoyle.add_particle(particle);
            }
            let e_star = maybe_hoyle.excitation_energy_MeV();

            e_star < 8.
        })
    }
}

impl AnlModule<Event> for AngularCorrelations {
    fn name(&self) -> String {
        "angular_correlations".to_string()
    }
    fn analyze_event(&mut self, event: &Event, _idx: usize) {
        let o16_thru_4a_sources_iter = AllCombinationsIter::new(&self.rules, event);

        for source in o16_thru_4a_sources_iter {
            let e_star = source.excitation_energy_MeV();
            self.e_ex.fill(e_star);

            for combo in (0..source.num_particles()).combinations(3) {
                let mut source_3a = Source::new();
                for idx in combo.iter() {
                    source_3a.add_particle(source.particles()[*idx]);
                }

                let e_star = source_3a.excitation_energy_MeV();

                let is_hoyle = e_star < 8.;
                if is_hoyle {
                    self.e_ex_thru_hoyle.fill(source.excitation_energy_MeV());

                    let mut alpha_idx: usize = 0;
                    for idx in 0..3 {
                        if !combo.contains(&idx) {
                            alpha_idx = idx;
                            break;
                        }
                    }

                    let alpha = source_3a.particles()[alpha_idx];

                    let p_alpha = alpha.momentum_MeV_per_c();
                    let p_o16 = source.momentum_MeV_per_c();
                    let theta = p_o16.inner_angle_deg(&p_alpha);

                    self.open_angle_thru_hoyle.fill(theta);
                }
            }
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
        self.e_ex_thru_hoyle.write();
        self.open_angle_thru_hoyle.write();
        file.close();
    }
}
