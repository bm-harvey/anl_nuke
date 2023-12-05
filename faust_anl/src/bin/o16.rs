use faust::event::Event;
use faust_anl::general_particle_selection::ParticleSelectionRule;
use faust_anl::relative_energy::AllCombinationsIter;
use faust_anl::source::Source;
use itertools::Itertools;
use nukers::anl_module::{Anl, AnlModule, EventFilter};

use cxx::SharedPtr;
use roost::hist::{RtH1D, RtH2D};
use roost::tree::{RtBranch, RtIntBranch, RtTree};

fn main() {
    Anl::new()
        .with_input_directory("K:\\tamu_data\\exp\\si28_c_35\\rkyv")
        .with_output_directory("K:\\tamu_data\\exp\\si28_c_35\\anl")
        .with_real_module(O16::new())
        //.with_max_real_events(5_000_000)
        .run();
}

struct O16 {
    e_rel_hist: SharedPtr<RtH1D>,
    e_rel_hist_3_v_4: SharedPtr<RtH2D>,
    e_rel_hist_2_v_4: SharedPtr<RtH2D>,
    e_ex_c12_he4: SharedPtr<RtH1D>,
    e_ex_c11_he4: SharedPtr<RtH1D>,
    e_ex_c13_he4: SharedPtr<RtH1D>,
    e_ex_thru_hoyle: SharedPtr<RtH1D>,
}

impl O16 {
    fn new() -> Self {
        Self {
            e_rel_hist: roost::hist::new_h1d(
                String::from("e_rel"),
                String::from(";E_{{rel}} [MeV];Yield"),
                4_096 * 2,
                -2.,
                100.,
            ),
            e_ex_thru_hoyle: roost::hist::new_h1d(
                String::from("e_ex_thru_hoyle"),
                String::from(";{}^{12}C + #alpha E_{x} [MeV];Yield"),
                4_096 * 2,
                -2.,
                100.,
            ),
            e_ex_c12_he4: roost::hist::new_h1d(
                String::from("e_ex_c12_he4"),
                String::from(";{}^{12}C + #alpha E_{x} [MeV];Yield"),
                4_096 * 2,
                -2.,
                100.,
            ),
            e_ex_c11_he4: roost::hist::new_h1d(
                String::from("e_ex_c11_he4"),
                String::from(";{}^{11}C + #alpha E_{x} [MeV];Yield"),
                4_096 * 2,
                -2.,
                100.,
            ),
            e_ex_c13_he4: roost::hist::new_h1d(
                String::from("e_ex_c13_he4"),
                String::from(";{}^{13}C + #alpha E_{x} [MeV];Yield"),
                4_096 * 2,
                -2.,
                100.,
            ),
            e_rel_hist_3_v_4: roost::hist::new_h2d(
                String::from("e_rel_3_v_4"),
                String::from(
                    ";#alpha#alpha#alpha#alpha E_{rel} [MeV];#alpha#alpha#alpha E_{rel} [MeV]",
                ),
                1_024,
                -2.,
                100.,
                1_024,
                -2.,
                100.,
            ),
            e_rel_hist_2_v_4: roost::hist::new_h2d(
                String::from("e_rel_2_v_4"),
                String::from(";#alpha#alpha#alpha#alpha E_{rel} [MeV];#alpha#alpha E_{rel} [MeV]"),
                1_024,
                -2.,
                100.,
                1_024,
                -2.,
                100.,
            ),
        }
    }
}

impl AnlModule<Event> for O16 {
    fn name(&self) -> String {
        "o16".to_string()
    }

    fn analyze_event(&mut self, event: &Event, _idx: usize) {
        if event.pid_mult(2, 4) >= 4 {
            let rules = [ParticleSelectionRule::new(4, 2, 4)];
            let combo_iter = AllCombinationsIter::new(&rules, &event);

            for combo_4_alpha in combo_iter {
                let e_rel_4_alpha = combo_4_alpha.relative_kinetic_energy_MeV();
                let e_ex = e_rel_4_alpha - combo_4_alpha.q_value_MeV();
                self.e_rel_hist.fill(e_rel_4_alpha);

                for combo_3_alpha in combo_4_alpha.iter().combinations(3) {
                    let mut source = Source::new();
                    for alpha in combo_3_alpha {
                        source.add_particle(alpha);
                    }

                    let e_rel_3_alpha = source.relative_kinetic_energy_MeV();
                    self.e_rel_hist_3_v_4.fill(e_rel_4_alpha, e_rel_3_alpha);

                    if e_rel_3_alpha < 0.6 {
                        //let e_ex = .excitation_energy_MeV();
                        self.e_ex_thru_hoyle.fill(e_ex);
                    }
                }
                for combo_2_alpha in combo_4_alpha.iter().combinations(2) {
                    let mut source = Source::new();
                    for alpha in combo_2_alpha {
                        source.add_particle(alpha);
                    }

                    let e_rel_2_alpha = source.relative_kinetic_energy_MeV();
                    self.e_rel_hist_2_v_4.fill(e_rel_4_alpha, e_rel_2_alpha);
                }
            }
        }

        /*
        let rules = [
            ParticleSelectionRule::new(1, 6, 12),
            ParticleSelectionRule::new(1, 2, 4),
        ];
        let combo_iter = AllCombinationsIter::new(&rules, &event);
        for combo in combo_iter {
            let e_ex = combo.excitation_energy_MeV();
            self.e_ex_c12_he4.fill(e_ex);
        }

        //
        let rules = [
            ParticleSelectionRule::new(1, 6, 11),
            ParticleSelectionRule::new(1, 2, 4),
        ];
        let combo_iter = AllCombinationsIter::new(&rules, &event);
        for combo in combo_iter {
            let e_ex = combo.excitation_energy_MeV();
            self.e_ex_c11_he4.fill(e_ex);
        }

        //
        let rules = [
            ParticleSelectionRule::new(1, 6, 13),
            ParticleSelectionRule::new(1, 2, 4),
        ];
        let combo_iter = AllCombinationsIter::new(&rules, &event);
        for combo in combo_iter {
            let e_ex = combo.excitation_energy_MeV();
            self.e_ex_c13_he4.fill(e_ex);
        }
        */
    }

    fn generate_output(&mut self, output_directory: &std::path::Path) {
        let file_name: String = output_directory
            .join(format!("{}.root", self.name()))
            .to_str()
            .unwrap()
            .into();
        println!("Creating output file: {}", file_name);
        let file = roost::file::create(file_name);
        self.e_rel_hist.write();
        self.e_ex_thru_hoyle.write();
        //self.e_ex_c11_he4.write();
        //self.e_ex_c12_he4.write();
        //self.e_ex_c13_he4.write();
        self.e_rel_hist_3_v_4.write();
        self.e_rel_hist_2_v_4.write();
        file.close();
    }
}
