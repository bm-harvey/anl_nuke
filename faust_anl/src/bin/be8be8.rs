use faust::event::{Event, Particle};
use faust_anl::relative_energy::AllCombinationsIter;
use faust_anl::source::Source;
use faust_anl::GeneralParticleMixer;
use faust_anl::{general_particle_selection::ParticleSelectionRule, GeneralParticleFilter};
use itertools::Itertools;
use nukers::anl_module::{Anl, AnlModule};
use rand::Rng;

use cxx::SharedPtr;
use roost::hist::RtH1D;

fn main() {
    let filter =
        GeneralParticleFilter::new(faust_anl::MatchingPattern::Minimum).with_particles(4, 2, 4);
    let mut mixer = GeneralParticleMixer::new().with_particles(4, 2, 4);
    mixer.set_faust_filter();

    Anl::new()
        .with_input_directory("K:\\tamu_data\\exp\\o16_c_35\\rkyv")
        .with_output_directory("K:\\tamu_data\\exp\\o16_c_35\\anl")
        //.with_input_directory("K:\\tamu_data\\exp\\c12_si_35\\rkyv")
        //.with_output_directory("K:\\tamu_data\\exp\\c12_si_35\\anl")
        .with_filter(filter)
        //.with_event_generator(nukers::anl_module::EventGenerator::Scrambler(Box::new(mixer)))
        .with_event_generator(nukers::anl_module::EventGenerator::Mixer(Box::new(mixer)))
        .with_real_module(Be8Correlations::new())
        .with_mixed_module(Be8Correlations::new())
        .with_max_mixed_events(nukers::anl_module::MixedEventMaximum::Absolute(100_000_000))
        .run();
}

struct Be8Correlations {
    e_ex: SharedPtr<RtH1D>,
    theta: SharedPtr<RtH1D>,
}

impl Be8Correlations {
    fn new() -> Self {
        Self {
            e_ex: roost::hist::new_h1d(
                String::from("e_ex"),
                String::from(";E^{*} [MeV];Yield"),
                4_096 * 2,
                -2.,
                100.,
            ),
            theta: roost::hist::new_h1d(
                String::from("theta"),
                String::from(";#theta [deg];Yield"),
                180,
                0.,
                180.,
            ),
        }
    }
}

impl AnlModule<Event> for Be8Correlations {
    fn name(&self) -> String {
        "be8be8".to_string()
    }

    fn analyze_event(&mut self, event: &Event, _idx: usize) {
        let rules = [ParticleSelectionRule::new(4, 2, 4)];

        let combo_iter = AllCombinationsIter::new(&rules, &event);

        let mut rng = rand::thread_rng();
        'alpha_4: for alphas_4_source in combo_iter {
            for pair in (0..4).combinations(2) {
                let mut be8_1 = Source::new();
                be8_1.add_particle(alphas_4_source.particles()[pair[0]]);
                be8_1.add_particle(alphas_4_source.particles()[pair[1]]);

                let other_indices = (0..4)
                    .filter(|idx| !pair.contains(idx))
                    .collect::<Vec<usize>>();

                let mut be8_2 = Source::new();
                be8_2.add_particle(alphas_4_source.particles()[other_indices[0]]);
                be8_2.add_particle(alphas_4_source.particles()[other_indices[1]]);

                let be8_1_e_ex = be8_1.excitation_energy_MeV();
                let be8_2_e_ex = be8_2.excitation_energy_MeV();

                if be8_1_e_ex < 0.3 && be8_2_e_ex < 0.3 {
                    let e_ex = alphas_4_source.excitation_energy_MeV();
                    self.e_ex.fill(e_ex);

                    let v_o16_lab = alphas_4_source.particles()[0].velocity_c();
                    let mut v_be8_1_o16 = &be8_1.velocity_c() - &alphas_4_source.velocity_c();
                    if rng.gen_bool(0.5) {
                        v_be8_1_o16.scale(-1.);
                    }
                    let theta = v_o16_lab.inner_angle_deg(&v_be8_1_o16);
                    self.theta.fill(theta);

                    break 'alpha_4;
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
        self.theta.write();
        file.close();
    }
}

enum Be8 {
    SS0p,
    SS2p,
    SS4p,
}

fn decay_pathway<'a>(c12: &'a Source<'a>) -> Option<(Be8, Source<'a>, &'a Particle)> {
    debug_assert_eq!(c12.particles().len(), 3);
    for alpha in c12.iter() {
        debug_assert!(alpha.is_pid(2, 4))
    }

    let mut be8_sources = c12
        .iter()
        .combinations(2)
        .map(|alphas_2| {
            let mut be8 = Source::<'a>::new();
            alphas_2.iter().for_each(|&alpha| {
                be8.add_particle(alpha);
            });
            be8
        })
        .collect::<Vec<_>>();

    be8_sources.sort_by(|a, b| {
        a.relative_kinetic_energy_MeV()
            .partial_cmp(&b.relative_kinetic_energy_MeV())
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    for (source, energy) in be8_sources.into_iter().map(|be8| {
        let e_rel = be8.relative_kinetic_energy_MeV();
        (be8, e_rel)
    }) {
        if (0.0..0.3).contains(&energy) {
            let alpha = c12
                .iter()
                .find(|&&alpha| {
                    !source
                        .iter()
                        .any(|&be8_alpha| std::ptr::eq(be8_alpha, alpha))
                })
                .unwrap();
            return Some((Be8::SS0p, source, &alpha));
        }
        if (2.0..4.0).contains(&energy) {
            let alpha = c12
                .iter()
                .find(|&&alpha| {
                    !source
                        .iter()
                        .any(|&be8_alpha| std::ptr::eq(be8_alpha, alpha))
                })
                .unwrap();
            return Some((Be8::SS2p, source, &alpha));
        }
        if (9.0..13.0).contains(&energy) {
            let alpha = c12
                .iter()
                .find(|&&alpha| {
                    !source
                        .iter()
                        .any(|&be8_alpha| std::ptr::eq(be8_alpha, alpha))
                })
                .unwrap();
            return Some((Be8::SS4p, source, &alpha));
        }
    }
    None
}
