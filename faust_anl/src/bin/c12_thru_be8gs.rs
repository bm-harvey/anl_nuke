use std::f64::consts::PI;

use faust::event::{Event, Particle};
use faust::phys_vec::PhysVec;
use faust_anl::relative_energy::AllCombinationsIter;
use faust_anl::source::{RelativisticSource, Source};
use faust_anl::GeneralParticleMixer;
use faust_anl::{general_particle_selection::ParticleSelectionRule, GeneralParticleFilter};
use itertools::Itertools;
use nukers::anl_module::{Anl, AnlModule};
use rand::seq::SliceRandom;

use cxx::SharedPtr;
use faust::nuclear_masses::NUCLEAR_DB;
use roost::hist::{RtH1D, RtH2D};
//use roost::tree::{RtBranch, RtIntBranch, RtTree};

fn main() {
    let filter =
        GeneralParticleFilter::new(faust_anl::MatchingPattern::Standard).with_particles(3, 2, 4);

    let mut mixer = GeneralParticleMixer::new().with_particles(3, 2, 4);
    mixer.set_faust_filter();

    Anl::new()
        .with_input_directory("/data/sjygroup/sjy20/bmharvey/acs/c12_si_35/rkyv")
        .with_output_directory("/data/sjygroup/sjy20/bmharvey/acs/c12_si_35/anl")
        .with_filter(filter)
        //.with_event_generator(nukers::anl_module::EventGenerator::Scrambler(Box::new(mixer)))
        .with_event_generator(nukers::anl_module::EventGenerator::Mixer(Box::new(mixer)))
        .with_real_module(AaaThruBe8gs::new())
        //.with_mixed_module(AaaThruBe8gs::new())
        //.with_max_mixed_events(nukers::anl_module::MixedEventMaximum::Absolute(5_000_000))
        //.with_max_mixed_events(nukers::anl_module::MixedEventMaximum::Absolute(500_000_000))
        .run();
}

struct AaaThruBe8gs {
    e_ex: SharedPtr<RtH1D>,
    e_ex_2a: SharedPtr<RtH1D>,
    e_ex_2a_be8gs_fills_once: SharedPtr<RtH1D>,
    e_ex_2a_from_c12_3m: SharedPtr<RtH1D>,
    e_ex_2a_not_from_c12_3m: SharedPtr<RtH1D>,
    e_ex_2a_thru_be8_any: SharedPtr<RtH1D>,
    e_ex_small_2a: SharedPtr<RtH1D>,
    e_ex_thru_be8_gs: SharedPtr<RtH1D>,
    e_ex_thru_be_2p: SharedPtr<RtH1D>,
    e_ex_2_v_3: SharedPtr<RtH2D>,

    dalitz_all: SharedPtr<RtH2D>,
    dalitz_thru_c12_0p: SharedPtr<RtH2D>,
    dalitz_thru_c12_3m: SharedPtr<RtH2D>,
    dalitz_thru_c12_low: SharedPtr<RtH2D>,

    e_missing_v_e_ex: SharedPtr<RtH2D>,

    phi_theta_1_3m: SharedPtr<RtH2D>,
    phi_theta_2_3m: SharedPtr<RtH2D>,
    phi_theta_0p: SharedPtr<RtH2D>,
    phi_theta_high_ex: SharedPtr<RtH2D>,

    e_ex_3_vs_theta: SharedPtr<RtH2D>,
    e_ex_3_vs_theta_thru_4p: SharedPtr<RtH2D>,
    ex_3_vs_theta_thru_2p: SharedPtr<RtH2D>,
    e_ex_3_vs_theta_thru_0p: SharedPtr<RtH2D>,
    e_ex_2_vs_theta: SharedPtr<RtH2D>,

    missing_energy_gate: f64,
}

//#[rustfmt::skip]
impl AaaThruBe8gs {
    fn new() -> Self {
        Self {
            missing_energy_gate: 15.,

            #[rustfmt::skip]
            e_ex: roost::hist::new_h1d( String::from("e_ex"), String::from(";E^{*} [MeV];Yield"), 4_096 * 2, -2., 100.,),
            #[rustfmt::skip]
            dalitz_all: roost::hist::new_h2d( String::from("dalitz_all"), String::from(";#sqrt{3}(E_{2} - E_{1});2E_{3} - E_{1} - E_{2}[MeV]"), 1_024, -1.1, 1.1, 1_024, -1.1, 1.1,),
            #[rustfmt::skip]
            dalitz_thru_c12_0p: roost::hist::new_h2d( String::from("dalitz_thru_c12_0p"), String::from(";#sqrt{3}(E_{2} - E_{1});2E_{3} - E_{1} - E_{2}[MeV]"), 1_024, -1.1, 1.1, 1_024, -1.1, 1.1,),
            #[rustfmt::skip]
            dalitz_thru_c12_low: roost::hist::new_h2d( String::from("dalitz_thru_c12_low"), String::from(";#sqrt{3}(E_{2} - E_{1});2E_{3} - E_{1} - E_{2}[MeV]"), 1_024, -1.1, 1.1, 1_024, -1.1, 1.1,),
            #[rustfmt::skip]
            dalitz_thru_c12_3m: roost::hist::new_h2d( String::from("dalitz_thru_c12_3m"), String::from(";#sqrt{3}(E_{2} - E_{1});2E_{3} - E_{1} - E_{2}[MeV]"), 1_024, -1.1, 1.1, 1_024, -1.1, 1.1,),
            #[rustfmt::skip]
            ex_3_vs_theta_thru_2p: roost::hist::new_h2d( String::from("ex_3_vs_theta_through_2p"), String::from(";#theta [deg];E^{*} [MeV]"), 180, 0., 180., 512, -2., 100.,),
            #[rustfmt::skip]
            e_missing_v_e_ex: roost::hist::new_h2d( String::from("e_missing_v_e_ex"), String::from(";E^{*} [MeV];E_{missing} [MeV]"), 256, -2., 100., 256, -2., 500.,),
            #[rustfmt::skip]
            phi_theta_1_3m: roost::hist::new_h2d( String::from("phi_theta_1_3m"), String::from(";#theta [deg];#phi [deg]"), 90, 0., 180., 90, 0., 360.,),
            #[rustfmt::skip]
            phi_theta_2_3m: roost::hist::new_h2d( String::from("phi_theta_2_3m"), String::from(";#theta [deg];#phi [deg]"), 90, 0., 180., 90, 0., 360.,),
            #[rustfmt::skip]
            phi_theta_high_ex: roost::hist::new_h2d( String::from("phi_theta_high_ex"), String::from(";#theta [deg];#phi [deg]"), 90, 0., 180., 90, 0., 360.,),
            #[rustfmt::skip]
            phi_theta_0p: roost::hist::new_h2d( String::from("phi_theta_0p"), String::from(";#theta [deg];#phi [deg]"), 90, 0., 180., 90, 0., 360.,),
            #[rustfmt::skip]
            e_ex_3_vs_theta_thru_0p: roost::hist::new_h2d( String::from("ex_3_vs_theta_through_0p"), String::from(";#theta [deg];E^{*} [MeV]"), 180, 0., 180., 512, -2., 100.,),
            #[rustfmt::skip]
            e_ex_3_vs_theta_thru_4p: roost::hist::new_h2d( String::from("ex_3_vs_theta_through_4p"), String::from(";#theta [deg];E^{*} [MeV]"), 180, 0., 180., 512, -2., 100.,),
            #[rustfmt::skip]
            e_ex_3_vs_theta: roost::hist::new_h2d( String::from("ex_3_vs_theta"), String::from(";#theta [deg];E^{*} [MeV]"), 180, 0., 180., 2_048, -2., 100.,),
            #[rustfmt::skip]
            e_ex_2_vs_theta: roost::hist::new_h2d( String::from("ex_2_vs_theta"), String::from(";#theta [deg];E^{*} [MeV]"), 180, 0., 180., 2_048, -2., 100.,),
            #[rustfmt::skip]
            e_ex_2a: roost::hist::new_h1d( String::from("e_ex_2a"), String::from(";E^{*} [MeV];Yield"), 4_096 * 2, -2., 100.,),
            #[rustfmt::skip]
            e_ex_2a_be8gs_fills_once: roost::hist::new_h1d( String::from("e_ex_2a_be8gs_fills_once"), String::from(";E^{*} [MeV];Yield"), 4_096 * 2, -2., 100.,),
            #[rustfmt::skip]
            e_ex_2a_from_c12_3m: roost::hist::new_h1d( String::from("e_ex_2a_from_c12_3m"), String::from(";E^{*} [MeV];Yield"), 4_096 * 2, -2., 100.,),
            #[rustfmt::skip]
            e_ex_2a_not_from_c12_3m: roost::hist::new_h1d( String::from("e_ex_2a_not_from_c12_3m"), String::from(";E^{*} [MeV];Yield"), 4_096 * 2, -2., 100.,),
            #[rustfmt::skip]
            e_ex_2a_thru_be8_any: roost::hist::new_h1d( String::from("e_ex_2a_be8_any"), String::from(";E^{*} [MeV];Yield"), 4_096 * 2, -2., 100.,),
            #[rustfmt::skip]
            e_ex_small_2a: roost::hist::new_h1d( String::from("e_ex_small_2a"), String::from(";E^{*} [MeV];Yield"), 4_096 * 2, -2., 100.,),
            #[rustfmt::skip]
            e_ex_thru_be8_gs: roost::hist::new_h1d( String::from("e_ex_thru_be8gs"), String::from("Through {}^{8}Be (g.s.);E^{*} [MeV];Yield"), 4_096 * 2, -2., 100.,),
            #[rustfmt::skip]
            e_ex_thru_be_2p: roost::hist::new_h1d( String::from("e_ex_thru_be_2p"), String::from("Through {}^{8}Be (2+);E^{*} [MeV];Yield"), 4_096 * 2, -2., 100.,),
            #[rustfmt::skip]
            e_ex_2_v_3: roost::hist::new_h2d( String::from("e_ex_2_v_3"), String::from(";#alpha#alpha#alpha E^{*} [MeV];#alpha#alpha E^{*} [MeV]"), 1_024, -2., 100., 1_024, -2., 100.,),
        }
    }

    #[allow(dead_code)]
    fn angles_method_1(c12: &Source, alpha: &Particle) -> (f64, f64) {
        let beam = PhysVec::from_cartesian(0., 0., 1.);

        let vel_c12_lab = c12.velocity_c();
        let vel_alpha_c_frame = &alpha.classical_velocity_c() - &vel_c12_lab;
        let ortho = beam.cross(&vel_c12_lab);
        let theta = { ortho.inner_angle_deg(&vel_alpha_c_frame) };

        let phi = {
            let projected_vel_alpha =
                ortho.as_normalized_to(vel_alpha_c_frame.mag() * theta.to_radians().cos());

            let transverse_alpha_vel = &vel_alpha_c_frame - &projected_vel_alpha;

            let phi = beam.inner_angle_deg(&transverse_alpha_vel);
            //let phi = vel_c12_lab.inner_angle_deg(&beam);

            if beam.dot(&vel_c12_lab.cross(&transverse_alpha_vel)) < 0. {
                //if vel_c12_lab.dot(&vel_c12_lab.cross(&beam)) < 0. {
                phi
            } else {
                360. - phi
            }
        };

        (theta, phi)
    }

    #[allow(dead_code)]
    fn angles_method_2(c12: &Source, alpha: &Particle) -> (f64, f64) {
        //let c12 = c12.relativistic_source();
        let mut alpha = alpha.clone();
        let vel_alpha_lab_frame = alpha.classical_velocity_c();
        alpha.classical_boost(&c12.velocity_c());

        let beam = PhysVec::from_cartesian(0., 0., 1.);
        let vel_c12_lab = c12.velocity_c();
        let vel_alpha_c_frame = alpha.classical_velocity_c();

        let theta = vel_c12_lab.inner_angle_deg(&vel_alpha_c_frame);

        let phi = {
            let norm_1 = beam
                .cross(&vel_c12_lab)
                .as_normalized_to(1.)
                .rotated_around_vec(&vel_c12_lab, PI / 2.);

            let projected_alpha = vel_c12_lab.as_normalized_to(
                vel_alpha_lab_frame.mag() * vel_c12_lab.inner_angle_rad(&vel_alpha_lab_frame).cos(),
            );

            let transverse_alpha_vel = &vel_alpha_lab_frame - &projected_alpha;
            //let transverse_alpha_vel = &vel_alpha_c_frame - &projected_alpha;

            let phi = norm_1.inner_angle_deg(&transverse_alpha_vel);
            if vel_c12_lab.dot(&norm_1.cross(&vel_alpha_c_frame)) < 0. {
                phi
            } else {
                360. - phi
            }

            //phi -= 90.;
            //if phi < 0. {
            //phi + 360.
            //} else {
            //phi
            //}
            //phi
        };

        (theta, phi)
    }
    fn kinematically_accept_c12(&mut self, c12: &RelativisticSource) -> bool {
        let e_ex = c12.excitation_energy_MeV();

        let c12_mass = NUCLEAR_DB.nuclear_mass_MeV_per_c2(6, 12);
        let si28_mass = NUCLEAR_DB.nuclear_mass_MeV_per_c2(14, 28);
        let beam_kinetic_energy = 35.0 * c12_mass / 931.5;

        let c12_momentum =
            PhysVec::from_cartesian(0., 0., (beam_kinetic_energy * 2. * c12_mass).sqrt());

        let c12x_momentum = c12.source_momentum_MeV_per_c();
        let c12x_kinetic_energy = c12.kinetic_energy_MeV();
        let si28x_momentum = &c12_momentum - &c12x_momentum;
        let si28x_kinetic_energy = si28x_momentum.mag_sqr() / (2. * si28_mass);

        let missing_energy =
            beam_kinetic_energy - e_ex - si28x_kinetic_energy - c12x_kinetic_energy;
        self.e_missing_v_e_ex.fill(e_ex, missing_energy);

        missing_energy < 50.
    }

    #[allow(dead_code)]
    fn set_kinematic_gate(&mut self, value: f64) -> &mut Self {
        self.missing_energy_gate = value;
        self
    }

    #[allow(dead_code)]
    fn with_kinematic_gate(mut self, value: f64) -> Self {
        self.missing_energy_gate = value;
        self
    }
}

impl AnlModule<Event> for AaaThruBe8gs {
    fn name(&self) -> String {
        "aaa_thru_be8gs".to_string()
    }

    fn analyze_event(&mut self, event: &Event, _idx: usize) {
        let rules = [ParticleSelectionRule::new(3, 2, 4)];
        let combo_iter = AllCombinationsIter::new(&rules, event);

        for source_3a in combo_iter {
            let source_3a_rel = source_3a.relativistic_source();
            let e_ex = source_3a_rel.excitation_energy_MeV();

            if !self.kinematically_accept_c12(&source_3a_rel) {
                continue;
            }

            self.e_ex.fill(e_ex);

            let decay_path = decay_pathway(&source_3a);
            let _clear_decay_path = decay_path.is_some();

            let mut rng = rand::thread_rng();
            let mut be8gs_found = false;
            match decay_path {
                Some((Be8::SS0p, be8, alpha)) => {
                    self.e_ex_thru_be8_gs.fill(e_ex);

                    let c12_is_3m = (9.0..10.2).contains(&e_ex);
                    let c12_is_0p = (0.0..8.0).contains(&e_ex);
                    let c12_is_high_energy = (25.0..100.0).contains(&e_ex);

                    let (theta, phi) = Self::angles_method_1(&source_3a, alpha);
                    //let (theta, phi) = Self::theta_phi_method_2(&source_3a, alpha);
                    if c12_is_3m {
                        self.phi_theta_1_3m.fill(theta, phi);
                    }
                    if c12_is_0p {
                        self.phi_theta_0p.fill(theta, phi);
                    }
                    if c12_is_high_energy {
                        self.phi_theta_high_ex.fill(theta, phi);
                    }
                    let (theta, phi) = Self::angles_method_2(&source_3a, alpha);
                    if c12_is_3m {
                        self.phi_theta_2_3m.fill(theta, phi);
                    }

                    self.e_ex_2a_be8gs_fills_once.fill(be8.excitation_energy_MeV());
                    be8gs_found = true;


                    self.e_ex_3_vs_theta_thru_0p.fill(theta, e_ex);
                    self.e_ex_2a_thru_be8_any.fill(be8.excitation_energy_MeV());
                }
                Some((Be8::SS2p, be8, alpha)) => {
                    let (theta, _phi) = Self::angles_method_2(&source_3a, alpha);
                    self.ex_3_vs_theta_thru_2p.fill(theta, e_ex);
                    self.e_ex_thru_be_2p.fill(e_ex);
                    self.e_ex_2a_thru_be8_any.fill(be8.excitation_energy_MeV());
                }
                Some((Be8::SS4p, be8, alpha)) => {
                    let (theta, _phi) = Self::angles_method_2(&source_3a, alpha);
                    self.e_ex_3_vs_theta_thru_4p.fill(theta, e_ex);
                    self.e_ex_2a_thru_be8_any.fill(be8.excitation_energy_MeV());
                }
                _ => {}
            }

            let mut be8_gs_found = false;
            let mut be8_2p_found = false;
            let mut indices = (0..3).collect_vec();

            let mut min_e_ex_2 = f64::MAX;
            indices.shuffle(&mut rng);
            indices.iter().combinations(2).for_each(|combo_2| {
                let combo_2_particles = combo_2
                    .iter()
                    .map(|&idx| source_3a.particles()[*idx])
                    .collect_vec();

                let mut be8 = Source::new();
                combo_2_particles.iter().for_each(|&particle| {
                    be8.add_particle(particle);
                });

                let e_ex_2 = be8.excitation_energy_MeV();
                let e_rel_2 = be8.relative_kinetic_energy_MeV();

                if e_ex_2 < min_e_ex_2 {
                    min_e_ex_2 = e_ex_2;
                }

                self.e_ex_2a.fill(e_ex_2);

                if (9.0..10.2).contains(&e_ex) {
                    self.e_ex_2a_from_c12_3m.fill(e_ex_2);
                } else {
                    self.e_ex_2a_not_from_c12_3m.fill(e_ex_2);
                }

                self.e_ex_2_v_3.fill(e_ex, e_ex_2);
                if !be8gs_found {
                    self.e_ex_2a_be8gs_fills_once.fill(e_ex_2);
                }

                if (0.0..0.3).contains(&e_rel_2) {
                    be8_gs_found = true;
                }
                if (2.0..4.0).contains(&e_rel_2) {
                    be8_2p_found = true;
                }
            });

            self.e_ex_small_2a.fill(min_e_ex_2);

            let c12_is_0p = (0.0..8.0).contains(&e_ex);
            let c12_is_3m = (9.0..10.2).contains(&e_ex);
            let c12_is_low = (8.0..9.0).contains(&e_ex);
            let mut energies = source_3a_rel
                .particles()
                .iter()
                .map(|p| p.kinetic_energy_MeV())
                .collect_vec();
            energies.shuffle(&mut rng);

            let e_sum = energies.iter().sum::<f64>();

            let (e_1, e_2, e_3) = (
                energies[0] / e_sum,
                energies[1] / e_sum,
                energies[2] / e_sum,
            );

            let dalitz_x = f64::sqrt(3.) * (e_2 - e_1);

            let dalitz_y = 2. * e_3 - e_1 - e_2;
            //if be8_gs_found {
            if true {
                if c12_is_0p {
                    self.dalitz_thru_c12_0p.fill(dalitz_x, dalitz_y);
                }
                if c12_is_3m {
                    self.dalitz_thru_c12_3m.fill(dalitz_x, dalitz_y);
                }
                if c12_is_low {
                    self.dalitz_thru_c12_low.fill(dalitz_x, dalitz_y);
                }
                self.dalitz_all.fill(dalitz_x, dalitz_y);
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
        self.e_missing_v_e_ex.write();
        self.e_ex_2a.write();
        self.e_ex_2a_be8gs_fills_once.write();
        self.e_ex_2a_from_c12_3m.write();
        self.e_ex_2a_not_from_c12_3m.write();
        self.e_ex_small_2a.write();
        self.dalitz_all.write();
        self.dalitz_thru_c12_low.write();
        self.dalitz_thru_c12_0p.write();
        self.dalitz_thru_c12_3m.write();
        self.e_ex_thru_be8_gs.write();
        self.e_ex_thru_be_2p.write();
        self.e_ex_3_vs_theta.write();
        self.e_ex_3_vs_theta_thru_4p.write();
        self.e_ex_2a_thru_be8_any.write();
        self.ex_3_vs_theta_thru_2p.write();
        self.e_ex_3_vs_theta_thru_0p.write();
        self.e_ex_2_vs_theta.write();
        self.phi_theta_1_3m.write();
        self.phi_theta_2_3m.write();
        self.phi_theta_0p.write();
        self.phi_theta_high_ex.write();

        self.e_ex_2_v_3.write();
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

    let mut be8_e_rels = be8_sources
        .into_iter()
        .map(|be8| {
            let e_rel = be8.relative_kinetic_energy_MeV();
            (be8, e_rel)
        })
        .collect::<Vec<_>>();

    be8_e_rels.sort_by(|(_, e1), (_, e2)| e1.partial_cmp(e2).unwrap());

    let mut found_be8: Option<Source> = None;
    let mut found_state: Option<Be8> = None;
    //let mut found_second_state: bool = false;

    for (source, energy) in be8_e_rels {
        let be0p = (0.0..0.3).contains(&energy);
        let be2p = (2.0..4.0).contains(&energy);
        let be4p = (9.0..13.0).contains(&energy);

        if be0p {
            found_state = Some(Be8::SS0p);
            found_be8 = Some(source);
            break;
        }
        if be2p {
            found_state = Some(Be8::SS2p);
            found_be8 = Some(source);
            break;
        }
        if be4p {
            found_state = Some(Be8::SS4p);
            found_be8 = Some(source);
            break;
        }
    }

    if found_be8.is_none() {
        None
    } else {
        let found_be8 = found_be8.unwrap();
        let ss = found_state.unwrap();
        let alpha = other_alpha(c12, &found_be8);
        Some((ss, found_be8, alpha))
    }
}

fn other_alpha<'a>(c12: &'a Source, be8: &Source) -> &'a Particle {
    c12.iter()
        .find(|&&alpha| !be8.iter().any(|&be8_alpha| std::ptr::eq(be8_alpha, alpha)))
        .unwrap()
}
