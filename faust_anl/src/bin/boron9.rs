use cxx::SharedPtr;
use faust::event::Event;
use faust_anl::general_particle_selection::GeneralParticleFilter;
use faust_anl::general_particle_selection::GeneralParticleMixer;
use faust_anl::general_particle_selection::ParticleSelectionRule;
use faust_anl::relative_energy::AllCombinationsIter;
use faust_anl::source::Source;
use nukers::anl_module::EventGenerator::Mixer;
use nukers::anl_module::MixedEventMaximum;
use nukers::anl_module::{Anl, AnlModule};
use rand::seq::SliceRandom;
use roost::hist::RtH1D;

fn main() {
    let filter = GeneralParticleFilter::new(faust_anl::MatchingPattern::Minimum)
        .with_particles(2, 2, 4)
        .with_particles(1, 1, 1);

    #[allow(unused_mut)]
    #[allow(unused_variables)]
    let mut mixer = GeneralParticleMixer::new()
        .with_particles(2, 2, 4)
        .with_particles(1, 1, 1);
    //mixer.set_faust_filter();

    let system = "c12_si_35";
    //let system = "si28_c_35";

    Anl::new()
        .with_filter(filter)
        .with_input_directory(&format!("/data/sjygroup/sjy20/bmharvey/acs/{}/rkyv", system))
        .with_output_directory(&format!("/data/sjygroup/sjy20/bmharvey/acs/{}/anl", system))
        .with_real_module(PaaThruBe8gs::new())
        .with_event_generator(Mixer(Box::new(mixer)))
        .with_mixed_module(PaaThruBe8gs::new())
        .with_update_interval(100_000)
        .with_max_mixed_events(MixedEventMaximum::Absolute(1_800_000))
        .run();
}

struct PaaThruBe8gs {
    rules: [ParticleSelectionRule; 2],
    e_ex: SharedPtr<RtH1D>,
    e_ex_thru_be8gs: SharedPtr<RtH1D>,
    e_ex_thru_li5gs: SharedPtr<RtH1D>,
    e_ex_aa: SharedPtr<RtH1D>,
    e_ex_pa_min: SharedPtr<RtH1D>,
    e_ex_pa_max: SharedPtr<RtH1D>,
    theta: SharedPtr<RtH1D>,
}

impl PaaThruBe8gs {
    fn new() -> Self {
        let rules = [
            ParticleSelectionRule::new(1, 1, 1),
            ParticleSelectionRule::new(2, 2, 4),
        ];

        Self {
            rules,
            #[rustfmt::skip]
            e_ex: roost::hist::new_h1d( String::from("e_ex"), String::from(";E_{x} [MeV];Yield"), 4_096, -2., 30.,),
            #[rustfmt::skip]
            theta: roost::hist::new_h1d( String::from("theta"), String::from(";#theta [deg];Yield"), 1_024, 0., 180.,),
            #[rustfmt::skip]
            e_ex_thru_be8gs: roost::hist::new_h1d( String::from("e_ex_thru_be8gs"), String::from(";E_{x} [MeV];Yield"), 4_096, -2., 30.,),
            #[rustfmt::skip]
            e_ex_pa_min: roost::hist::new_h1d( String::from("e_ex_pa_min"), String::from(";E_{x} [MeV];Yield"), 4_096, -2., 30.,),
            #[rustfmt::skip]
            e_ex_pa_max: roost::hist::new_h1d( String::from("e_ex_pa_max"), String::from(";E_{x} [MeV];Yield"), 4_096, -2., 30.,),
            #[rustfmt::skip]
            e_ex_aa: roost::hist::new_h1d( String::from("e_ex_aa"), String::from(";E_{x} [MeV];Yield"), 4_096, -2., 30.,),
            #[rustfmt::skip]
            e_ex_thru_li5gs: roost::hist::new_h1d( String::from("e_ex_thru_li5gs"), String::from(";E_{x} [MeV];Yield"), 4_096, -2., 30.,),
        }
    }
}

impl AnlModule<Event> for PaaThruBe8gs {
    fn name(&self) -> String {
        "paa_thru_be8gs".to_string()
    }

    fn analyze_event(&mut self, event: &Event, _idx: usize) {
        if event.pid_mult(1, 1) < 1 || event.pid_mult(2, 4) < 2 {
            return;
        }

        let combo_iter = AllCombinationsIter::new(&self.rules, event);

        for combo in combo_iter {
            let ex = combo.excitation_energy_MeV();

            let alphas = combo
                .iter()
                .filter(|particle| particle.is_pid(2, 4))
                .collect::<Vec<_>>();

            let proton = combo.iter().find(|particle| particle.is_pid(1, 1)).unwrap();

            let mut source_2a = Source::new();
            for particle in alphas.iter() {
                source_2a.add_particle(particle);
            }

            let e_rel_2a = source_2a.relative_kinetic_energy_MeV();

            self.e_ex_aa.fill(e_rel_2a);

            let mut source_pa = Source::new();
            let li5gs = alphas.iter().any(|&particle| {
                source_pa.clear();
                source_pa.add_particle(proton);
                source_pa.add_particle(particle);

                source_pa.relative_kinetic_energy_MeV() < 0.3
            });

            if 2. < ex && ex < 3. {
                let v = combo.velocity_c();

                let partciles_new_frame = combo
                    .iter()
                    .map(|particle| {
                        let mut particle_new_frame = (*particle).clone();
                        particle_new_frame.boost(&v.as_scaled_by(-1.0));
                        particle_new_frame
                    })
                    .collect::<Vec<_>>();

                //let proton = partciles_new_frame
                //.iter()
                //.find(|particle| particle.is_pid(1, 1))
                //.unwrap();
                let mut alphas = partciles_new_frame
                    .iter()
                    .filter(|particle| particle.is_pid(2, 4))
                    .collect::<Vec<_>>();
                alphas.shuffle(&mut rand::thread_rng());

                let v_rel = &alphas[0].velocity_c() - &alphas[1].velocity_c();
                let v_cm = &alphas[0].velocity_c() + &alphas[1].velocity_c();
                let theta = v_rel.inner_angle_deg(&v_cm);

                self.theta.fill(theta);
            }

            let be8gs = e_rel_2a < 0.2;
            if be8gs {
                self.e_ex_thru_be8gs.fill(ex);
            }
            if li5gs {
                self.e_ex_thru_li5gs.fill(ex);
            }

            self.e_ex.fill(ex);

            if (2.0..3.0).contains(&ex) {
                let pa_e_rel = alphas
                    .iter()
                    .map(|alpha| {
                        source_pa.clear();
                        source_pa.add_particle(proton);
                        source_pa.add_particle(alpha);
                        source_pa.relative_kinetic_energy_MeV()
                    })
                    //.min_by(|a, b| a.total_cmp(b))
                    .max_by(|a, b| a.total_cmp(b))
                    .unwrap();

                self.e_ex_pa_max.fill(pa_e_rel);

                let pa_e_rel = alphas
                    .iter()
                    .map(|alpha| {
                        source_pa.clear();
                        source_pa.add_particle(proton);
                        source_pa.add_particle(alpha);
                        source_pa.relative_kinetic_energy_MeV()
                    })
                    .min_by(|a, b| a.total_cmp(b))
                    .unwrap();
                self.e_ex_pa_min.fill(pa_e_rel);
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
        self.e_ex_pa_max.write();
        self.e_ex_aa.write();
        self.e_ex_pa_min.write();
        self.theta.write();
        self.e_ex_thru_be8gs.write();
        self.e_ex_thru_li5gs.write();
        file.close();
    }
}
