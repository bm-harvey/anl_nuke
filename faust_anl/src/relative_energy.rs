use crate::general_particle_selection::ParticleSelectionRule;
use crate::source::Source;
use cxx::{SharedPtr, UniquePtr};
use faust::event::Event;
use faust::event::Particle;
use itertools::Itertools;
use nukers::anl_module::AnlModule;
use rand::seq::SliceRandom;
use roost::file::RtFile;
use roost::hist::{RtH1D, RtH2D};
use roost::tree::{RtBranch, RtIntBranch, RtTree};
use std::f64::consts::PI;
use std::path::Path;

#[derive(Clone, Debug, Default)]
pub struct RelativeEnergyConfig {
    count_states: bool,
    shape: bool,
    inner_angle: bool,
}

impl RelativeEnergyConfig {
    pub fn new() -> Self {
        Self::default()
    }
    pub fn with_count_states(mut self, count_states: bool) -> Self {
        self.count_states = count_states;
        self
    }
    pub fn with_shape(mut self, shape: bool) -> Self {
        self.shape = shape;
        println!("shape: {}", self.shape);
        self
    }
}

pub struct RelativeEnergy {
    particle_rules: Vec<ParticleSelectionRule>,

    //file
    out_file: Option<SharedPtr<RtFile>>,

    // spectra for output
    e_rel_hist: SharedPtr<RtH1D>,
    inner_angle_avg_hist: SharedPtr<RtH1D>,
    source_particles_rel_ke_hist: SharedPtr<RtH1D>,
    source_particles_rel_p_hist: SharedPtr<RtH1D>,
    splatter_hist: SharedPtr<RtH2D>,

    //tree
    tree: SharedPtr<RtTree>,
    e_rel_branch: UniquePtr<RtBranch>,
    source_lab_ke: UniquePtr<RtBranch>,

    inner_angle_avg_branch: Option<UniquePtr<RtBranch>>,

    sphericity: Option<UniquePtr<RtBranch>>,
    coplanarity: Option<UniquePtr<RtBranch>>,
    shape_hist: Option<SharedPtr<RtH2D>>,

    be8gs_count: Option<UniquePtr<RtIntBranch>>,
    hoyle_count: Option<UniquePtr<RtIntBranch>>,
    c12_3m_count: Option<UniquePtr<RtIntBranch>>,

    source_particles_rel_p_w_1_8begs: Option<SharedPtr<RtH1D>>,
    source_particles_rel_p_w_1_hoyle: Option<SharedPtr<RtH1D>>,
    source_particles_rel_p_w_2_8begs: Option<SharedPtr<RtH1D>>,
    source_particles_rel_p_w_only_2_8begs: Option<SharedPtr<RtH1D>>,
    source_particles_rel_p_w_only_2_8begs_v_e_rel: Option<SharedPtr<RtH2D>>,
    source_particles_rel_p_w_1_12c3m: Option<SharedPtr<RtH1D>>,

    config: RelativeEnergyConfig,
}

impl RelativeEnergy {
    const BE8GS_THRESHOLD_REL: f64 = 0.20;
    const HOYLE_THRESHOLD_REL: f64 = 0.61;
    const C12_3M_LOWER_REL: f64 = 2.1;
    const C12_3M_UPPER_REL: f64 = 3.;

    pub fn new() -> Self {
        Self::from_config(RelativeEnergyConfig::default())
    }

    pub fn from_config(config: RelativeEnergyConfig) -> Self {
        let tree = roost::tree::new_tree("T".into(), "T".into());

        let config_clone = config.clone();
        Self {
            config,

            particle_rules: Vec::new(),

            // hist
            e_rel_hist: roost::hist::new_h1d(
                "e_rel_hist".into(),
                format!("{};Erel (MeV);Yield", ""),
                4 * 4_096,
                -5.,
                200.,
            ),

            inner_angle_avg_hist: roost::hist::new_h1d(
                "phi_rel_hist".into(),
                format!("{};avg. inner angle #phi [MeV] ;Yield", ""),
                2_048,
                0.,
                360.,
            ),

            source_particles_rel_ke_hist: roost::hist::new_h1d(
                "source_particles_rel_ke_hist".into(),
                format!("{};K.E. of all Particles relative to COM [MeV];Yield", ""),
                2_048,
                0.,
                100.,
            ),

            source_particles_rel_p_hist: roost::hist::new_h1d(
                "source_particles_rel_p_hist".into(),
                format!(
                    "{};Momentum of all particles relative to COM [MeV];Yield",
                    ""
                ),
                2_048,
                0.,
                600.,
            ),
            source_particles_rel_p_w_1_8begs: if config_clone.count_states {
                Some(roost::hist::new_h1d(
                    "source_particles_rel_p_w_1_8begs".into(),
                    format!(
                        "1 {{}}^{{8}}Be g.s.;Momentum of all particles relative to COM [MeV];Yield"
                    ),
                    2_048,
                    0.,
                    600.,
                ))
            } else {
                None
            },

            source_particles_rel_p_w_2_8begs: if config_clone.count_states {
                Some(roost::hist::new_h1d(
                    "source_particles_rel_p_w_2_8begs".into(),
                    format!(
                        "2 {{}}^{{8}}Be g.s.;Momentum of all particles relative to COM [MeV];Yield"
                    ),
                    2_048,
                    0.,
                    600.,
                ))
            } else {
                None
            },

            source_particles_rel_p_w_only_2_8begs: if config_clone.count_states {
                Some(roost::hist::new_h1d(
                    "source_particles_rel_p_w_only_2_8begs".into(),
                    format!(
                        "2 {{}}^{{8}}Be g.s. && no Hoyle or {{}}^{{12}}C(3-);Momentum of all particles relative to COM [MeV];Yield"
                    ),
                    2_048,
                    0.,
                    600.,
                ))
            } else {
                None
            },

            source_particles_rel_p_w_only_2_8begs_v_e_rel: if config_clone.count_states {
                Some(roost::hist::new_h2d(
                    "source_particles_rel_p_w_only_2_8begs_v_e_rel".into(),
                    format!(
                        "2 {{}}^{{8}}Be g.s. && no Hoyle or {{}}^{{12}}C(3-);E_{{rel}} [MeV];Momentum of all particles relative to COM [MeV];Yield"
                    ),
                    1_024,
                    0.,
                    200.,
                    1_024,
                    0.,
                    600.,
                ))
            } else {
                None
            },

            source_particles_rel_p_w_1_12c3m: if config_clone.count_states {
                Some(roost::hist::new_h1d(
                    "source_particles_rel_p_w_1_c12_3m".into(),
                    format!(
                        "1 {{}}^{{12}}C(3-).;Momentum of all particles relative to COM [MeV];Yield"
                    ),
                    2_048,
                    0.,
                    600.,
                ))
            } else {
                None
            },

            source_particles_rel_p_w_1_hoyle: if config_clone.count_states {
                Some(roost::hist::new_h1d(
                    "source_particles_rel_p_w_1_hoyle".into(),
                    format!(
                        "1 {{}}^{{12}}C(0+).;Momentum of all particles relative to COM [MeV];Yield"
                    ),
                    2_048,
                    0.,
                    600.,
                ))
            } else {
                None
            },

            splatter_hist: roost::hist::new_h2d(
                "splatter".into(),
                format!("{};#thetacos#varphi;#thetasin#varphi", ""),
                1_024,
                -45.,
                45.,
                1_024,
                -45.,
                45.,
            ),

            shape_hist: if config_clone.shape {
                Some(roost::hist::new_h2d(
                    "shape".into(),
                    format!("{};Sphericity;Coplanarity", ""),
                    256,
                    0.,
                    1.,
                    256,
                    0.,
                    3_f64.sqrt() / 4.,
                ))
            } else {
                None
            },

            e_rel_branch: tree.make_branch("e_rel".into()),
            sphericity: if config_clone.shape {
                Some(tree.make_branch("sphericity".into()))
            } else {
                None
            },
            coplanarity: if config_clone.shape {
                Some(tree.make_branch("coplanarity".into()))
            } else {
                None
            },

            source_lab_ke: tree.make_branch("source_lab_ke".into()),

            inner_angle_avg_branch: if config_clone.inner_angle {
                Some(tree.make_branch("inner_angle_avg".into()))
            } else {
                None
            },

            be8gs_count: if config_clone.count_states {
                Some(tree.make_int_branch("be8gs_count".into()))
            } else {
                None
            },

            hoyle_count: if config_clone.count_states {
                Some(tree.make_int_branch("hoyle_count".into()))
            } else {
                None
            },

            c12_3m_count: if config_clone.count_states {
                Some(tree.make_int_branch("c12_3m_count".into()))
            } else {
                None
            },
            tree,
            out_file: None,
        }
    }

    pub fn add_particles(&mut self, mult: usize, z: usize, a: usize) -> &mut Self {
        for rule in self.particle_rules.iter_mut() {
            if rule.z == z && rule.a == a {
                rule.mult += mult;
                return self;
            }
        }

        self.particle_rules
            .push(ParticleSelectionRule::new(mult, z, a));
        self
    }

    pub fn with_particles(mut self, mult: usize, z: usize, a: usize) -> Self {
        self.add_particles(mult, z, a);
        self
    }

    pub fn event_passes_rules(&self, event: &Event) -> bool {
        self.particle_rules
            .iter()
            .all(|rule| rule.exactly_satisfied(event))
    }
}
pub fn event_passes_rules(event: &Event, rules: &[ParticleSelectionRule]) -> bool {
    rules.iter().all(|rule| rule.exactly_satisfied(event))
}

pub fn construct_source<'a>(
    event: &'a Event,
    particle_rules: &'a Vec<ParticleSelectionRule>,
) -> Source<'a> {
    let mut counters = particle_rules
        .iter()
        .map(|rule| rule.mult)
        .collect::<Vec<_>>();

    let mut source = Source::new();

    let mut indices = (0..event.mult()).collect::<Vec<_>>();
    indices.shuffle(&mut rand::thread_rng());

    'particles: for idx in indices.iter() {
        let particle = &event.particles()[*idx];

        for rule_idx in 0..particle_rules.len() {
            let counter = &mut counters[rule_idx];
            if *counter == 0 {
                continue;
            }

            let rule = &particle_rules[rule_idx];

            if particle.is_pid(rule.z, rule.a) {
                *counter -= 1;
                source.add_particle(particle);
                let remaining_particles = counters.iter().sum::<usize>();
                if remaining_particles == 0 {
                    break 'particles;
                } else {
                    continue 'particles;
                }
            }
        }
    }
    source
}

impl<'a> RelativeEnergy {
    pub fn construct_source<'b>(
        rules: &'b [ParticleSelectionRule],
        event: &'a Event,
    ) -> Source<'a> {
        let mut counters = rules.iter().map(|rule| rule.mult).collect::<Vec<_>>();

        let mut source = Source::new();

        let mut indices = (0..event.mult()).collect::<Vec<_>>();
        indices.shuffle(&mut rand::thread_rng());

        'particles: for idx in indices.iter() {
            let particle = &event.particles()[*idx];

            for rule_idx in 0..rules.len() {
                let counter = &mut counters[rule_idx];
                if *counter == 0 {
                    continue;
                }

                let rule = &rules[rule_idx];

                if particle.is_pid(rule.z, rule.a) {
                    *counter -= 1;
                    source.add_particle(particle);
                    let remaining_particles = counters.iter().sum::<usize>();
                    if remaining_particles == 0 {
                        break 'particles;
                    } else {
                        continue 'particles;
                    }
                }
            }
        }
        source
    }

    fn shape_analysis(&mut self, source: &Source) {
        let (sphericity, coplanarity) = if self.config.shape {
            source.shape()
        } else {
            (0., 0.)
        };
        if self.config.shape {
            self.shape_hist
                .as_ref()
                .unwrap()
                .fill(sphericity, coplanarity);
            self.sphericity
                .as_mut()
                .unwrap()
                .pin_mut()
                .set_value(sphericity);
            self.coplanarity
                .as_mut()
                .unwrap()
                .pin_mut()
                .set_value(coplanarity);
        }

        for energy in source.relative_kinetic_energies_MeV() {
            self.source_particles_rel_ke_hist.fill(energy);
        }
        for momentum in source.relative_momenta_MeV_per_c() {
            self.source_particles_rel_p_hist.fill(momentum);
        }
    }
    fn count_states_analysis(&mut self, source: &Source) {
        if self.config.count_states {
            // count 8Begs
            let be8gs_count = source
                .iter()
                .combinations(2)
                .filter(|combo| {
                    let mut potential_be8gs = Source::new();
                    for particle in combo {
                        potential_be8gs.add_particle(particle);
                    }
                    potential_be8gs.relative_kinetic_energy_MeV()
                        < RelativeEnergy::BE8GS_THRESHOLD_REL
                })
                .count();

            // count Hoyle
            let hoyle_count = source
                .iter()
                .combinations(3)
                .filter(|combo| {
                    let mut potential_hoyle = Source::new();
                    for particle in combo {
                        potential_hoyle.add_particle(particle);
                    }

                    potential_hoyle.relative_kinetic_energy_MeV()
                        < RelativeEnergy::HOYLE_THRESHOLD_REL
                })
                .count();

            // count 12 C (3-)
            let c12_3m_count = source
                .iter()
                .combinations(3)
                .filter(|combo| {
                    let mut potential_state = Source::new();
                    for particle in combo {
                        potential_state.add_particle(particle);
                    }
                    let e_rel = potential_state.relative_kinetic_energy_MeV();
                    RelativeEnergy::C12_3M_LOWER_REL < e_rel
                        && e_rel < RelativeEnergy::C12_3M_UPPER_REL
                })
                .count();

            self.be8gs_count
                .as_mut()
                .unwrap()
                .pin_mut()
                .set_value(be8gs_count as i64);

            self.hoyle_count
                .as_mut()
                .unwrap()
                .pin_mut()
                .set_value(hoyle_count as i64);

            self.c12_3m_count
                .as_mut()
                .unwrap()
                .pin_mut()
                .set_value(c12_3m_count as i64);

            if be8gs_count == 1 {
                for &momentum in source.relative_momenta_MeV_per_c().iter() {
                    self.source_particles_rel_p_w_1_8begs
                        .as_mut()
                        .unwrap()
                        .fill(momentum);
                }
            }

            if be8gs_count == 2 {
                for &momentum in source.relative_momenta_MeV_per_c().iter() {
                    self.source_particles_rel_p_w_2_8begs
                        .as_mut()
                        .unwrap()
                        .fill(momentum);
                }

                if c12_3m_count == 0 && hoyle_count == 0 {
                    for &momentum in source.relative_momenta_MeV_per_c().iter() {
                        self.source_particles_rel_p_w_only_2_8begs
                            .as_mut()
                            .unwrap()
                            .fill(momentum);
                        self.source_particles_rel_p_w_only_2_8begs_v_e_rel
                            .as_mut()
                            .unwrap()
                            .fill(source.relative_kinetic_energy_MeV(), momentum);
                    }
                }
            }

            if c12_3m_count == 1 {
                for &momentum in source.relative_momenta_MeV_per_c().iter() {
                    self.source_particles_rel_p_w_1_12c3m
                        .as_mut()
                        .unwrap()
                        .fill(momentum);
                }
            }
            if hoyle_count == 1 {
                for &momentum in source.relative_momenta_MeV_per_c().iter() {
                    self.source_particles_rel_p_w_1_hoyle
                        .as_mut()
                        .unwrap()
                        .fill(momentum);
                }
            }
        }
    }

    fn inner_angle_analysis(&mut self, source: &Source) {
        if self.config.inner_angle {
            // averaging the inner angles between every particle in the source
            let mut inner_angle_avg = 0.;
            for idx_1 in 0..source.mult() {
                for idx_2 in (idx_1 + 1)..source.mult() {
                    let p_1 = source.particles()[idx_1].momentum_MeV_per_c();
                    let p_2 = source.particles()[idx_2].momentum_MeV_per_c();
                    inner_angle_avg += p_1.inner_angle_deg(p_2);
                }
            }
            let num_angles = source.mult() * (source.mult() - 1) / 2;
            inner_angle_avg /= num_angles as f64;

            self.inner_angle_avg_hist.fill(inner_angle_avg);

            self.inner_angle_avg_branch
                .as_mut()
                .unwrap()
                .pin_mut()
                .set_value(inner_angle_avg);
        }
    }
}

pub struct AllCombinationsIter<'a> {
    combos: Vec<Vec<Vec<&'a Particle>>>,
    current_indices: Vec<usize>,
    current_status: bool,
}

impl<'a> AllCombinationsIter<'a> {
    pub fn new(rules: &[ParticleSelectionRule], event: &'a Event) -> Self {
        let mut combos = Vec::with_capacity(rules.len());
        let mut current_indices = Vec::with_capacity(rules.len());

        let mut current_status = true;
        for rule in rules {
            current_indices.push(0);
            let pid_combo = event
                .iter_by_pid(rule.z, rule.a)
                .combinations(rule.mult)
                .collect_vec();
            if pid_combo.is_empty() {
                current_status = false;
            }
            combos.push(pid_combo);
        }

        Self {
            combos,
            current_indices,
            current_status,
        }
    }

    fn increment_index(&mut self, idx: usize) -> bool {
        if idx >= self.current_indices.len() {
            return false;
        }
        if self.current_indices[idx] + 1 >= self.combos[idx].len() {
            self.current_indices[idx] = 0;
            self.increment_index(idx + 1)
        } else {
            self.current_indices[idx] += 1;

            true
        }
    }
}

impl<'a> Iterator for AllCombinationsIter<'a> {
    type Item = Source<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.current_status {
            let mut result = Source::new();
            for idx in 0..self.current_indices.len() {
                let combo_idx = self.current_indices[idx];

                for particle in self.combos[idx][combo_idx].iter() {
                    result.add_particle(particle);
                }
            }

            self.current_status = self.increment_index(0);
            return Some(result);
        }
        None
    }
}

impl Default for RelativeEnergy {
    fn default() -> Self {
        Self::new()
    }
}

impl AnlModule<Event> for RelativeEnergy {
    fn name(&self) -> String {
        let mut name = String::from("relative_energy");
        for rule in self.particle_rules.iter() {
            name.push_str(format!("_{}", rule.to_string()).as_str());
        }
        name
    }
    fn initialize(&mut self, _output_directory: &Path) {
        let file_name: String = _output_directory
            .join(format!("{}.root", self.name()))
            .to_str()
            .unwrap()
            .into();
        println!("Creating output file: {}", file_name);
        self.out_file = Some(roost::file::create(file_name));
    }

    fn analyze_event(&mut self, event: &Event, _idx: usize) {
        let rules = &self.particle_rules;
        let all_combos_iter = AllCombinationsIter::new(&rules, event);
        for source in all_combos_iter {
            self.shape_analysis(&source);
            self.count_states_analysis(&source);
            self.inner_angle_analysis(&source);

            // basic properties of the source
            let e_rel = source.relative_kinetic_energy_MeV();
            let source_kinetic_energy = source.relative_kinetic_energy_MeV();

            self.e_rel_hist.fill(e_rel);

            for particle in source.particles() {
                let theta = particle.momentum_MeV_per_c().theta_deg();
                let phi = particle.momentum_MeV_per_c().phi_rad();
                self.splatter_hist
                    .fill(theta * phi.cos(), theta * phi.sin());
            }

            self.e_rel_branch.pin_mut().set_value(e_rel);
            self.source_lab_ke
                .pin_mut()
                .set_value(source_kinetic_energy);

            self.out_file.as_ref().unwrap().cd();
            self.tree.fill();
        }
    }

    fn generate_output(&mut self, _output_directory: &Path) {
        self.out_file.as_ref().unwrap().cd();
        self.e_rel_hist.write();
        self.inner_angle_avg_hist.write();
        self.splatter_hist.write();
        if self.config.shape {
            self.shape_hist.as_mut().unwrap().write();
            self.source_particles_rel_ke_hist.write();
            self.source_particles_rel_p_hist.write();
            if self.config.shape && self.config.count_states {
                self.source_particles_rel_p_w_1_12c3m
                    .as_ref()
                    .unwrap()
                    .write();
                self.source_particles_rel_p_w_1_8begs
                    .as_ref()
                    .unwrap()
                    .write();
                self.source_particles_rel_p_w_2_8begs
                    .as_ref()
                    .unwrap()
                    .write();
                self.source_particles_rel_p_w_only_2_8begs
                    .as_ref()
                    .unwrap()
                    .write();
                self.source_particles_rel_p_w_only_2_8begs_v_e_rel
                    .as_ref()
                    .unwrap()
                    .write();
                self.source_particles_rel_p_w_1_hoyle
                    .as_ref()
                    .unwrap()
                    .write();
            }
        }
        self.tree.write();
        self.out_file.as_ref().unwrap().close();
    }
}

pub struct RelativeEnergyAfterScrambling {
    particle_rules: Vec<ParticleSelectionRule>,
    e_rel_real_v_rndm_phi: SharedPtr<RtH2D>,
    out_file: Option<SharedPtr<RtFile>>,
    e_rel_real_hist: SharedPtr<RtH1D>,
    e_rel_scrambled_hist: SharedPtr<RtH1D>,
    e_rel_scrambled_similar_hist: SharedPtr<RtH1D>,
    e_rel_real_similar_hist: SharedPtr<RtH1D>,
    e_rel_rescrambled_match_hist: SharedPtr<RtH1D>,
}

impl RelativeEnergyAfterScrambling {
    pub fn new() -> Self {
        let e_rel_real_v_rndm_phi = roost::hist::new_h2d(
            "e_rel_real_v_rndm_phi".into(),
            "e_rel_real_v_rndm_phi".into(),
            1_024,
            0.,
            200.,
            1_024,
            0.,
            200.,
        );

        let e_rel_real_hist = roost::hist::new_h1d(
            "e_rel_real".into(),
            "e_rel_real".into(),
            1_024 * 8,
            0.,
            200.,
        );

        let e_rel_scrambled_hist = roost::hist::new_h1d(
            "e_rel_mixed".into(),
            "e_rel_mixed".into(),
            1_024 * 8,
            0.,
            200.,
        );
        let e_rel_real_similar_hist = roost::hist::new_h1d(
            "e_rel_real_similar".into(),
            "e_rel_real_similar".into(),
            1_024 * 8,
            0.,
            200.,
        );

        let e_rel_scrambled_similar_hist = roost::hist::new_h1d(
            "e_rel_scrambled_similar".into(),
            "e_rel_scrambled_similar".into(),
            1_024 * 8,
            0.,
            200.,
        );

        let e_rel_rescrambled_match_hist = roost::hist::new_h1d(
            "e_rel_rescrambled_match".into(),
            "e_rel_rescrambled_match".into(),
            1_024 * 8,
            0.,
            200.,
        );

        Self {
            particle_rules: Vec::new(),
            e_rel_real_v_rndm_phi,
            out_file: None,
            e_rel_real_hist,
            e_rel_scrambled_hist,
            e_rel_scrambled_similar_hist,
            e_rel_real_similar_hist,
            e_rel_rescrambled_match_hist,
        }
    }

    pub fn add_rule(&mut self, mult: usize, z: usize, a: usize) -> &mut Self {
        for rule in self.particle_rules.iter_mut() {
            if rule.z == z && rule.a == a {
                rule.mult += mult;
                return self;
            }
        }

        self.particle_rules
            .push(ParticleSelectionRule::new(mult, z, a));
        self
    }

    pub fn with_particles(mut self, mult: usize, z: usize, a: usize) -> Self {
        self.add_rule(mult, z, a);
        self
    }
}

impl AnlModule<Event> for RelativeEnergyAfterScrambling {
    fn name(&self) -> String {
        let mut name = String::from("relative_energy_after_scrambling");
        for rule in self.particle_rules.iter() {
            name.push_str(format!("_{}", rule.to_string()).as_str());
        }
        name
    }

    fn initialize(&mut self, _output_directory: &Path) {
        let file_name: String = _output_directory
            .join(format!("{}.root", self.name()))
            .to_str()
            .unwrap()
            .into();
        println!("Creating output file: {}", file_name);
        self.out_file = Some(roost::file::create(file_name));
    }

    fn analyze_event(&mut self, _event: &Event, _idx: usize) {
        let rules = &self.particle_rules;

        let all_combos_iter = AllCombinationsIter::new(&rules, _event);
        let mut e_rels_real = Vec::new();
        for source in all_combos_iter {
            let e_rel = source.relative_kinetic_energy_MeV();
            e_rels_real.push(e_rel);
            self.e_rel_real_hist.fill(e_rel);
        }

        let num_iterations = 1;
        for _idx in 0..num_iterations {
            // rndm phi
            let mut scrambled_event = _event.clone();
            for particle in scrambled_event.iter_mut() {
                particle.set_phi_rad(rand::random::<f64>() * 2. * PI);
            }

            let mut rescrambled_event = scrambled_event.clone();
            for particle in rescrambled_event.iter_mut() {
                particle.set_phi_rad(rand::random::<f64>() * 2. * PI);
            }

            let all_combos_iter = AllCombinationsIter::new(&rules, &scrambled_event);
            let mut e_rels_mixed = Vec::new();
            for source in all_combos_iter {
                let e_rel = source.relative_kinetic_energy_MeV();
                e_rels_mixed.push(e_rel);
            }

            let all_combos_iter = AllCombinationsIter::new(&rules, &rescrambled_event);
            let mut e_rels_remixed = Vec::new();
            for source in all_combos_iter {
                let e_rel = source.relative_kinetic_energy_MeV();
                e_rels_remixed.push(e_rel);
            }

            for (&e_rel_real, &e_rel_mixed) in e_rels_real.iter().zip(e_rels_mixed.iter()) {
                self.e_rel_real_v_rndm_phi.fill(e_rel_mixed, e_rel_real);

                let range_low = e_rel_real * (0.99);
                let range_up = e_rel_real * (1.01);

                self.e_rel_scrambled_hist.fill(e_rel_mixed);
                if e_rel_mixed > range_low && e_rel_mixed < range_up {
                    self.e_rel_scrambled_similar_hist.fill(e_rel_mixed);
                    self.e_rel_real_similar_hist.fill(e_rel_real);
                }
            }
            for (&e_rel_remixed, &e_rel_mixed) in e_rels_remixed.iter().zip(e_rels_mixed.iter()) {
                let range_low = e_rel_remixed * (0.99);
                let range_up = e_rel_remixed * (1.01);

                if e_rel_mixed > range_low && e_rel_mixed < range_up {
                    self.e_rel_rescrambled_match_hist.fill(e_rel_remixed);
                }
            }
        }
    }

    fn generate_output(&mut self, _output_directory: &Path) {
        self.out_file.as_ref().unwrap().cd();
        self.e_rel_real_v_rndm_phi.write();
        self.e_rel_scrambled_hist.write();
        self.e_rel_real_hist.write();
        self.e_rel_scrambled_similar_hist.write();
        self.e_rel_rescrambled_match_hist.write();
        self.out_file.as_ref().unwrap().close();
    }
}
