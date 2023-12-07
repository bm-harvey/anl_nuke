use faust::event::Event;
use faust_filter::faust_filter::FaustFilter;
use itertools::Itertools;
use nukers::anl_module::EventScrambler;
use nukers::anl_module::{EventFilter, EventMixer};
use nukers::data_set::ArchivedData;
use rand;
use rand::seq::SliceRandom;
use rand::Rng;
use rand_distr::{Distribution, Normal};
use std::f64::consts::PI;
use std::path::Path;

use crate::source::Source;

///
/// Forms the foundation for selecting and creating events
///
#[derive(Clone)]
pub struct ParticleSelectionRule {
    pub z: usize,
    pub a: usize,
    pub mult: usize,
}

impl ParticleSelectionRule {
    pub fn new(mult: usize, z: usize, a: usize) -> Self {
        Self { z, a, mult }
    }

    pub fn exactly_satisfied(&self, event: &Event) -> bool {
        event.pid_mult(self.z, self.a) == self.mult
    }

    pub fn at_least_satisfied(&self, event: &Event) -> bool {
        event.pid_mult(self.z, self.a) >= self.mult
    }

    pub fn z(&self) -> usize {
        self.z
    }
    pub fn a(&self) -> usize {
        self.a
    }
    pub fn mult(&self) -> usize {
        self.mult
    }
}

impl ToString for ParticleSelectionRule {
    fn to_string(&self) -> String {
        let particle_id = match (self.z, self.a) {
            (1, 1) => "p".to_string(),
            (1, 2) => "d".to_string(),
            (1, 3) => "t".to_string(),
            (2, 3) => "h".to_string(),
            (2, 4) => "a".to_string(),
            _ => format!("z{}a{}", self.z.clone(), self.a.clone()),
        };

        let result = format!("{}{}", self.mult, particle_id);

        result
    }
}

/// The pattern that the event must match to be accepted by the filter
pub enum MatchingPattern {
    /// The rules must be at least satisfied, and particles not referred to by the rules are ignored
    Minimum,
    /// The rules must be exactly satisfied, but other particles not specified by the rules are
    /// ignored
    Standard,
    /// The event must exactly match the set rule with no other measured particles
    Strict,
}

/// A filter that can be applied to events to remove events that do not meet certain criteria
/// There is a small "gotcha" here. `UniqueDetectors` will only check that the tagged detectors are
/// different. If angles of the particles are changed, consider use the `Faust` filter which
/// calculates .
pub enum GeometricFilter {
    None,
    Faust(FaustFilter),
    InnerAngleCut(f64),
    UniqueDetectors,
}
impl GeometricFilter {
    fn filter_event(&self, event: &Event) -> bool {
        match self {
            GeometricFilter::None => true,
            GeometricFilter::Faust(filter) => {
                let hit_detectors = event
                    .iter()
                    .map(|particle| filter.detector_hit(particle.p()))
                    .collect::<Vec<_>>();

                if hit_detectors.iter().any(|det| det.is_none()) {
                    return false;
                }

                hit_detectors.iter().all_unique()
            }
            GeometricFilter::InnerAngleCut(cutoff) => event.iter().combinations(2).all(|pair| {
                let (p1, p2) = (pair[0], pair[1]);
                let angle = p1.p().inner_angle_deg(p2.p());
                angle >= *cutoff
            }),
            GeometricFilter::UniqueDetectors => event
                .iter()
                .map(|particle| particle.detector())
                .all_unique(),
        }
    }
}

///
/// Gate on events with a certain selection of particles.  
///
pub struct GeneralParticleFilter {
    particle_rules: Vec<ParticleSelectionRule>,
    matching_pattern: MatchingPattern,
}

impl GeneralParticleFilter {
    pub fn new(matching_pattern: MatchingPattern) -> Self {
        Self {
            particle_rules: Vec::new(),
            matching_pattern,
        }
    }
    pub fn with_particles(mut self, mult: usize, z: usize, a: usize) -> Self {
        self.add_particles(mult, z, a);
        self
    }
    pub fn with_particles_tuple(mut self, properties: (usize, usize, usize)) -> Self {
        self.add_particles(properties.0, properties.1, properties.2);
        self
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
}

impl EventFilter<Event> for GeneralParticleFilter {
    fn name(&self) -> String {
        let mut name = String::from("filt");
        for rule in self.particle_rules.iter() {
            name.push_str(format!("_{}", rule.to_string()).as_str());
        }

        let name_suffix = match self.matching_pattern {
            MatchingPattern::Strict => "_strict",
            MatchingPattern::Standard => "_std",
            MatchingPattern::Minimum => "_min",
        };

        let name = format!("{}{}", name, name_suffix);
        name
    }

    fn filter_event(&self, event: &Event, _idx: usize) -> bool {
        let rules_mult = self.particle_rules.iter().map(|rule| rule.mult).sum();
        if event.mult() < rules_mult {
            return false;
        }

        match self.matching_pattern {
            MatchingPattern::Minimum => self
                .particle_rules
                .iter()
                .all(|rule| rule.at_least_satisfied(event)),

            MatchingPattern::Standard => self
                .particle_rules
                .iter()
                .all(|rule| rule.exactly_satisfied(event)),

            MatchingPattern::Strict => {
                self.particle_rules
                    .iter()
                    .all(|rule| rule.exactly_satisfied(event))
                    && self
                        .particle_rules
                        .iter()
                        .map(|rule| rule.mult)
                        .sum::<usize>()
                        == event.mult()
            }
        }
    }
}

/// The standard mixing technique. Grab each of the target particles from seperate events and
/// return that as the new event.
pub struct GeneralParticleMixer {
    particle_rules: Vec<ParticleSelectionRule>,
    filter: GeometricFilter,
}

impl Default for GeneralParticleMixer {
    fn default() -> Self {
        Self::new()
    }
}

impl GeneralParticleMixer {
    pub fn new() -> Self {
        Self {
            particle_rules: Vec::new(),
            filter: GeometricFilter::None,
        }
    }

    pub fn add_rule(&mut self, mult: usize, z: usize, a: usize) -> &mut Self {
        self.particle_rules
            .push(ParticleSelectionRule::new(mult, z, a));
        self
    }

    pub fn with_particles(mut self, mult: usize, z: usize, a: usize) -> Self {
        self.particle_rules
            .push(ParticleSelectionRule::new(mult, z, a));
        self
    }

    pub fn set_faust_filter(&mut self) -> &mut Self {
        self.filter = GeometricFilter::UniqueDetectors;
        self
    }

    pub fn set_inner_angle_filter(&mut self, inner_angle_cutoff: f64) -> &mut Self {
        self.filter = GeometricFilter::InnerAngleCut(inner_angle_cutoff);
        self
    }

    pub fn event_passes_rules(&self, event: &Event) -> bool {
        self.particle_rules
            .iter()
            .all(|rule| rule.exactly_satisfied(event))
    }

    pub fn detectors_are_unique(event: &Event) -> bool {
        let detectors = event.iter().map(|p| p.detector()).collect::<Vec<_>>();
        !(1..detectors.len()).any(|idx| detectors[idx..].contains(&detectors[idx - 1]))
    }
}

impl EventMixer<Event> for GeneralParticleMixer {
    fn name(&self) -> String {
        let mut name = String::from("mix");
        for rule in self.particle_rules.iter() {
            name.push_str(format!("_{}", rule.to_string()).as_str());
        }
        match self.filter {
            GeometricFilter::None => {}
            GeometricFilter::Faust(_) => name.push_str("_faust"),
            GeometricFilter::InnerAngleCut(value) => {
                name.push_str(&format!("_inner_angle_{:.2}", value))
            }
            GeometricFilter::UniqueDetectors => name.push_str("_unique"),
        }
        name
    }

    fn mix_events(&mut self, data_collection: &ArchivedData<Event>, _idx: usize) -> Option<Event> {
        let mut counters = self
            .particle_rules
            .iter()
            .map(|rule| rule.mult)
            .collect::<Vec<_>>();

        let mut result = Event::default();

        'counters: while counters.iter().sum::<usize>() > 0 {
            let event = data_collection.random_event();

            let mut indices = (0..event.mult()).collect::<Vec<_>>();

            indices.shuffle(&mut rand::thread_rng());
            for idx in indices.iter() {
                let particle = &event.particles()[*idx];

                for (rule_idx, rule) in self.particle_rules.iter().enumerate() {
                    let counter = &mut counters[rule_idx];
                    if *counter == 0 {
                        continue;
                    }

                    if particle.is_pid(rule.z, rule.a) {
                        *counter -= 1;
                        result.add_particle_by_clone(particle);
                        continue 'counters;
                    }
                }
            }
        }

        match self.filter.filter_event(&result) {
            true => Some(result),
            false => None,
        }
    }
}

pub struct RandomizeLabPhiAngles {
    name: String,
    particle_rules: Vec<ParticleSelectionRule>,
    enforce_unique_detectors: bool,
    use_filter: bool,
    filter: GeometricFilter,
}

impl RandomizeLabPhiAngles {
    pub fn new() -> Self {
        Self {
            name: "rndm_phi".into(),
            particle_rules: Vec::new(),
            enforce_unique_detectors: true,
            use_filter: false,
            filter: GeometricFilter::None,
        }
    }

    pub fn set_faust_filter(&mut self) -> &mut Self {
        self.filter = GeometricFilter::Faust(FaustFilter::new(Path::new(
            "..\\faust_filter\\param\\detnum_corner_z_x_y_r.txt",
        )));
        self
    }

    pub fn set_inner_angle_filter(&mut self, inner_angle_cutoff: f64) -> &mut Self {
        self.filter = GeometricFilter::InnerAngleCut(inner_angle_cutoff);
        self
    }

    pub fn add_rule(&mut self, mult: usize, z: usize, a: usize) -> &mut Self {
        self.particle_rules
            .push(ParticleSelectionRule::new(mult, z, a));
        self
    }

    pub fn with_rule(mut self, mult: usize, z: usize, a: usize) -> Self {
        self.particle_rules
            .push(ParticleSelectionRule::new(mult, z, a));
        self
    }
    pub fn use_filter(mut self, use_filter: bool) -> Self {
        self.use_filter = use_filter;
        self
    }

    pub fn enforce_unique_detectors(&mut self, enforce: bool) -> &mut Self {
        self.enforce_unique_detectors = enforce;
        self
    }

    pub fn event_passes_rules(&self, event: &Event) -> bool {
        self.particle_rules
            .iter()
            .map(|rule| rule.exactly_satisfied(event))
            .reduce(|r1, r2| r1 && r2)
            .unwrap()
    }

    pub fn detectors_are_unique(event: &Event) -> bool {
        let detectors = event.iter().map(|p| p.detector()).collect::<Vec<_>>();
        !(1..detectors.len()).any(|idx| detectors[idx..].contains(&detectors[idx - 1]))
    }
}

impl EventScrambler<Event> for RandomizeLabPhiAngles {
    fn name(&self) -> String {
        let mut result = self.name.clone();
        match self.filter {
            GeometricFilter::None => {}
            GeometricFilter::Faust(_) => result.push_str("_faust"),
            GeometricFilter::InnerAngleCut(value) => {
                result.push_str(&format!("_inner_angle_{:.2}", value))
            }
            GeometricFilter::UniqueDetectors => result.push_str("_unique"),
        }
        result
    }

    fn scramble_event(&mut self, event: &Event , _idx: usize) -> Option<Event> {
        let mut result = event.clone();

        result.iter_mut().for_each(|particle| {
            particle.set_phi_rad(rand::thread_rng().gen_range(0_f64..(2. * PI)));
        });

        let passes_filter: bool = self.filter.filter_event(&result);

        match passes_filter {
            true => Some(result),
            false => None,
        }
    }
}

pub struct ShuffledPhiMixer {
    name: String,
}

impl ShuffledPhiMixer {
    pub fn new() -> Self {
        Self {
            name: "shuffle_phi".into(),
        }
    }
}

impl EventMixer<Event> for ShuffledPhiMixer {
    fn name(&self) -> String {
        let result = self.name.clone();
        result
    }

    fn mix_events(&mut self, data_collection: &ArchivedData<Event>, idx: usize) -> Option<Event> {
        let real_event = data_collection
            .event_by_idx(idx % data_collection.len())
            .unwrap();
        let mut result = Event::default();

        let mut particles = real_event.particles().clone();

        particles.shuffle(&mut rand::thread_rng());

        for idx in 0..particles.len() {
            let mut new_particle = particles[idx].clone();
            let next_particle = &particles[(idx + 1) % particles.len()];
            new_particle.set_phi_rad(next_particle.p().phi_rad());
            result.add_particle(new_particle);
        }
        Some(result)
    }
}

pub struct MoveGaussianMixer {
    filter: GeometricFilter,
    distribution: Normal<f64>,
}

impl MoveGaussianMixer {
    pub fn new() -> Self {
        Self {
            filter: GeometricFilter::None,
            distribution: Normal::new(0., 0.1).unwrap(),
        }
    }

    pub fn set_faust_filter(&mut self) -> &mut Self {
        self.filter = GeometricFilter::Faust(FaustFilter::new(Path::new(
            "..\\faust_filter\\param\\detnum_corner_z_x_y_r.txt",
        )));
        self
    }

    pub fn set_inner_angle_filter(&mut self, inner_angle_cutoff: f64) -> &mut Self {
        self.filter = GeometricFilter::InnerAngleCut(inner_angle_cutoff);
        self
    }
}

impl EventMixer<Event> for MoveGaussianMixer {
    fn name(&self) -> String {
        let mut result = String::from("move_gauss");
        match self.filter {
            GeometricFilter::None => {}
            GeometricFilter::Faust(_) => result.push_str("_faust"),
            GeometricFilter::InnerAngleCut(value) => {
                result.push_str(&format!("_inner_angle_{:.2}", value))
            }
            GeometricFilter::UniqueDetectors => result.push_str("_unique"),
        }
        result
    }

    fn mix_events(&mut self, data_collection: &ArchivedData<Event>, idx: usize) -> Option<Event> {
        let real_event = data_collection
            .event_by_idx(idx % data_collection.len())
            .unwrap();
        let mut result = real_event.clone();

        result.iter_mut().for_each(|particle| {
            let p = particle.momentum_MeV_per_c().mag();
            let mut rng = rand::thread_rng();

            let dp: [f64; 3] = [
                self.distribution.sample(&mut rng) * 0.0 * p,
                self.distribution.sample(&mut rng) * 0.0 * p,
                self.distribution.sample(&mut rng) * 0.0 * p,
            ];

            let dp = faust::phys_vec::PhysVec::new(dp);

            let momentum = particle.momentum_MeV_per_c() + &dp;

            particle.set_momentum(&momentum);
        });

        let passes_filter: bool = self.filter.filter_event(&result);

        match passes_filter {
            true => Some(result),
            false => None,
        }
    }
}

pub struct RndmAngleInSourceFrameScrambler {}

impl RndmAngleInSourceFrameScrambler {
    pub fn new() -> Self {
        Self {}
    }
}

impl EventScrambler<Event> for RndmAngleInSourceFrameScrambler {
    fn name(&self) -> String {
        String::from("rndm_angle_in_source_frame")
    }
    fn scramble_event(&mut self, event: &Event, _idx: usize) -> Option<Event> {
        let mut rng = rand::thread_rng();
        let mut source = Source::new();
        let mut result = event.clone();

        result.iter().for_each(|particle| {
            source.add_particle(particle);
        });

        let source_velocity = source.velocity_c();

        result.iter_mut().for_each(|particle| {
            let rndm_phi = rng.gen_range(0_f64..(2. * PI));

            let unif_number = rng.gen_range(0_f64..1_f64);
            let rndm_theta = (1. - unif_number).acos();

            let mut momentum = particle.momentum_MeV_per_c() - &source_velocity;
            momentum.set_phi_rad(rndm_phi);
            momentum.set_theta_rad(rndm_theta);
            let momentum = &momentum + &source_velocity;
            particle.set_momentum(&momentum);
        });
        Some(result)
    }
}
