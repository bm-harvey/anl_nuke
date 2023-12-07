use clap::Parser;
use faust::event::Event;
use faust_anl::general_particle_selection::MatchingPattern;
use faust_anl::ParticleSelectionRule;
use faust_anl::RelativeEnergy;
use faust_anl::RelativeEnergyConfig;
use faust_anl::{
    GeneralParticleFilter, GeneralParticleMixer, RandomizeLabPhiAngles, ShuffledPhiMixer,
};
use itertools::Itertools;
use nukers::anl_module::EventGenerator;

use nukers::anl_module::{Anl, MixedEventMaximum};

/// Relative energy analysis for looking at ensembles of alpha particles
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    #[arg(short, long, default_value_t = String::from("3_2_4"))]
    particle_rules: String,

    /// Maximum number of real events to look at
    #[arg(short, long, default_value_t=usize::MAX)]
    events: usize,

    /// Update interval
    #[arg(short, long, default_value_t = 100_000)]
    update_interval: usize,

    /// Number of mixed events to look at. If no number is passed, a 1:1 scaling of mixed:real will
    /// be calculated
    #[arg(short, long)]
    mixed: Option<usize>,

    /// Count low lying substates that are potentially populated
    #[arg(long)]
    count_states: bool,

    /// Count low lying substates that are potentially populated
    #[arg(long)]
    shape: bool,

    /// Do the inner angle analysis  
    #[arg(long)]
    inner_angle: bool,

    #[arg(long)]
    skip_real: bool,

    #[arg(short, long, default_value_t = String::from("si28_c_35"))]
    system: String,

    /// Mixer method
    #[arg(long, default_value_t=String::from("none"))]
    mixer: String,

    /// Detector filter method
    #[arg(long, default_value_t=String::from("faust"))]
    det_filter: String,

    /// Filtering policy
    #[arg(short, long, default_value_t=String::from("min"))] 
    filtering_policy: String,
}

fn main() {
    let args = Args::parse();
    let start = std::time::Instant::now();

    let filter_method = filter_method_from_args(&args);

    let mut particle_rules = Vec::<ParticleSelectionRule>::new();
    args.particle_rules.split("|").for_each(|rule_str| {
        let rule = rule_str.split("_").collect_vec();

        particle_rules.push(ParticleSelectionRule::new(
            rule[0].parse::<usize>().unwrap(),
            rule[1].parse::<usize>().unwrap(),
            rule[2].parse::<usize>().unwrap(),
        ));
    });

    let mut filter = GeneralParticleFilter::new(filter_method);
    for rule in particle_rules.iter() {
        filter.add_particles(rule.mult(), rule.z(), rule.a());
    }

    let config = RelativeEnergyConfig::new()
        .with_count_states(args.count_states)
        .with_shape(args.shape);
    let mut relative_energy = RelativeEnergy::from_config(config);
    for rule in particle_rules.iter() {
        relative_energy.add_particles(rule.mult(), rule.z(), rule.a());
    }

    let config = RelativeEnergyConfig::new()
        .with_count_states(args.count_states)
        .with_shape(args.shape);
    let mut relative_energy_mixed = RelativeEnergy::from_config(config);
    for rule in particle_rules.iter() {
        relative_energy_mixed.add_particles(rule.mult(), rule.z(), rule.a());
    }

    let anl = Anl::<faust::event::Event>::new()
        .with_input_directory(&format!("/data/sjygroup/sjy20/bmharvey/acs/{}/rkyv", args.system))
        .with_output_directory(&format!("/data/sjygroup/sjy20/bmharvey/acs/{}/anl", args.system))
        .with_filter(filter)
        .with_filtered_output_size(2_000_000)
        .with_mixed_module(relative_energy_mixed)
        .use_existing_filtered(true)
        .with_update_interval(args.update_interval)
        .with_max_real_events(args.events);

    let anl = assign_mixer(&args, anl, &particle_rules);

    let anl = if args.skip_real {
        anl
    } else {
        anl.with_real_module(relative_energy)
    };

    let mut anl = match args.mixed {
        None => anl.with_max_mixed_events(MixedEventMaximum::Factor(1.0)),
        Some(value) => anl.with_max_mixed_events(MixedEventMaximum::Absolute(value)),
    };

    // run the analysis
    anl.run();

    println!("Analysis took {} seconds", start.elapsed().as_secs_f32());
}

fn assign_mixer<'a>(args: &'a Args, anl: Anl<'a, Event>, particle_rules : &[ParticleSelectionRule]) -> Anl<'a, Event> {
    match args.mixer.as_str() {
        "std" => {
            let mut mixer = GeneralParticleMixer::new();
            for rule in particle_rules.iter() {
                mixer.add_rule(rule.mult(), rule.z(), rule.a());
            }

            match args.det_filter.as_str() {
                "faust" => {
                    mixer.set_faust_filter();
                }
                "none" => {}
                cutoff => {
                    let cutoff = cutoff.parse::<f64>().unwrap();
                    mixer.set_inner_angle_filter(cutoff);
                }
            }
            anl.with_event_generator(EventGenerator::Mixer(Box::new(mixer)))
        }
        "phi" => {
            let mut mixer = RandomizeLabPhiAngles::new();
            match args.det_filter.as_str() {
                "faust" => {
                    mixer.set_faust_filter();
                }
                "none" => {}
                cutoff => {
                    let cutoff = cutoff.parse::<f64>().unwrap();
                    mixer.set_inner_angle_filter(cutoff);
                }
            }
            anl.with_event_generator(EventGenerator::Scrambler(Box::new(mixer)))
        }
        "sh_phi" => {
            let mixer = ShuffledPhiMixer::new();
            anl.with_event_generator(EventGenerator::Mixer(Box::new(mixer)))
        }
        _ => anl,
    }
}

fn filter_method_from_args(args: &Args) -> MatchingPattern {
    match args.filtering_policy.as_str() {
        "std" => MatchingPattern::Standard,
        "strict" => MatchingPattern::Strict,
        "min" => MatchingPattern::Minimum,
        _ => {
            panic!("The set filtering policy is not valid. Options are `std`, `strict`, `min`");
        }
    }
}
