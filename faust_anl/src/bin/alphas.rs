use clap::Parser;
use faust::event::Event;
use nukers::anl_module::EventGenerator;
use faust_anl::general_particle_selection::MatchingPattern;
use faust_anl::RelativeEnergy;
use faust_anl::RelativeEnergyConfig;
use faust_anl::{GeneralParticleFilter, GeneralParticleMixer, RandomizeLabPhiAngles, ShuffledPhiMixer};

use nukers::anl_module::{Anl, MixedEventMaximum};

/// Relative energy analysis for looking at ensembles of alpha particles
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Number of alphas particles to look for in the events
    #[arg(short, long, default_value_t = 4)]
    num_alphas: usize,

    /// Number of alphas particles to study from the events
    #[arg(short = 'a', default_value_t = 4)]
    num_alphas_anl: usize,

    /// Filtering policy (std, strict, min)
    #[arg(short, long, default_value_t=String::from("std"))]
    filtering_policy: String,

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

    /// Mixer method
    #[arg(long, default_value_t=String::from("none"))]
    mixer: String,

    /// Detector filter method
    #[arg(long, default_value_t=String::from("faust"))]
    det_filter: String,
}

fn main() {
    let args = Args::parse();
    let start = std::time::Instant::now();

    let filter_method = filter_method_from_args(&args);

    let num = args.num_alphas;
    let num_anl = args.num_alphas_anl;

    let filter = GeneralParticleFilter::new(filter_method).with_particles(num, 2, 4);

    let config = RelativeEnergyConfig::new()
        .with_count_states(args.count_states)
        .with_shape(args.shape);
    let relative_energy = RelativeEnergy::from_config(config).with_particles(num_anl, 2, 4);

    let config = RelativeEnergyConfig::new()
        .with_count_states(args.count_states)
        .with_shape(args.shape);
    let relative_energy_mixed = RelativeEnergy::from_config(config).with_particles(num_anl, 2, 4);

    let alphas_anl = Anl::<Event>::new()
        // i/o management
        .with_input_directory("D:\\tamu_data\\exp\\si28_c_35\\rkyv")
        .with_output_directory("K:\\tamu_data\\exp\\si28_c_35\\anl")
        // filter
        .with_filter(filter)
        .with_filtered_output_size(2_000_000)
        // modules
        .with_mixed_module(relative_energy_mixed)
        // settings
        .use_existing_filtered(true)
        .with_update_interval(args.update_interval)
        .with_max_real_events(args.events);

    let alphas_anl = assign_mixer(&args, alphas_anl);

    let alphas_anl = if args.skip_real {
        alphas_anl
    } else {
        alphas_anl.with_real_module(relative_energy)
    };

    let mut alphas_anl = match args.mixed {
        None => alphas_anl.with_max_mixed_events(MixedEventMaximum::Factor(num_integer::binomial(
            args.num_alphas,
            args.num_alphas_anl,
        ) as f64)),
        Some(value) => alphas_anl.with_max_mixed_events(MixedEventMaximum::Absolute(value)),
    };

    // run the analysis
    alphas_anl.run();

    println!("Analysis took {} seconds", start.elapsed().as_secs_f32());
}

fn assign_mixer<'a>(args: &'a Args, anl: Anl<'a, Event>) -> Anl<'a, Event> {
    match args.mixer.as_str() {
        "std" => {
            let mut mixer = GeneralParticleMixer::new().with_particles(args.num_alphas_anl, 2, 4);
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
