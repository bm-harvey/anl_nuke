use amd::event::Event;
use clap::Parser;
use n_z_equilibration::n_z_equilibration::NZEquilibration;
use n_z_equilibration::basic_stats::BasicStats;
use nukers::anl_module::Anl;

///Analysis scripte for N/Z Equilibration
#[derive(Parser, Debug)]
#[command(version)]
struct Args {
    ///Maximum of events to analyze
    #[arg(short = 'e', long, default_value_t=usize::MAX)]
    max_events: usize,
}

fn main() {
    let start = std::time::Instant::now();

    let args = Args::parse();

    let _nze_module = NZEquilibration::new()
        .with_combined_nucleus_min_mass(100)
        .with_plf_min_charge(15)
        .with_tlf_min_charge(15)
        .with_heavy_frag_min_charge(11)
        .with_heavy_frag_max_charge(35)
        .with_light_frag_min_charge(4)
        .with_light_frag_max_charge(20)
        .with_time_buffer(70.)
        .track_final_fragments_to_end_of_event(false)
        .boxed();

    let basic_stats = BasicStats::default();
    
    Anl::<Event>::new()
        .with_input_directory("../data/amd/zn70_zn70_35_g_rkyv/")
        .with_output_directory("./post_analysis/data/")
        //.with_module(nze_module)
        .with_real_module(basic_stats)
        .with_max_real_events(args.max_events)
        .run();

    println!("Analysis took {} s", start.elapsed().as_secs_f32());
}
