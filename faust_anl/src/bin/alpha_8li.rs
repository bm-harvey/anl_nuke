use faust::event::Event;
use faust_anl::relative_energy::RelativeEnergy;
use nukers::anl_module::Anl;
use nukers::anl_module::EventGenerator;

use faust_anl::GeneralParticleFilter;
use faust_anl::GeneralParticleMixer;
use faust_anl::MatchingPattern;

fn main() {
    let filter = GeneralParticleFilter::new(MatchingPattern::Standard)
        .with_particles(1, 3, 8)
        .with_particles(1, 2, 4);

    let mixer = GeneralParticleMixer::new()
        .with_particles(1, 3, 8)
        .with_particles(1, 2, 4);

    let generator = EventGenerator::Mixer(Box::new(mixer));

    let e_rel_mod = RelativeEnergy::new()
        .with_particles(1, 3, 8)
        .with_particles(1, 2, 4);

    let _e_rel_mod_mixed = RelativeEnergy::new()
        .with_particles(1, 3, 8)
        .with_particles(1, 2, 4);

    Anl::<Event>::new()
        // i/o management
        .with_input_directory("D:\\tamu_data\\exp\\si28_c_35\\rkyv")
        .with_output_directory("K:\\tamu_data\\exp\\si28_c_35\\anl")
        // filter
        .with_filter(filter)
        .with_filtered_output_size(2_000_000)
        // event mixer
        .with_event_generator(generator)
        // modules
        .with_real_module(e_rel_mod)
        // settings
        .use_existing_filtered(true)
        .with_update_interval(10_000)
        .run();
    let filter = GeneralParticleFilter::new(MatchingPattern::Strict)
        .with_particles(1, 3, 8)
        .with_particles(1, 2, 4);

    let mixer = GeneralParticleMixer::new()
        .with_particles(1, 3, 8)
        .with_particles(1, 2, 4);

    let e_rel_mod = RelativeEnergy::new()
        .with_particles(1, 3, 8)
        .with_particles(1, 2, 4);

    let _e_rel_mod_mixed = RelativeEnergy::new()
        .with_particles(1, 3, 8)
        .with_particles(1, 2, 4);

    Anl::<Event>::new()
        // i/o management
        .with_input_directory("D:\\tamu_data\\exp\\si28_c_35\\rkyv")
        .with_output_directory("K:\\tamu_data\\exp\\si28_c_35\\anl")
        // filter
        .with_filter(filter)
        .with_filtered_output_size(2_000_000)
        // event mixer
        .with_event_generator(EventGenerator::Mixer(Box::new(mixer)))
        // modules
        .with_real_module(e_rel_mod)
        // settings
        .use_existing_filtered(true)
        .with_update_interval(10_000)
        .run();

    let filter = GeneralParticleFilter::new(MatchingPattern::Minimum)
        .with_particles(1, 3, 8)
        .with_particles(1, 2, 4);

    let mixer = GeneralParticleMixer::new()
        .with_particles(1, 3, 8)
        .with_particles(1, 2, 4);

    let e_rel_mod = RelativeEnergy::new()
        .with_particles(1, 3, 8)
        .with_particles(1, 2, 4);

    let _e_rel_mod_mixed = RelativeEnergy::new()
        .with_particles(1, 3, 8)
        .with_particles(1, 2, 4);

    Anl::<Event>::new()
        // i/o management
        .with_input_directory("D:\\tamu_data\\exp\\si28_c_35\\rkyv")
        .with_output_directory("K:\\tamu_data\\exp\\si28_c_35\\anl")
        // filter
        .with_filter(filter)
        .with_filtered_output_size(2_000_000)
        // event mixer
        .with_event_generator(EventGenerator::Mixer(Box::new(mixer)))
        // modules
        .with_real_module(e_rel_mod)
        // settings
        .use_existing_filtered(true)
        .with_update_interval(10_000)
        .run();
}
