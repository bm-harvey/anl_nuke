use faust::event::Event;
use faust_anl::RandomizeLabPhiAngles;
use faust_anl::relative_energy::RelativeEnergy;
use nukers::anl_module::Anl;

use faust_anl::GeneralParticleFilter;
use faust_anl::GeneralParticleMixer;
use faust_anl::MatchingPattern;
//use nukers::anl_module::MixedEventMaximum;

fn main() {
    let filter = GeneralParticleFilter::new(MatchingPattern::Strict).with_particles(2, 1, 2);

    let mut mixer = GeneralParticleMixer::new().with_particles(2, 1, 2);
    mixer.set_faust_filter();

    let e_rel_mod = RelativeEnergy::new().with_particles(2, 1, 2);

    let e_rel_mod_mixed = RelativeEnergy::new().with_particles(2, 1, 2);

    Anl::<Event>::new()
        .with_input_directory("D:\\tamu_data\\exp\\si28_c_35\\rkyv")
        .with_output_directory("K:\\tamu_data\\exp\\si28_c_35\\anl")
        .with_filter(filter)
        .with_filtered_output_size(2_000_000)
        .with_event_generator(nukers::anl_module::EventGenerator::Mixer(Box::new(mixer)))
        .with_mixed_module(e_rel_mod_mixed)
        .with_real_module(e_rel_mod)
        .use_existing_filtered(true)
        .with_update_interval(100_000)
        .run();

    let filter = GeneralParticleFilter::new(MatchingPattern::Strict).with_particles(2, 1, 2);

    let mut mixer = RandomizeLabPhiAngles::new();
    mixer.set_faust_filter();

    let e_rel_mod_mixed = RelativeEnergy::new().with_particles(2, 1, 2);

    Anl::<Event>::new()
        .with_input_directory("D:\\tamu_data\\exp\\si28_c_35\\rkyv")
        .with_output_directory("K:\\tamu_data\\exp\\si28_c_35\\anl")
        .with_filter(filter)
        .with_filtered_output_size(2_000_000)
        .with_event_generator(nukers::anl_module::EventGenerator::Scrambler(Box::new(mixer)))
        .with_mixed_module(e_rel_mod_mixed)
        .use_existing_filtered(true)
        .with_update_interval(100_000)
        .run();

    let filter = GeneralParticleFilter::new(MatchingPattern::Minimum).with_particles(2, 1, 2);

    let mut mixer = GeneralParticleMixer::new().with_particles(2, 1, 2);
    mixer.set_faust_filter();

    let e_rel_mod = RelativeEnergy::new().with_particles(2, 1, 2);

    let e_rel_mod_mixed = RelativeEnergy::new().with_particles(2, 1, 2);

    Anl::<Event>::new()
        .with_input_directory("D:\\tamu_data\\exp\\si28_c_35\\rkyv")
        .with_output_directory("K:\\tamu_data\\exp\\si28_c_35\\anl")
        .with_filter(filter)
        .with_filtered_output_size(2_000_000)
        .with_event_generator(nukers::anl_module::EventGenerator::Mixer(Box::new(mixer)))
        .with_mixed_module(e_rel_mod_mixed)
        .with_real_module(e_rel_mod)
        .use_existing_filtered(true)
        .with_update_interval(100_000)
        .run();

    let filter = GeneralParticleFilter::new(MatchingPattern::Minimum).with_particles(2, 1, 2);

    let mut mixer = RandomizeLabPhiAngles::new();
    mixer.set_faust_filter();

    let e_rel_mod_mixed = RelativeEnergy::new().with_particles(2, 1, 2);

    Anl::<Event>::new()
        .with_input_directory("D:\\tamu_data\\exp\\si28_c_35\\rkyv")
        .with_output_directory("K:\\tamu_data\\exp\\si28_c_35\\anl")
        .with_filter(filter)
        .with_filtered_output_size(2_000_000)
        .with_event_generator(nukers::anl_module::EventGenerator::Scrambler(Box::new(mixer)))
        .with_mixed_module(e_rel_mod_mixed)
        .use_existing_filtered(true)
        .with_update_interval(100_000)
        .run();
    

}
