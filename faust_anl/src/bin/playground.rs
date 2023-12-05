use faust::event::Event;
use faust_anl::relative_energy::RelativeEnergy;
use nukers::anl_module::Anl;
use nukers::anl_module::MixedEventMaximum::{Absolute, Factor};

//use faust_anl::GeneralParticleFilter;
//use faust_anl::GeneralParticleMixer;
//use faust_anl::MatchingPattern;
use faust_anl::*;
//use nukers::anl_module::MixedEventMaximum::{Absolute, Factor};
//use faust_anl::RandomizedPhiMixer;
fn main() {
    //let mult = 3;
    //let z = 2;
    //let a = 4;

    let filter = GeneralParticleFilter::new(MatchingPattern::Minimum)
        .with_particles(2, 2, 4)
        .with_particles(1, 2, 3);

    //let mut mixer = GeneralParticleMixer::new()
    //.with_particles(2, 1, 2);
    //mixer.set_faust_filter();

    //let mut mixer = MoveGaussianMixer::new();
    let mut mixer = RandomizeLabPhiAngles::new();
    mixer.set_faust_filter();

    //let e_rel_mod = RelativeEnergy::new().with_particles(2, 1, 2);

    let e_rel_mod_rndm_phi = RelativeEnergyAfterScrambling::new()
        .with_particles(2, 2, 4)
        .with_particles(1, 2, 3);

    //let e_rel_mod_mixed = RelativeEnergy::new()
    //.with_particles(2, 1, 2);

    Anl::<Event>::new()
        .with_input_directory("D:\\tamu_data\\exp\\si28_c_35\\rkyv")
        .with_output_directory("K:\\tamu_data\\exp\\si28_c_35\\anl")
        .with_filter(filter)
        .with_filtered_output_size(10_000_000)
        //.with_mixer(mixer)
        //.with_real_module(e_rel_md)
        .with_real_module(e_rel_mod_rndm_phi)
        //.with_mixed_module(e_rel_mod_mixed)
        //.with_max_mixed_events(Absolute(1_000_000))
        .with_max_real_events(10_000_000)
        .use_existing_filtered(true)
        .with_update_interval(100_000)
        .run();
}
