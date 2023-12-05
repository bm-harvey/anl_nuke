use faust_anl::GeneralParticleFilter;
use faust_anl::relative_energy::RelativeEnergy;
use faust_anl::MatchingPattern;

use faust_anl::alphas_7_thru::Si28AlphasThrough0Be8gs;
use faust_anl::alphas_7_thru::Si28AlphasThrough2Be8gs;
use faust_anl::alphas_7_thru::Si28AlphasThroughBe8gs;
use nukers::anl_module::Anl;

fn main() {
    //
    // All 7 alpha
    //
    let anl = RelativeEnergy::new().with_particles(7, 2, 4);
    let filter = GeneralParticleFilter::new(MatchingPattern::Standard);
    Anl::new()
        .with_input_directory("D:\\tamu_data\\exp\\si28_c_35\\rkyv")
        .with_output_directory("D:\\tamu_data\\exp\\si28_c_35\\anl")
        .with_filter(filter)
        .with_real_module(anl)
        .run();
    //
    // 7a -> not Be8gs
    //
    let anl = RelativeEnergy::new().with_particles(7, 2, 4);
    let filter = Si28AlphasThrough0Be8gs::default();
    Anl::new()
        .with_input_directory("D:\\tamu_data\\exp\\si28_c_35\\rkyv")
        .with_output_directory("D:\\tamu_data\\exp\\si28_c_35\\anl")
        .with_filter(filter)
        .with_real_module(anl)
        .run();

    //
    // 7a -> Be8gs
    //
    let anl = RelativeEnergy::new().with_particles(7, 2, 4);
    let filter = Si28AlphasThroughBe8gs::default();
    Anl::new()
        .with_input_directory("D:\\tamu_data\\exp\\si28_c_35\\rkyv")
        .with_output_directory("D:\\tamu_data\\exp\\si28_c_35\\anl")
        .with_filter(filter)
        .with_real_module(anl)
        .run();

    //
    // 7a -> 2 Be8gs
    //
    let anl = RelativeEnergy::new().with_particles(7, 2, 4);
    let filter = Si28AlphasThrough2Be8gs::default();
    Anl::new()
        .with_input_directory("D:\\tamu_data\\exp\\si28_c_35\\rkyv")
        .with_output_directory("D:\\tamu_data\\exp\\si28_c_35\\anl")
        .with_filter(filter)
        .with_real_module(anl)
        .run();
}
