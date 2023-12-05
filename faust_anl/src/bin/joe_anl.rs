use faust_anl::MatchingPattern;
use faust_anl::{joe_anl, GeneralParticleFilter};
use nukers::anl_module::Anl;

fn main() {
    let anl = joe_anl::JoeAnlShift8Alpha::new("joes_8_alphas");

    Anl::new()
        .with_input_directory("D:\\tamu_data\\exp\\si28_c_35\\rkyv")
        .with_output_directory("D:\\tamu_data\\exp\\si28_c_35\\anl")
        .with_filter(GeneralParticleFilter::new(MatchingPattern::Standard).with_particles(8, 2, 4))
        .with_real_module(anl)
        .run();
}
