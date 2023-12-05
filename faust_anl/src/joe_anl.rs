use crate::general_particle_selection::ParticleSelectionRule;
use crate::relative_energy;
use crate::source::Source;
use cxx::{SharedPtr, UniquePtr};
use faust::event::Event;
use nukers::anl_module::AnlModule;
use roost::{
    file::RtFile,
    tree::{RtBranch, RtTree},
};

pub struct JoeAnlShift8Alpha {
    name: String,
    particle_rules: Vec<ParticleSelectionRule>,

    out_file: Option<SharedPtr<RtFile>>,

    tree: SharedPtr<RtTree>,
    alphas_8_erel: UniquePtr<RtBranch>,
    alphas_8_erel_shifted: UniquePtr<RtBranch>,
    alphas_7_erel: UniquePtr<RtBranch>,
}

impl JoeAnlShift8Alpha {
    pub fn new(name: &str) -> Self {
        let tree = roost::tree::new_tree("T".into(), "T".into());
        Self {
            name: name.into(),
            alphas_8_erel: tree.make_branch("erel_8".into()),
            alphas_7_erel: tree.make_branch("erel_7".into()),
            alphas_8_erel_shifted: tree.make_branch("erel_8_sft".into()),
            tree,
            particle_rules: vec![ParticleSelectionRule::new(8, 2, 4)],
            out_file: None,
        }
    }
    pub fn strip_target_like_alpha(source: &mut Source) {
        let idx = (0..source.mult()).min_by(|&idx1, &idx2| {
            let p1 = source.particles()[idx1].momentum_MeV_per_c();
            let p2 = source.particles()[idx2].momentum_MeV_per_c();
            p1.z().total_cmp(&p2.z())
        });
        source.swap_remove_particle(idx.unwrap());
    }
}

impl AnlModule<Event> for JoeAnlShift8Alpha {
    fn name(&self) -> String {
        self.name.clone()
    }

    fn initialize(&mut self, _output_directory: &std::path::Path) {
        let file_name: String = _output_directory
            .join(self.name.clone())
            .to_str()
            .unwrap()
            .into();

        self.out_file = Some(roost::file::create(file_name));
    }

    fn filter_event(&mut self, _event: &Event, _idx: usize) -> bool {
        relative_energy::event_passes_rules(_event, &self.particle_rules)
    }

    fn analyze_event(&mut self, event: &Event, _idx: usize) {
        let mut source = relative_energy::construct_source(event, &self.particle_rules);

        let erel_8 = source.relative_kinetic_energy_MeV();

        JoeAnlShift8Alpha::strip_target_like_alpha(&mut source);

        let erel_7 = source.relative_kinetic_energy_MeV();

        self.alphas_7_erel.pin_mut().set_value(erel_7);
        self.alphas_8_erel.pin_mut().set_value(erel_8);
        self.alphas_8_erel_shifted
            .pin_mut()
            .set_value(erel_8 - 37.4);

        self.out_file.as_ref().unwrap().cd();
        self.tree.fill();
    }

    fn finalize(&mut self) {
        self.out_file.as_ref().unwrap().cd();
        self.tree.write();
    }
}
