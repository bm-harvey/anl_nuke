use amd::event::Event;
use std::path::Path;

use nukers::anl_module::*;

#[derive(Default)]
pub struct BasicStats {
    avg_impact_param: f32,
    avg_mult: f32,
    num_events: usize,
}

impl BasicStats {
    pub fn boxed(self) -> Box<Self> {
        Box::new(self)
    }
}

impl AnlModule<Event> for BasicStats {
    fn name(&self) -> String {
        "basic_stats".into()
    }

    fn initialize(&mut self, _input_directory: &Path) {
        self.avg_impact_param = 0.;
    }

    fn analyze_event(&mut self, _event: &Event, _event_idx: usize) {
        self.avg_impact_param += _event.impact_parameter();
        self.avg_mult += _event.last_time_step().multiplicity() as f32;
        self.num_events += 1;
    }

    fn finalize(&mut self) {
        self.avg_impact_param /= self.num_events as f32;
        self.avg_mult /= self.num_events as f32;
    }

    fn generate_output(&mut self, _output_directory: &Path) {
        println!("Avg. Impact Parameter : {} fm", self.avg_impact_param);
        println!("Avg. Multiplicity : {}", self.avg_mult);
    }
}
