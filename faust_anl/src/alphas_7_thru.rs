use crate::source::Source;
use faust::event::{Event, Particle};
use itertools::Itertools;
use nukers::anl_module::{EventFilter, EventMixer};
use nukers::data_set::ArchivedData;

// THROUGH HOYLE STATE
pub struct Si28AlphasThroughHoyle {
    name: String,
}

impl Si28AlphasThroughHoyle {
    const HOYLE_THRESHOLD_MEV: f64 = 0.6;
    pub fn new(name: &str) -> Self {
        Self { name: name.into() }
    }
}

impl Default for Si28AlphasThroughHoyle {
    fn default() -> Self {
        Si28AlphasThroughHoyle::new("7a_through_hoyle")
    }
}

impl EventFilter<Event> for Si28AlphasThroughHoyle {
    fn name(&self) -> String {
        self.name.clone()
    }

    fn filter_event(&self, event: &Event, _idx: usize) -> bool {
        if event.pid_mult(2, 4) != 7 {
            return false;
        }

        for combo_alphas_3 in event.iter_by_pid(2, 4).combinations(3) {
            let mut potential_hoyle = crate::source::Source::new();

            for alpha in combo_alphas_3 {
                potential_hoyle.add_particle(alpha);
            }

            if potential_hoyle.relative_kinetic_energy_MeV()
                < Si28AlphasThroughHoyle::HOYLE_THRESHOLD_MEV
            {
                return true;
            }
        }

        false
    }
}

/// Mix 7a particles, but require that 3 of them are hoyle state
pub struct Si28AlphasThroughHoyleMixer {
    name: String,
}

impl Default for Si28AlphasThroughHoyleMixer {
    fn default() -> Self {
        Self {
            name: "7a_thru_hoyle".into(),
        }
    }
}

impl EventMixer<Event> for Si28AlphasThroughHoyleMixer {
    fn name(&self) -> String {
        self.name.clone()
    }

    /// Grabs 3 alphas that are consistent with a hoyle state from one event, and a random alpha
    /// from 4 other independent events
    fn mix_events(&mut self, data: &ArchivedData<Event>, _idx: usize) -> Option<Event> {
        let mut result = Event::default();

        // Get the hoyle-alphas
        let real_event_for_hoyle = data.random_event();

        let mut hoyle = None;
        for combo in real_event_for_hoyle.iter_by_pid(2, 4).combinations(3) {
            let mut source = Source::new();
            for alpha in combo {
                source.add_particle(alpha);
            }
            if source.relative_kinetic_energy_MeV() < Si28AlphasThroughHoyle::HOYLE_THRESHOLD_MEV {
                hoyle = Some(source)
            }
        }

        hoyle.as_ref()?;

        for alpha in hoyle.unwrap().iter() {
            result.add_particle_by_clone(alpha);
        }

        // Get the 4 random alphas
        let real_events_for_alphas = (0..4).map(|_| data.random_event()).collect::<Vec<_>>();
        for event in real_events_for_alphas.iter() {
            let rndm_alpha = event.random_particle_by_pid(2, 4);
            rndm_alpha?;
            result.add_particle_by_clone(rndm_alpha.unwrap());
        }

        // Final check
        let detectors_are_unique = result
            .iter()
            .map(|particle| particle.detector())
            .all_unique();

        match detectors_are_unique {
            true => Some(result),
            false => None,
        }
    }
}

// THROUGH SINGLE 8Be g.s.
/// 7a -> at least one 8Be g.s.
pub struct Si28AlphasThroughBe8gs {
    name: String,
}

impl Si28AlphasThroughBe8gs {
    const BE8GS_THRESHOLD_MEV: f64 = 0.2;
    pub fn new(name: &str) -> Self {
        Self { name: name.into() }
    }
}

impl Default for Si28AlphasThroughBe8gs {
    fn default() -> Self {
        Si28AlphasThroughBe8gs::new("7a_through_be8gs")
    }
}

impl EventFilter<Event> for Si28AlphasThroughBe8gs {
    fn name(&self) -> String {
        self.name.clone()
    }

    fn filter_event(&self, event: &Event, _idx: usize) -> bool {
        if event.pid_mult(2, 4) != 7 {
            return false;
        }

        for combo_alphas_2 in event.iter_by_pid(2, 4).combinations(2) {
            let mut potential_be8gs = crate::source::Source::new();

            for alpha in combo_alphas_2 {
                potential_be8gs.add_particle(alpha);
            }

            if potential_be8gs.relative_kinetic_energy_MeV()
                < Si28AlphasThroughBe8gs::BE8GS_THRESHOLD_MEV
            {
                return true;
            }
        }

        false
    }
}

/// Mix 7a particles, but require that 2 of them are 8Be gs state
pub struct Si28AlphasThroughBe8gsMixer {
    name: String,
}

impl Default for Si28AlphasThroughBe8gsMixer {
    fn default() -> Self {
        Self {
            name: "7a_thru_be8gs".into(),
        }
    }
}

impl EventMixer<Event> for Si28AlphasThroughBe8gsMixer {
    fn name(&self) -> String {
        self.name.clone()
    }

    /// Grabs 2 alphas that are consistent with a 8Be gs state from one event, and a random alpha
    /// from 5 other independent events
    fn mix_events(&mut self, data: &ArchivedData<Event>, _idx: usize) -> Option<Event> {
        let mut result = Event::default();

        // Get the 8Be-alphas
        let real_event_for_8begs = data.random_event();

        let mut be8gs = None;
        for combo in real_event_for_8begs.iter_by_pid(2, 4).combinations(2) {
            let mut source = Source::new();
            for alpha in combo {
                source.add_particle(alpha);
            }
            if source.relative_kinetic_energy_MeV() < Si28AlphasThroughBe8gs::BE8GS_THRESHOLD_MEV {
                be8gs = Some(source)
            }
        }

        be8gs.as_ref()?;

        for alpha in be8gs.unwrap().iter() {
            result.add_particle_by_clone(alpha);
        }

        // Get the 5 random alphas
        let real_events_for_alphas = (0..5).map(|_| data.random_event()).collect::<Vec<_>>();
        for event in real_events_for_alphas.iter() {
            let rndm_alpha = event.random_particle_by_pid(2, 4);
            rndm_alpha?;
            result.add_particle_by_clone(rndm_alpha.unwrap());
        }

        // Final check
        let detectors_are_unique = result
            .iter()
            .map(|particle| particle.detector())
            .all_unique();

        match detectors_are_unique {
            true => Some(result),
            false => None,
        }
    }
}

// THROUGH TWO 8Be g.s.
/// 7a -> at least two 8Be g.s.
pub struct Si28AlphasThrough2Be8gs {
    name: String,
}

impl Si28AlphasThrough2Be8gs {
    const BE8GS_THRESHOLD_MEV: f64 = 0.2;
    pub fn new(name: &str) -> Self {
        Self { name: name.into() }
    }
}

impl Default for Si28AlphasThrough2Be8gs {
    fn default() -> Self {
        Si28AlphasThrough2Be8gs::new("7a_through_2be8gs")
    }
}

impl EventFilter<Event> for Si28AlphasThrough2Be8gs {
    fn name(&self) -> String {
        self.name.clone()
    }

    fn filter_event(&self, event: &Event, _idx: usize) -> bool {
        if event.pid_mult(2, 4) != 7 {
            return false;
        }

        let mut be8gs_list = Vec::<Source>::new();
        for combo_alphas_2 in event.iter_by_pid(2, 4).combinations(2) {
            let mut potential_be8gs = crate::source::Source::new();

            for alpha in combo_alphas_2 {
                potential_be8gs.add_particle(alpha);
            }

            if potential_be8gs.relative_kinetic_energy_MeV()
                < Si28AlphasThrough2Be8gs::BE8GS_THRESHOLD_MEV
            {
                be8gs_list.push(potential_be8gs.clone());
            }
        }

        if be8gs_list.len() < 2 {
            return false;
        }

        // this for loop checks if any of the pairs of Be8 in the system have a set of 4 unique
        // alphas. The if statement is some pointer and iterator magic, but all it's is checking is
        // that if one loops through all of the alphas of a pair of sources, are all of the alpkhas
        // unique (by pointer)
        for pair in be8gs_list.iter().combinations(2) {
            if pair
                .iter()
                .flat_map(|&source| source.iter())
                .map(|&p| p as *const Particle)
                .all_unique()
            {
                return true;
            }
        }

        false
    }
}

/// Mix 7a particles, but require that 2 of them are 8Be gs state
pub struct Si28AlphasThrough2Be8gsMixer {
    name: String,
}

impl Default for Si28AlphasThrough2Be8gsMixer {
    fn default() -> Self {
        Self {
            name: "7a_thru_2be8gs".into(),
        }
    }
}

impl EventMixer<Event> for Si28AlphasThrough2Be8gsMixer {
    fn name(&self) -> String {
        self.name.clone()
    }

    /// Grabs 2 alphas that are consistent with a 8Be gs state from one event (twice), and a random alpha
    /// from 3 other independent events
    fn mix_events(&mut self, data: &ArchivedData<Event>, _idx: usize) -> Option<Event> {
        let mut result = Event::default();

        // Get the 8Be-alphas
        for _ in 0..2 {
            let real_event_for_8begs = data.random_event();

            let mut be8gs = None;
            for combo in real_event_for_8begs.iter_by_pid(2, 4).combinations(2) {
                let mut source = Source::new();
                for alpha in combo {
                    source.add_particle(alpha);
                }
                if source.relative_kinetic_energy_MeV() < Si28AlphasThroughBe8gs::BE8GS_THRESHOLD_MEV {
                    be8gs = Some(source)
                }
            }

            be8gs.as_ref()?;

            for alpha in be8gs.unwrap().iter() {
                result.add_particle_by_clone(alpha);
            }
        }

        // Get the 3 random alphas
        let real_events_for_alphas = (0..3).map(|_| data.random_event()).collect::<Vec<_>>();
        for event in real_events_for_alphas.iter() {
            let rndm_alpha = event.random_particle_by_pid(2, 4);
            rndm_alpha?;
            result.add_particle_by_clone(rndm_alpha.unwrap());
        }

        // Final check
        let detectors_are_unique = result
            .iter()
            .map(|particle| particle.detector())
            .all_unique();

        match detectors_are_unique {
            true => Some(result),
            false => None,
        }
    }
}

// THROUGH No 8Be g.s.
pub struct Si28AlphasThrough0Be8gs {
    name: String,
}

impl Si28AlphasThrough0Be8gs {
    const BE8GS_THRESHOLD_MEV: f64 = 0.2;
    pub fn new(name: &str) -> Self {
        Self { name: name.into() }
    }
}

impl Default for Si28AlphasThrough0Be8gs {
    fn default() -> Self {
        Si28AlphasThrough0Be8gs::new("7a_not_through_be8gs")
    }
}

impl EventFilter<Event> for Si28AlphasThrough0Be8gs {
    fn name(&self) -> String {
        self.name.clone()
    }

    fn filter_event(&self, event: &Event, _idx: usize) -> bool {
        if event.pid_mult(2, 4) != 7 {
            return false;
        }

        for combo_alphas_2 in event.iter_by_pid(2, 4).combinations(2) {
            let mut potential_be8gs = crate::source::Source::new();

            for alpha in combo_alphas_2 {
                potential_be8gs.add_particle(alpha);
            }

            if potential_be8gs.relative_kinetic_energy_MeV()
                < Si28AlphasThrough0Be8gs::BE8GS_THRESHOLD_MEV
            {
                return false;
            }
        }

        true 
    }
}

