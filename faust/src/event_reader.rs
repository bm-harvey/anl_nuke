use crate::event::Event;

use std::fs::File;
use std::io::{BufRead, BufReader};

use crate::nuclear_masses::NUCLEAR_DB;

pub struct Reader {
    buffer_reader: BufReader<File>,
}

impl Reader {
    pub fn builder() -> ReaderBuilder {
        ReaderBuilder::default()
    }

    pub fn next_event(&mut self) -> Option<Event> {
        if self
            .buffer_reader
            .fill_buf()
            .expect("Buffer failed to read data")
            .is_empty()
        {
            return None;
        }

        // actually reading data
        let mut line = String::new();
        self.buffer_reader.read_line(&mut line).unwrap();

        let multiplicity = line
            .trim()
            .parse::<usize>()
            .expect("Failed to parse a multiplicity");

        let mut event = Event::default();

        for _ in 0..multiplicity {
            line.clear();
            self.buffer_reader.read_line(&mut line).unwrap();
            let words = line.split_whitespace().collect::<Vec<_>>();

            let charge: usize = words[0].parse().unwrap();
            let mass: usize = words[1].parse().unwrap();
            let px: f64 = words[2].parse::<f64>().unwrap()
                * (NUCLEAR_DB.nuclear_mass_MeV_per_c2(charge, mass) / (mass as f64) / 931.5).sqrt();
            let py: f64 = words[3].parse::<f64>().unwrap()
                * (NUCLEAR_DB.nuclear_mass_MeV_per_c2(charge, mass) / (mass as f64) / 931.5).sqrt();
            let pz: f64 = words[4].parse::<f64>().unwrap()
                * (NUCLEAR_DB.nuclear_mass_MeV_per_c2(charge, mass) / (mass as f64) / 931.5).sqrt();
            let detector: usize = words[5].parse().unwrap();

            // The following is converting some classical approximations back to the original values then converting to the values used in the `Event`
            // struct accounting for special relativity.
            // classical energy that would have been stored in the original file
            //let kinetic_energy =
            //(px.powi(2) + py.powi(2) + pz.powi(2)) / 2. / 931.5 / (mass as f64);
            //let nuclear_mass = NUCLEAR_DB.nuclear_mass_MeV_per_c2(charge, mass);
            //let lorentz_factor = kinetic_energy / nuclear_mass + 1.;
            //let p_mag = (lorentz_factor.powi(2) - 1.).sqrt() * nuclear_mass;
            //let momentum = PhysVec::from_cartesian(px, py, pz).as_normalized_to(p_mag);

            event.add_particle_by_properties(charge, mass, px, py, pz, detector);
        }
        Some(event)
    }
}

#[derive(Default)]
pub struct ReaderBuilder {
    input_file: String,
    capacity: usize,
}

impl ReaderBuilder {
    pub fn build(self) -> Reader {
        Reader {
            buffer_reader: BufReader::with_capacity(8_096 * 4, File::open(self.input_file).expect("unable to open file")),
        }
    }

    pub fn with_input_file(mut self, input_file: &str) -> Self {
        self.input_file = input_file.into();
        self
    }

    pub fn with_capacity(mut self, capacity: usize) -> Self {
        self.capacity = capacity;
        self
    }
}
