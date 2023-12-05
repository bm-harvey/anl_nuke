use crate::event::*;
use crate::phys::PhysVec;

use std::{
    fs::File,
    io::{BufRead, BufReader},
};

pub struct EventReader {
    timestep_readers: Vec<TimeStepReader>,
    file_patterns: Vec<String>,
    expected_number_of_timesteps: usize,
    num_nucleons: usize,
}

impl EventReader {
    pub fn new(
        data_directory: &str,
        expected_number_of_timesteps: usize,
        num_nucleons: usize,
    ) -> Result<Self, &'static str> {
        let mut search_pattern = String::from(data_directory);

        search_pattern.push_str("*_T000.dat");
        let mut file_patterns = Vec::new();
        for file_pattern in glob::glob(&search_pattern).unwrap() {
            match file_pattern {
                Ok(path) => {
                    let path = String::from(path.to_str().unwrap());
                    let stripped_path = path.strip_suffix("_T000.dat").unwrap();
                    file_patterns.push(stripped_path.into());
                }

                Err(_) => {
                    return Err("error in glob");
                }
            }
        }

        if file_patterns.is_empty() {
            return Err("no files matching pattern found");
        }

        Ok(EventReader {
            // timestep_readers will get generated on an as-needed basis
            timestep_readers: Vec::new(),
            file_patterns,
            expected_number_of_timesteps,
            num_nucleons,
        })
    }

    /// returns Some(Event) if there is an event to be constucted still
    /// returns None if there are no more events represented in the files
    pub fn next_event(&mut self) -> Option<Event> {
        // look to see if at the end of file

        if self.timestep_readers.is_empty() {
            let result = self.establish_next_readers();
            result?;
        }

        let mut time_steps: Vec<TimeStep> = Vec::new();

        for reader in self.timestep_readers.iter_mut() {
            let time_step = reader.next_time_step();

            match time_step {
                None => {
                    self.timestep_readers.clear();
                    let result = self.establish_next_readers();
                    match result {
                        None => {
                            return None;
                        }
                        Some(_) => {
                            return self.next_event();
                        }
                    }
                }
                Some(time_step) => {
                    time_steps.push(time_step);
                }
            }
        }

        // need to get the impact parameter from the first step.
        let projectile = time_steps[0].fragment(0);
        let target = time_steps[0].fragment(1);
        let impact_parameter = projectile.r_vec().transverse() + target.r_vec().transverse();
        Some(Event::new(impact_parameter, time_steps))
    }

    pub fn establish_next_readers(&mut self) -> Option<()> {
        // need to create new buffer readers based on another set of runs
        while !self.file_patterns.is_empty() && self.timestep_readers.is_empty() {
            let search_pattern = self.file_patterns.last();
            search_pattern?;

            let mut search_pattern = search_pattern.unwrap().clone();
            self.file_patterns.pop();
            search_pattern.push('*');

            println!("opening file pattern : {}", search_pattern);

            for path in glob::glob(&search_pattern).unwrap().flatten() {
                self.timestep_readers.push(TimeStepReader::new(
                    path.to_str().unwrap(),
                    self.num_nucleons,
                ));
            }

            if self.timestep_readers.len() != self.expected_number_of_timesteps {
                self.timestep_readers.clear();
                println!("did not find the correct number of files for the requested pattern");
            } else {
                self.timestep_readers
                    .sort_by(|ts1, ts2| ts1.time().total_cmp(&ts2.time()));
                return Some(());
            }
        }
        if self.timestep_readers.is_empty() {
            None
        } else {
            Some(())
        }
    }

    pub fn print_file_patterns(&self) {
        for pattern in self.file_patterns.iter() {
            println!("{}", pattern);
        }
    }
}

pub struct TimeStepReader {
    /// The buffer reader used to read in the large files line by line
    buffer_reader: BufReader<File>,
    /// Time that the timestep reader is responsible for reading
    time: f32,
    /// Number of nucleons in the simulation... needs to be passed because AMD output is dumb.
    num_nucleons: usize,
}

impl TimeStepReader {
    pub fn new(file_path: &str, num_nucleons: usize) -> Self {
        // the file_path should look like ..._TXXX.dat
        let time_start = file_path.len() - 7;
        let time_end = file_path.len() - 4;
        let time_str: String = file_path[time_start..time_end].to_owned();
        let time = time_str.parse::<f32>().expect("failed to parse time");

        let file = File::open(file_path).unwrap();
        let buffer_reader = BufReader::with_capacity(100_000, file);
        //let buffer_reader = BufReader::new(file);

        TimeStepReader {
            buffer_reader,
            time,
            num_nucleons,
        }
    }

    fn next_time_step(&mut self) -> Option<TimeStep> {
        let mut line = String::new();

        self.buffer_reader.read_line(&mut line).unwrap();

        // early return None for the End of File tag
        if line.contains("0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0") {
            return None;
        }

        let mut accounted_for_nucleons = 0;

        let mut time_step = TimeStep::new(self.time);
        while accounted_for_nucleons < self.num_nucleons {
            if accounted_for_nucleons != 0 {
                line.clear();
                self.buffer_reader.read_line(&mut line).unwrap();
            }
            let words = line.split_whitespace().collect::<Vec<_>>();
            let fragment_charge = words[0]
                .parse::<u8>()
                .expect("failed to parse a fragment charge");
            let fragment_mass = fragment_charge
                + words[1]
                    .parse::<u8>()
                    .expect("failed to parse a fragment nucleons");
            let fragment_px = words[2]
                .parse::<f32>()
                .expect("failed to parse a fragment px")
                * (fragment_mass as f32);
            let fragment_py = words[3]
                .parse::<f32>()
                .expect("failed to parse a fragment py")
                * (fragment_mass as f32);

            accounted_for_nucleons += fragment_mass as usize;

            line.clear();
            self.buffer_reader.read_line(&mut line).unwrap();
            let words = line.split_whitespace().collect::<Vec<_>>();
            let fragment_pz = words[0]
                .parse::<f32>()
                .expect("failed to parse a fragment pz")
                * (fragment_mass as f32);

            let fragment_rx = words[1]
                .parse::<f32>()
                .expect("failed to read a fragment rx");

            let fragment_ry = words[2]
                .parse::<f32>()
                .expect("failed to read a fragment ry");

            line.clear();
            self.buffer_reader.read_line(&mut line).unwrap();
            let words = line.split_whitespace().collect::<Vec<_>>();
            let fragment_rz = words[0]
                .parse::<f32>()
                .expect("failed to read a fragment rz");
            let fragment_internal_energy = words[1]
                .parse::<f32>()
                .expect("failed to read a fragment internal energy");
            let fragment_total_l = words[2]
                .parse::<f32>()
                .expect("failed to read a fragment L");

            line.clear();
            self.buffer_reader.read_line(&mut line).unwrap();
            let words = line.split_whitespace().collect::<Vec<_>>();
            let fragment_jx = words[0]
                .parse::<f32>()
                .expect("failed to read a fragment Jx");
            let fragment_jy = words[1]
                .parse::<f32>()
                .expect("failed to read a fragment Jy");

            let fragment_jz = words[2]
                .parse::<f32>()
                .expect("failed to read a fragment Jz");
            line.clear();
            self.buffer_reader.read_line(&mut line).unwrap();
            let words = line.split_whitespace().collect::<Vec<_>>();
            //let rmsf = words[0].parse::<f32>().expect("failed to read a fragment rmsf");
            let fragment_max_density = words[1].parse::<f32>().expect("failed to read max density");

            let fragment_p_vec = PhysVec::new([fragment_px, fragment_py, fragment_pz]);
            let fragment_r_vec = PhysVec::new([fragment_rx, fragment_ry, fragment_rz]);
            let fragment_j_vec = PhysVec::new([fragment_jx, fragment_jy, fragment_jz]);

            let mut fragment = Fragment::new(
                fragment_charge,
                fragment_mass,
                fragment_p_vec,
                fragment_r_vec,
                fragment_j_vec,
                fragment_total_l,
                fragment_internal_energy,
                fragment_max_density,
                Vec::new(),
            );
            self.buffer_reader.read_line(&mut line).unwrap();

            for _ in 0..fragment_mass {
                line.clear();
                self.buffer_reader.read_line(&mut line).unwrap();
                let words = line.split_whitespace().collect::<Vec<_>>();

                let identifier = words[0]
                    .parse::<usize>()
                    .expect("failed to read nucleon identifier");

                let id = words[1]
                    .parse::<u8>()
                    .expect("failed to read nucleon identifier");

                let (proton, spin_up) = TimeStepReader::nucleon_type_from_int(identifier);

                let nucleon = Nucleon::new(id, proton, spin_up);

                self.buffer_reader.read_line(&mut line).unwrap();
                self.buffer_reader.read_line(&mut line).unwrap();

                fragment.add_nucleon(nucleon);
            }
            time_step.add_fragment(fragment);
        }
        Some(time_step)
    }

    fn nucleon_type_from_int(identifier: usize) -> (bool, bool) {
        let proton = identifier == 1 || identifier == 2;
        let spin_up = identifier == 1 || identifier == 3;

        (proton, spin_up)
    }

    fn time(&self) -> f32 {
        self.time
    }
}

//#[derive(Clone, Copy)]
//pub enum TimeStepReaderStatus {
//FileHeader,
//EventHeader,
//ImpactParameter,
//FragmentCount,
//FragmentHeader,
//NucleusInfo,
//NucleonInfo,
//}

//pub struct TimeStepClustReader {
///// The status of what line the reader should be on in the current file
//status: TimeStepReaderStatus,
///// The buffer reader used to read in the large files line by line
//buffer_reader: BufReader<File>,
///// The time in fm / c that the timestep represents
//time: f32,
//}

//impl TimeStepClustReader {
//pub fn new(file_path: &str) -> Self {
//let file = File::open(file_path).unwrap();
//let buffer_reader = BufReader::with_capacity(100_000, file);
////let buffer_reader = BufReader::new(file);

//let mut reader = TimeStepClustReader {
////file,
//status: TimeStepReaderStatus::FileHeader,
//buffer_reader,
//time: 0.,
//};

//reader.extract_time();
//reader
//}

//fn extract_time(&mut self) {
//match self.status {
//TimeStepReaderStatus::FileHeader => {
//let mut line = String::new();
//self.buffer_reader.read_line(&mut line).unwrap();
//self.time = line.trim().parse().unwrap();
//self.status = TimeStepReaderStatus::EventHeader;
//}
//_ => {
//println!("Time tried to be extracted at the wrong part of the file");
//}
//}
//}

//fn next_time_step(&mut self) -> (TimeStep, f32) {
//let mut line = String::new();

//self.buffer_reader.read_line(&mut line).unwrap();

//self.status = TimeStepReaderStatus::ImpactParameter;
//line.clear();
//self.buffer_reader
//.read_line(&mut line)
//.expect("Failed to read in impact parameter");
//let impact_parameter: f32 = line
//.trim()
//.parse()
//.expect("Failed to parse impact parameter");

//self.status = TimeStepReaderStatus::FragmentCount;
//line.clear();
//self.buffer_reader
//.read_line(&mut line)
//.expect("Failed to read in fragment count");
//let num_frags: usize = line
//.trim()
//.parse()
//.expect("Failed to parse number of fragments");

//let mut frag_list: Vec<Fragment> = Vec::new();
//for _ in 0..num_frags {
//self.status = TimeStepReaderStatus::NucleusInfo;
//line.clear();
//self.buffer_reader
//.read_line(&mut line)
//.expect("Failed to read in fragment info");
//let words = line.split_whitespace().collect::<Vec<_>>();
//let charge_num: u8 = words[0].parse::<u8>().unwrap() + words[1].parse::<u8>().unwrap();
//let neutron_num: u8 = words[2].parse::<u8>().unwrap() + words[3].parse::<u8>().unwrap();
//let mass_num = charge_num + neutron_num;

//line.clear();
//self.buffer_reader
//.read_line(&mut line)
//.expect("Failed to read in fragment info");

//let chars: Vec<char> = line.chars().collect();
//let vector_info: String = chars[4..].iter().collect();
//let vector_info = vector_info.split_whitespace().collect::<Vec<_>>();
//let frag_coordinate = PhysVec::new([
//vector_info[0]
//.parse::<f32>()
//.expect("Failed to parse a coordinate"),
//vector_info[2]
//.parse::<f32>()
//.expect("Failed to parse a coordinate"),
//vector_info[4]
//.parse::<f32>()
//.expect("Failed to parse a coordinate"),
//]);
//let factor = 931.5 / 0.4 * 2.;
//let frag_momentum = PhysVec::new([
//vector_info[1]
//.parse::<f32>()
//.expect("Failed to parse a momentum")
//* factor,
//vector_info[3]
//.parse::<f32>()
//.expect("Failed to parse a momentum")
//* factor,
//vector_info[5]
//.parse::<f32>()
//.expect("Failed to parse a momentum")
//* factor,
//]);

//let mut nucleons: Vec<Nucleon> = Vec::new();

//for _ in 0..mass_num {
//line.clear();
//self.buffer_reader
//.read_line(&mut line)
//.expect("Failed to read in fragment info");

//let chars: Vec<char> = line.chars().collect();
//let vector_info: String = chars[5..].iter().collect();
//let vector_info = vector_info.split_whitespace().collect::<Vec<_>>();
//let nucleon_coordinate = PhysVec::new([
//vector_info[0]
//.parse::<f32>()
//.expect("Failed to parse a coordinate"),
//vector_info[2]
//.parse::<f32>()
//.expect("Failed to parse a coordinate"),
//vector_info[4]
//.parse::<f32>()
//.expect("Failed to parse a coordinate"),
//]);
//let nucleon_momentum = PhysVec::new([
//vector_info[1]
//.parse::<f32>()
//.expect("Failed to parse a momentum"),
//vector_info[3]
//.parse::<f32>()
//.expect("Failed to parse a momentum"),
//vector_info[5]
//.parse::<f32>()
//.expect("Failed to parse a momentum"),
//]);
//let particle_info: usize = chars[1..2]
//.iter()
//.collect::<String>()
//.parse()
//.expect("Failed to read particle type info");

//let (proton, spin_up) = TimeStepClustReader::nucleon_type_from_int(particle_info);

//let id: u8 = vector_info.last().unwrap().parse().unwrap();
//let nucleon =
//Nucleon::new(id, nucleon_coordinate, nucleon_momentum, proton, spin_up);

//nucleons.push(nucleon);
//}
////FIXME
//frag_list.push(Fragment::new(
//charge_num,
//mass_num,
//frag_momentum,
//frag_coordinate,
//nucleons,
//));
//}

//let result = TimeStep::new(self.time, frag_list);

//(result, impact_parameter)
//}

//fn nucleon_type_from_int(identifier: usize) -> (bool, bool) {
//let proton = identifier == 1 || identifier == 2;
//let spin_up = identifier == 1 || identifier == 3;

//(proton, spin_up)
//}

//pub fn status(&self) -> TimeStepReaderStatus {
//self.status
//}
//}

//pub struct EventClustReader {
//timestep_readers: Vec<TimeStepClustReader>,
//file_patterns: Vec<String>,
//expected_number_of_timesteps: usize,
//}

//impl EventClustReader {
//pub fn new(
//data_directory: &str,
//expected_number_of_timesteps: usize,
//) -> Result<Self, &'static str> {
//let mut search_pattern = String::from(data_directory);

//search_pattern.push_str("*_T000_clust.dat");
//let mut file_patterns = Vec::new();
//for file_pattern in glob::glob(&search_pattern).unwrap() {
//match file_pattern {
//Ok(path) => {
//let path = String::from(path.to_str().unwrap());
//let stripped_path = path.strip_suffix("_T000_clust.dat").unwrap();
//file_patterns.push(stripped_path.into());
//}
//Err(e) => {
//println!("{}", e);
//return Err("error in glob");
//}
//}
//}

//if file_patterns.is_empty() {
//return Err("no files matching pattern found");
//}

//Ok(EventClustReader {
//// timestep_readers will get generated on an as-needed basis
//timestep_readers: Vec::new(),
//file_patterns,
//expected_number_of_timesteps,
//})
//}

///// returns Some(Event) if there is an event to be constucted still
///// returns None if there are no more events represented in the files
//pub fn next_event(&mut self) -> Option<Event> {
//// look to see if at the end of file
//let mut end_of_files = false;

//for reader in self.timestep_readers.iter_mut() {
//if reader.buffer_reader.fill_buf().unwrap().is_empty() {
//end_of_files = true;
//break;
//}
//}

//if end_of_files {
//self.timestep_readers.clear();
//}

//// return None if there are no buffers currently open, and no more to open
//if self.timestep_readers.is_empty() && self.file_patterns.is_empty() {
//return None;
//}

//// need to create new buffer readers based on another set of runs

//while !self.file_patterns.is_empty() && self.timestep_readers.is_empty() {
//let search_pattern = self.file_patterns.last();
//search_pattern?;

//let mut search_pattern = search_pattern.unwrap().clone();
//self.file_patterns.pop();
//search_pattern.push('*');

//println!("opening file pattern : {}", search_pattern);

//for path in glob::glob(&search_pattern).unwrap().flatten() {
//self.timestep_readers
//.push(TimeStepClustReader::new(path.to_str().unwrap()));
//}

//if self.timestep_readers.len() != self.expected_number_of_timesteps {
//self.timestep_readers.clear();
//println!("did not find the correct number of files for the requested pattern");
//}
//}

//let mut time_steps: Vec<TimeStep> = Vec::new();
//let mut impacts: Vec<f32> = Vec::new();

//for reader in self.timestep_readers.iter_mut() {
//let (step, impact) = reader.next_time_step();
//time_steps.push(step);
//impacts.push(impact);
//}

//let impact_parameter = impacts.get(0).unwrap();
//for b in impacts.iter() {
//if b != impact_parameter {
//println!("WARNING: IMPACT PARAMS BETWEEN TIMESTEPS DO NOT MATCH");
//}
//}

//Some(Event::new(*impact_parameter, time_steps))
//}

//pub fn print_file_patterns(&self) {
//for pattern in self.file_patterns.iter() {
//println!("{}", pattern);
//}
//}
//}
