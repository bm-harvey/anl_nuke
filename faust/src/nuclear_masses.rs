use std::collections::HashMap;

#[derive(Debug)]
pub struct NucleusData {
    _symbol: String,
    mass_excess: f64,
    atomic_mass: f64,
    binding_energy_per_a: f64,
}
impl NucleusData {
    pub fn new(
        symbol: String,
        mass_excess: f64,
        atomic_mass: f64,
        binding_energy_per_a: f64,
    ) -> NucleusData {
        NucleusData {
            _symbol: symbol,
            mass_excess,
            atomic_mass,
            binding_energy_per_a,
        }
    }
}

lazy_static::lazy_static! {
    pub static ref NUCLEAR_DB: NuclearDB = NuclearDB::new();
}

#[derive(Default)]
pub struct NuclearDB {
    data: HashMap<(usize, usize), NucleusData>,
}

impl NuclearDB {
    const AMU: f64 = 931.5;
    //const SPEED_OF_LIGHT: f64 = 299_792_458.0;


    #[allow(non_snake_case)]
    pub fn atomic_mass_amu(&self, charge: usize, mass_number: usize) -> f64 {
        let entry = self.data.get(&(charge, mass_number)).unwrap();
        entry.atomic_mass
    }
    
    #[allow(non_snake_case)]
    pub fn nuclear_mass_MeV_per_c2(&self, charge: usize, mass_number: usize) -> f64 {
        let mass = Self::AMU * mass_number as f64;
        let entry = self.data.get(&(charge, mass_number)).unwrap();
        mass + entry.mass_excess
    }

    #[allow(non_snake_case)]
    pub fn binding_energy_per_nucleon_MeV(&self, charge: usize, mass_number: usize) -> f64 {
        let entry = self.data.get(&(charge, mass_number)).unwrap();
        entry.binding_energy_per_a
    }

    pub fn new() -> NuclearDB {
        let raw_data = include_str!("nuclear_masses.dat"); 

        let mut data = NuclearDB {
            data: HashMap::new(),
        };

        let lines = raw_data.lines().skip(36);

        for line in lines {
            //let current_line = line;
            let chars = line.chars().collect::<Vec<char>>();

            let mass_num_str = chars[16..20].iter().collect::<String>();
            let mass_number = mass_num_str.trim().parse::<usize>().unwrap();

            let charge_num_str = chars[11..15].iter().collect::<String>();
            let charge = charge_num_str.trim().parse::<usize>().unwrap();

            let symbol = chars[20..22].iter().collect::<String>().trim().to_string();

            let mass_excess_str = chars[29..40].iter().collect::<String>().trim().to_string();
            let mass_excess_str = mass_excess_str
                .chars()
                .map(|c| if c == '#' { '.' } else { c })
                .collect::<String>();
            let mass_excess = mass_excess_str.parse::<f64>().unwrap_or(0.) / 1_000.;

            let binding_energy_per_a = chars[57..67].iter().collect::<String>().trim().to_string();
            let binding_energy_per_a = binding_energy_per_a
                .chars()
                .map(|c| if c == '#' { '.' } else { c })
                .collect::<String>();
            let binding_energy_per_a = binding_energy_per_a.parse::<f64>().unwrap_or(0.) / 1_000.;

            let atomic_mass_str = chars[106..123]
                .iter()
                .collect::<String>()
                //.trim()
                .to_string();
            let atomic_mass_str = atomic_mass_str
                .chars()
                .filter(|c| *c != ' ')
                .map(|c| if c == '#' { '.' } else { c })
                .collect::<String>();
            let atomic_mass = atomic_mass_str.parse::<f64>().unwrap_or(0.) / 1_000.;

            let entry = NucleusData::new(symbol, mass_excess, atomic_mass, binding_energy_per_a);

            data.data.insert((charge, mass_number), entry);
        }
        data
    }

    pub fn mass(&self) -> f64 {
        0.0
    }
}
