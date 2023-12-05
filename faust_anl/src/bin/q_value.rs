use faust::nuclear_masses::NUCLEAR_DB;

fn main() {
    let mut input = String::new();
    let mut mass_entrance_channel = 0.;
    println!("Enter entrance channel nuclei.");
    loop {
        input.clear();
        std::io::stdin().read_line(&mut input).unwrap();
        if input.trim() == "" {
            break;
        }
        mass_entrance_channel += mass(&input);
    }

    println!("Enter exit channel nuclei.");
    let mut mass_exit_channel = 0.;
    loop {
        input.clear();
        std::io::stdin().read_line(&mut input).unwrap();
        if input.trim() == "" {
            break;
        }
        mass_exit_channel += mass(&input);
    }

    let q_value = mass_entrance_channel - mass_exit_channel;
    println!("Q-value: {} MeV", q_value);
}

fn mass(input_line: &str) -> f64 {
    let line = input_line.trim().split_whitespace().collect::<Vec<&str>>();
    let num_nuclei = line[0]
        .parse::<usize>()
        .expect("multilicity value could not be parsed");
    let charge_num = line[1]
        .parse::<usize>()
        .expect("proton number could not be parsed");
    let mass_num = line[2]
        .parse::<usize>()
        .expect("mass number could not be parsed");

    NUCLEAR_DB.nuclear_mass_MeV_per_c2(charge_num, mass_num) * num_nuclei as f64
}
