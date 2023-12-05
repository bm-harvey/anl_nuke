use faust::event::Particle;
use faust::nuclear_masses::NUCLEAR_DB;
use faust::phys_vec::PhysVec;
use nalgebra::Matrix3;

#[derive(Clone, Default)]
pub struct Source<'a> {
    particles: Vec<&'a Particle>,
}

impl<'a> Source<'a> {
    pub fn new() -> Self {
        Self {
            particles: Vec::new(),
        }
    }
    pub fn from(particles: &'a [Particle]) -> Self {
        let mut source = Source::new();
        for particle in particles {
            source.add_particle(particle);
        }
        source
    }

    pub fn clear(&mut self) -> &mut Self {
        self.particles.clear();
        self
    }

    pub fn num_particles(&self) -> usize {
        self.particles.len()
    }

    pub fn sort_by_momentum_magnitude(&mut self) -> &mut Self {
        self.particles
            .sort_by(|p1, p2| p1.p().mag().total_cmp(&p2.p().mag()));
        self
    }

    pub fn add_particle(&mut self, particle: &'a Particle) -> &mut Self {
        self.particles.push(particle);
        self
    }

    pub fn swap_remove_particle(&mut self, idx: usize) -> &mut Self {
        self.particles.swap_remove(idx);
        self
    }

    pub fn remove_particle(&mut self, idx: usize) -> &mut Self {
        self.particles.remove(idx);
        self
    }

    pub fn charge_number(&self) -> usize {
        self.particles
            .iter()
            .map(|particle| particle.proton_num() as usize)
            .sum()
    }

    pub fn mass_number(&self) -> usize {
        self.particles
            .iter()
            .map(|particle| particle.mass_num() as usize)
            .sum()
    }

    #[allow(non_snake_case)]
    pub fn children_mass_sum_MeV_per_c2(&self) -> f64 {
        self.particles
            .iter()
            .map(|particle| particle.mass_MeV_per_c2())
            .sum()
    }

    #[allow(non_snake_case)]
    pub fn mass_MeV_per_c2(&self) -> f64 {
        NUCLEAR_DB
            .nuclear_mass_MeV_per_c2(self.charge_number() as usize, self.mass_number() as usize)
    }

    #[allow(non_snake_case)]
    pub fn q_value_MeV(&self) -> f64 {
        self.mass_MeV_per_c2() - self.children_mass_sum_MeV_per_c2()
    }

    #[allow(non_snake_case)]
    pub fn excitation_energy_MeV(&self) -> f64 {
        self.relative_kinetic_energy_MeV() - self.q_value_MeV()
    }

    pub fn iter(&self) -> std::slice::Iter<'_, &Particle> {
        self.particles.iter()
    }

    pub fn particles(&self) -> &Vec<&'a Particle> {
        &self.particles
    }

    #[allow(non_snake_case)]
    pub fn velocity_c(&self) -> PhysVec {
        self.momentum_MeV_per_c()
            .as_scaled_by(1. / self.children_mass_sum_MeV_per_c2())
    }

    #[allow(non_snake_case)]
    pub fn momentum_MeV_per_c(&self) -> PhysVec {
        let mut source_momentum = PhysVec::default();
        self.particles
            .iter()
            .for_each(|particle| source_momentum += particle.classical_momentum_MeV_per_c());

        source_momentum
    }

    #[allow(non_snake_case)]
    pub fn source_kinetic_energy_MeV(&self) -> f64 {
        0.5 * self.mass_MeV_per_c2() * self.velocity_c().mag_sqr()
    }

    #[allow(non_snake_case)]
    pub fn relative_kinetic_energies_MeV(&self) -> Vec<f64> {
        let velocity = self.velocity_c();
        self.particles
            .iter()
            .map(|particle| {
                0.5 * particle.mass_MeV_per_c2()
                    * (&particle.classical_velocity_c() - &velocity).mag_sqr()
            })
            .collect()
    }

    #[allow(non_snake_case)]
    pub fn relative_momenta_MeV_per_c(&self) -> Vec<f64> {
        let velocity = self.velocity_c();
        self.particles
            .iter()
            .map(|particle| particle.mass_MeV_per_c2() * (&particle.velocity_c() - &velocity).mag())
            .collect()
    }

    #[allow(non_snake_case)]
    pub fn relative_kinetic_energy_MeV(&self) -> f64 {
        let source_velocity = self.velocity_c();
        self.particles
            .iter()
            .map(|particle| {
                let relative_velocity = &particle.classical_velocity_c() - &source_velocity;
                0.5 * particle.mass_MeV_per_c2() * relative_velocity.mag_sqr()
            })
            .sum()
    }

    pub fn mult(&self) -> usize {
        self.particles.len()
    }

    pub fn shape(&self) -> (f64, f64) {
        let source_velocity = self.velocity_c();

        let mut momentum_tensor = Matrix3::<f64>::default();
        for particle in self.iter() {
            let p = (&particle.velocity_c() - &source_velocity)
                .as_scaled_by(particle.mass_MeV_per_c2());

            for i in 0..3 {
                for j in 0..3 {
                    momentum_tensor[(i, j)] += p.at(i) * p.at(j);
                }
            }
        }

        let eigen_values = momentum_tensor.eigenvalues();

        if eigen_values.is_none() {
            return (0., 0.);
        }

        let eigen_values = eigen_values.unwrap();

        let norm = eigen_values.iter().sum::<f64>();

        let mut normalized_eigen_values = eigen_values
            .iter()
            .map(|lambda| lambda / norm)
            .collect::<Vec<_>>();

        normalized_eigen_values.sort_by(|a, b| a.total_cmp(b));

        let lambda_1 = normalized_eigen_values[0];
        let lambda_2 = normalized_eigen_values[1];
        let lambda_3 = normalized_eigen_values[2];

        //let normalized_eigen_values;
        let sphericity = 1.5 * (1. - lambda_3);
        let coplanarity = 3_f64.sqrt() * 0.5 * (lambda_2 - lambda_1);

        (sphericity, coplanarity)
    }

    pub fn relativistic_source(&self) -> RelativisticSource {
        RelativisticSource::new(self)
    }
}

pub struct RelativisticSource {
    particles: Vec<Particle>,
    velocity: PhysVec,
}

impl RelativisticSource {
    pub fn new(source: &Source) -> Self {
        let energy = source
            .particles()
            .iter()
            .map(|&p| p.energy_MeV())
            .sum::<f64>();

        let source_velocity = source
            .particles()
            .iter()
            .fold(PhysVec::default(), |total, particle| {
                total + particle.momentum_MeV_per_c()
            })
            .as_scaled_by(1.0 / energy);

        let particles = source
            .particles()
            .iter()
            .map(|&p| {
                let mut boosted_particle = p.clone();
                boosted_particle.boost(&source_velocity);
                boosted_particle
            })
            .collect::<Vec<_>>();

        Self {
            particles,
            velocity: source_velocity,
        }
    }

    pub fn particles(&self) -> &Vec<Particle> {
        &self.particles
    }

    pub fn velocity_c(&self) -> &PhysVec {
        &self.velocity
    }

    pub fn mass_number(&self) -> usize {
        self.particles
            .iter()
            .map(|particle| particle.mass_num())
            .sum()
    }

    pub fn proton_number(&self) -> usize {
        self.particles
            .iter()
            .map(|particle| particle.proton_num())
            .sum()
    }

    #[allow(non_snake_case)]
    pub fn reconstructed_mass_MeV_per_c2(&self) -> f64 {

        NUCLEAR_DB.nuclear_mass_MeV_per_c2(self.proton_number(), self.mass_number())
    }

    #[allow(non_snake_case)]
    pub fn children_mass_MeV_per_c2(&self) -> f64 {
        self.particles.iter().map(|p| p.mass_MeV_per_c2()).sum()
    }

    pub fn lorentz_factor(&self) -> f64 {
        1. / (1. - self.velocity.mag_sqr()).sqrt()
    }

    #[allow(non_snake_case)]
    pub fn sum_of_energies_MeV(&self) -> f64 {
        self.particles
            .iter()
            .map(|particle| particle.energy_MeV())
            .sum()
    }

    #[allow(non_snake_case)]
    pub fn source_momentum_MeV_per_c(&self) -> PhysVec {
        self.velocity_c()
            .as_scaled_by(self.reconstructed_mass_MeV_per_c2() * self.lorentz_factor())
    }

    #[allow(non_snake_case)]
    pub fn relative_kinetic_energy_MeV(&self) -> f64 {

        self.particles
            .iter()
            .map(|particle| particle.kinetic_energy_MeV())
            .sum::<f64>()
    }

    #[allow(non_snake_case)]
    pub fn kinetic_energy_MeV(&self) -> f64 {
        (self.lorentz_factor() - 1.) * self.reconstructed_mass_MeV_per_c2()
    }

    #[allow(non_snake_case)]
    pub fn q_value_MeV_per_c2(&self) -> f64 {
        self.reconstructed_mass_MeV_per_c2() - self.children_mass_MeV_per_c2()
    }

    #[allow(non_snake_case)]
    pub fn excitation_energy_MeV(&self) -> f64 {
        self.relative_kinetic_energy_MeV() - self.q_value_MeV_per_c2()
    }
}
