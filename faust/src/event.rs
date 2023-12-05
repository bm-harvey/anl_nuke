use crate::nuclear_masses::NUCLEAR_DB;
use crate::phys_vec::PhysVec;
use rand::seq::{IteratorRandom, SliceRandom};
use rkyv::{Archive, Deserialize, Serialize};
use std::slice::{Iter, IterMut};

#[derive(Clone, Default, Debug, Archive, Serialize, Deserialize)]
#[archive(check_bytes)]
/// Represents an nuclear event as measured by FAUST. Contains a list of the particles measured.
/// Most of the analysis of the events is done by analyzing the particles in the event.
pub struct Event {
    particles: Vec<Particle>,
}

impl Event {
    /// Shared reference to the particle list.
    pub fn particles(&self) -> &Vec<Particle> {
        &self.particles
    }

    /// Iterate over the particles in the event.
    pub fn iter(&self) -> Iter<'_, Particle> {
        self.particles.iter()
    }

    /// Iterate mutably over the particles in the event.
    pub fn iter_mut(&mut self) -> IterMut<'_, Particle> {
        self.particles.iter_mut()
    }

    /// Returns the number of particles in the event
    pub fn mult(&self) -> usize {
        self.particles.len()
    }

    /// Iterate over the particles in the event with the specified Z and A values.
    pub fn iter_by_pid(&self, z: usize, a: usize) -> impl Iterator<Item = &Particle> {
        self.particles
            .iter()
            .filter(move |particle| particle.is_pid(z, a))
    }

    /// Count the number of particles in the event with the specified Z and A values.
    pub fn pid_mult(&self, z: usize, a: usize) -> usize {
        self.particles.iter().filter(|p| p.is_pid(z, a)).count()
    }

    pub fn z_mult(&self, z: usize) -> usize {
        self.particles.iter().filter(|p| p.Z() == z).count()
    }

    /// Add a particle to the list by specifying the particle's charge number, mass number, momentum, and which detector measured it.
    pub fn add_particle_by_properties(
        &mut self,
        charge_num: usize,
        mass_num: usize,
        px: f64,
        py: f64,
        pz: f64,
        detector: usize,
    ) {
        self.particles.push(Particle {
            charge_num,
            mass_num,
            momentum: PhysVec::new([px, py, pz]),
            detector,
        })
    }

    /// Clone a particle into the particle list.
    pub fn add_particle_by_clone(&mut self, particle: &Particle) -> &mut Self {
        self.particles.push(particle.clone());
        self
    }

    /// Move a particle to the list.
    pub fn add_particle(&mut self, particle: Particle) -> &mut Self {
        self.particles.push(particle);
        self
    }

    /// Retrieve a random particle from the event.
    pub fn random_particle(&self) -> &Particle {
        self.particles.choose(&mut rand::thread_rng()).unwrap()
    }

    /// Retrieve a random particle matching the specified Z and A values.
    pub fn random_particle_by_pid(&self, z: usize, a: usize) -> Option<&Particle> {
        self.iter_by_pid(z, a).choose(&mut rand::thread_rng())
    }

    /// Create a collected vector of particle references matching the specified Z and A values.
    pub fn particles_by_pid(&self, z: usize, a: usize) -> Vec<&Particle> {
        self.particles.iter().filter(|p| p.is_pid(z, a)).collect()
    }
}

/// Represents a particle in an event. Contains the particle's momentum, charge number, mass number, and detector number.
/// The momentum is stored as a PhysVec, which is a 3D vector with units of MeV/c.
#[derive(Default, Debug, Clone, Archive, Serialize, Deserialize)]
pub struct Particle {
    momentum: PhysVec,
    charge_num: usize,
    mass_num: usize,
    detector: usize,
}

impl Particle {
    /// Charge number of the particle.
    pub fn proton_num(&self) -> usize {
        self.charge_num
    }

    /// Charge number of the particle.
    #[allow(non_snake_case)]
    pub fn Z(&self) -> usize {
        self.proton_num()
    }

    /// Mass number of the particle.
    pub fn mass_num(&self) -> usize {
        self.mass_num
    }
    /// Mass number of the particle.
    #[allow(non_snake_case)]
    pub fn A(&self) -> usize {
        self.mass_num()
    }

    /// Neutron number of the particle.
    pub fn neutron_number(&self) -> usize {
        self.mass_num - self.charge_num
    }
    /// Neutron number of the particle.
    #[allow(non_snake_case)]
    pub fn N(&self) -> usize {
        self.neutron_number()
    }

    /// Mass of the particle in MeV/c^2.
    #[allow(non_snake_case)]
    pub fn mass_MeV_per_c2(&self) -> f64 {
        NUCLEAR_DB.nuclear_mass_MeV_per_c2(self.charge_num, self.mass_num)
    }
    /// Mass of the particle in MeV/c^2.
    #[allow(non_snake_case)]
    pub fn m(&self) -> f64 {
        self.mass_MeV_per_c2()
    }

    /// Detector index that the particle was measured by.
    pub fn detector(&self) -> usize {
        self.detector
    }

    /// Momentum of the particle in MeV/c.
    #[allow(non_snake_case)]
    pub fn momentum_MeV_per_c(&self) -> &PhysVec {
        &self.momentum
    }

    #[allow(non_snake_case)]
    pub fn classical_momentum_MeV_per_c(&self) -> &PhysVec {
        &self.momentum
    }

    /// Momentum of the particle in MeV/c.
    pub fn p(&self) -> &PhysVec {
        &self.momentum
    }

    /// Relativistic energy of the particle in MeV.
    #[allow(non_snake_case)]
    pub fn energy_MeV(&self) -> f64 {
        (self.momentum_MeV_per_c().mag_sqr() + self.mass_MeV_per_c2().powi(2)).sqrt()
    }

    /// Lorentz factor of the particle, often denoted as gamma.
    pub fn lorentz_factor(&self) -> f64 {
        (1. + self.momentum_MeV_per_c().mag_sqr() / self.mass_MeV_per_c2().powi(2)).sqrt()
    }

    /// Velocity of the particle in c.
    pub fn velocity_c(&self) -> PhysVec {
        self.momentum
            .as_scaled_by(1. / self.lorentz_factor() / self.mass_MeV_per_c2())
    }

    pub fn classical_velocity_c(&self) -> PhysVec {
        self.momentum.as_scaled_by(1. / self.mass_MeV_per_c2())
    }

    /// Velocity of the particle in c.
    pub fn v(&self) -> PhysVec {
        self.velocity_c()
    }

    /// Kinetic energy of the particle in MeV.
    #[allow(non_snake_case)]
    pub fn kinetic_energy_MeV(&self) -> f64 {
        (self.lorentz_factor() - 1.) * self.mass_MeV_per_c2()
    }

    #[allow(non_snake_case)]
    pub fn classical_kinetic_energy_MeV(&self) -> f64 {
        0.5 / self.mass_MeV_per_c2() * self.classical_momentum_MeV_per_c().mag_sqr()
    }

    /// Boost the particles by a velocity vector.
    pub fn boost(&mut self, velocity: &PhysVec) -> &mut Self {
        // gamma of the frame shift, not of the particle in the original frame
        let gamma = 1. / (1. - velocity.mag_sqr()).sqrt();

        //let v_old_frame = self.velocity_c();
        let u = &self.velocity_c();
        let v = velocity;

        let u_prime = 1. / (1. - v.dot(u)) * (u / gamma - v + gamma / (gamma + 1.) * v.dot(u) * v);

        //let v_new_frame = (&v_old_frame / gamma - velocity
        //+ velocity.dot(&v_old_frame) * gamma / (gamma + 1.) * velocity)
        // * 1. / (1. - velocity.dot(&v_old_frame));
        let new_gamma = 1. / (1. - u_prime.mag_sqr()).sqrt();
        self.momentum = u_prime.as_scaled_by(new_gamma * self.mass_MeV_per_c2());

        self
    }

    pub fn classical_boost(&mut self, velocity: &PhysVec) -> &mut Self {
        let v_old_frame = self.classical_velocity_c();
        let v_new_frame = &v_old_frame - velocity;
        self.momentum = v_new_frame.as_scaled_by(self.mass_MeV_per_c2());
        self
    }

    /// Checks if the particle identification matches the passed Z and A values.
    pub fn is_pid(&self, z: usize, a: usize) -> bool {
        self.charge_num == z && self.mass_num == a
    }

    /// returns the particle's particle identification as (Z, A)
    pub fn pid(&self) -> (usize, usize) {
        (self.charge_num, self.mass_num)
    }

    /// Change the angle of the momentum vector in the x-y plane
    pub fn set_phi_rad(&mut self, phi: f64) -> &mut Self {
        self.momentum.set_phi_rad(phi);
        self
    }

    /// Change the angle of the momentum vector in the x-y plane
    pub fn set_phi_deg(&mut self, phi: f64) -> &mut Self {
        self.momentum.set_phi_deg(phi);
        self
    }

    /// Set the momentum vector  
    pub fn set_momentum(&mut self, momentum: &PhysVec) -> &mut Self {
        self.momentum = PhysVec::new([momentum.x(), momentum.y(), momentum.z()]);
        self
    }
}
