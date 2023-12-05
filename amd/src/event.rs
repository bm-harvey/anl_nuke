//use std::path::Path;

use std::cmp::Ordering;

use crate::phys::PhysVec;

use rkyv::{Archive, Deserialize, Serialize};

#[derive(Archive, Debug, Serialize, Deserialize)]
#[archive(check_bytes)]
pub struct Nucleon {
    /// The unique specifier for the nucleon. This identifier connects the nucleon through
    /// `TimeSteps`.
    id: u8,
    /// Every `Nucleon` is either a proton or a neutron.
    proton: bool,
    /// Every `Nucleon is either spin up or down.
    spin_up: bool,
}

impl Nucleon {
    /// New requires explicit passing of all of the required data to fully describe the nucleon
    pub fn new(id: u8, proton: bool, spin_up: bool) -> Self {
        Self {
            id,
            proton,
            spin_up,
        }
    }

    pub fn is_proton(&self) -> bool {
        self.proton
    }
    pub fn is_neutron(&self) -> bool {
        !self.is_proton()
    }
    pub fn is_spin_up(&self) -> bool {
        self.spin_up
    }
    pub fn is_spin_down(&self) -> bool {
        !self.is_spin_up()
    }

    /// Unique identifier for the nucleon. Might not be the only nucleon with this ID. Nucleons of
    /// different spins and types can have overlapping IDs. This value will be shared across timesteps.
    pub fn id(&self) -> u8 {
        self.id
    }

    /// In some readouts of AMD, the id is only unique for a given nucleon type. This checks
    /// everything in case one is analyszing one of those formats.If you know the nucleon ID is
    /// perfectly unique within a time step, comparing `id` is sufficient.
    pub fn is_same_nucleon(&self, other: &Nucleon) -> bool {
        (self.id() == other.id())
            && (self.is_spin_up() == other.is_spin_up())
            && (self.is_proton() == other.is_proton())
    }
}

/// Essentially a storage struct for a group of `Nucleon`. Additionally stores basic information
/// that can be calculated at run time, but easier and faster to just store in this case.
#[derive(Archive, Debug, Serialize, Deserialize)]
#[archive(check_bytes)]
pub struct Fragment {
    charge_num: u8,
    mass_num: u8,
    momentum: PhysVec,
    location: PhysVec,
    angular_momentum: PhysVec,
    total_angular_momentum: f32,
    internal_energy: f32,
    max_density: f32,
    nucleons: Vec<Nucleon>,
}

impl Fragment {
    /// This estimates the converstion between mass number and mass in MeV / c.
    const AMU_TO_MEV_PER_C: f32 = 931.5;

    /// ctor for Fragment requires all information of Fragment to be passed
    pub fn new(
        charge_num: u8,
        mass_num: u8,
        momentum: PhysVec,
        location: PhysVec,
        angular_momentum: PhysVec,
        total_angular_momentum: f32,
        internal_energy: f32,
        max_density: f32,
        nucleons: Vec<Nucleon>,
    ) -> Fragment {
        Fragment {
            charge_num,
            mass_num,
            momentum,
            location,
            nucleons,
            max_density,
            total_angular_momentum,
            angular_momentum,
            internal_energy,
        }
    }

    /// Z number of the nucleon
    pub fn charge_num(&self) -> u8 {
        self.charge_num
    }

    /// A number of the nucleon
    pub fn mass_num(&self) -> u8 {
        self.mass_num
    }

    /// N number of a nucleon
    pub fn neutron_num(&self) -> u8 {
        self.mass_num - self.charge_num
    }

    /// Estimated mass of the Fragment in MeV
    pub fn mass(&self) -> f32 {
        (self.mass_num() as f32) * Fragment::AMU_TO_MEV_PER_C
    }

    /// `pid` stands for Particle Identification. Returns a tuple (Z,A). This function is useful
    /// for checking the type of a fragment
    pub fn pid(&self) -> (u8, u8) {
        (self.charge_num, self.mass_num)
    }

    /// Checks whether the passed particle type info matches the stored info. Useful for counting
    /// the number of particles of a certain type in a timestep
    pub fn is_pid(&self, z: u8, a: u8) -> bool {
        (z, a) == self.pid()
    }

    /// Refrence to the stored momentum vector in MeV / c.
    pub fn momentum_vector(&self) -> &PhysVec {
        &self.momentum
    }

    /// See `momentum_vector`
    pub fn p_vec(&self) -> &PhysVec {
        self.momentum_vector()
    }

    /// Reference to the stored coordinate vector in fm.
    pub fn coordinate_vector(&self) -> &PhysVec {
        &self.location
    }

    /// See `coordinate_vector`
    pub fn r_vec(&self) -> &PhysVec {
        self.coordinate_vector()
    }

    /// Constructs and returns a velocity vector
    pub fn v_vec(&self) -> PhysVec {
        self.momentum_vector().as_normalized_by(self.mass())
    }

    /// Checks if a nucleon is in the fragment by comparing ID's
    pub fn contains_nucleon_by_id(&self, nucleon: &Nucleon) -> bool {
        for nuc in self.nucleons.iter() {
            if nuc.is_same_nucleon(nucleon) {
                return true;
            }
        }
        false
    }

    /// Classical kinetic energy of the particle in MeV
    pub fn kinetic_energy(&self) -> f32 {
        0.5 / self.mass() * self.p_vec().mag_sqr()
    }

    /// Kinetic energy per nucleon in MeV
    pub fn kinetic_energy_per_u(&self) -> f32 {
        self.kinetic_energy() / (self.mass_num() as f32)
    }

    /// adjusts the z-coordinate of the momentum by a specified amount.
    fn shift_momentum_z(&mut self, delta_pz: f32) {
        self.momentum.set_z(self.momentum.z() + delta_pz)
    }

    /// add a nucleon
    pub fn add_nucleon(&mut self, nucleon: Nucleon) {
        self.nucleons.push(nucleon);
    }

    pub fn mass_cmp(&self, other:&Fragment) ->Ordering {
        if self.charge_num() != other.charge_num(){
           return self.charge_num().cmp(&other.charge_num());
        }
        if self.mass_num() != other.mass_num(){
           return self.mass_num().cmp(&other.mass_num());
        }

        return self.p_vec().z().total_cmp(&other.p_vec().z());
    }


    pub fn nucleons(&self) -> &Vec<Nucleon> {
        &self.nucleons
    }
}


/// Acts as a storage container for `Fragments`.  
#[derive(Archive, Debug, Serialize, Deserialize)]
#[archive(check_bytes)]
pub struct TimeStep {
    /// List of `Fragment`
    fragments: Vec<Fragment>,
    /// Time in fm / c that this timestep represents
    time: f32,
}

impl TimeStep {
    /// Creates a `Timestep` with a fixed time and fragment list
    pub fn new(time: f32) -> TimeStep {
        Self {
            fragments: Vec::new(),
            time,
        }
    }

    /// Counts the number of `Fragment`s in the timestep
    pub fn multiplicity(&self) -> usize {
        self.fragments.len()
    }

    /// returns a reference to the largest fragment in the timestep
    pub fn largest_fragment(&self) -> &Fragment {
        self.fragment_iter().max_by(|f1, f2| f1.mass_cmp(f2)).unwrap()
    }

    /// Returns the stored time in fm / c
    pub fn time(&self) -> f32 {
        self.time
    }

    /// returns the fragment list by reference.
    pub fn fragments(&self) -> &Vec<Fragment> {
        &self.fragments
    }

    /// Directly returns an iterator over the fragments.
    pub fn fragment_iter(&self) -> std::slice::Iter<'_, Fragment> {
        self.fragments().iter()
    }

    /// Returns a fragment by index; will panic if `idx > self.fragmetnts.len()`. For most
    /// applications, iterating over the list of the fragments is safer and probably faster, but
    /// this is useful if one is extracting information about a specific fragment (and useful for
    /// testing).
    pub fn fragment(&self, idx: usize) -> &Fragment {
        &self.fragments[idx]
    }

    /// Sum of the kinetic energy of all fragments in MeV
    pub fn total_kinetic_energy(&self) -> f32 {
        self.fragment_iter().map(|frag| frag.kinetic_energy()).sum()
    }

    /// Sum of the z-axis momentum in MeV/c. If momentum is being properly conserved, this should
    /// be very close to 0.
    pub fn total_parallel_momentum(&self) -> f32 {
        self.fragment_iter().map(|frag| frag.p_vec().z()).sum()
    }

    /// Sum of the z-axis momentum magnitudes in MeV/c. Even if momentum is well conserved, this value might be very far from 0.
    pub fn total_parallel_momentum_mag(&self) -> f32 {
        self.fragment_iter()
            .map(|frag| frag.p_vec().z().abs())
            .sum()
    }

    /// Loop over all particles and adjust their momentum by an amount that corresponds to the
    /// passed velocity (in units of c)
    pub fn shift_frame(&mut self, beta_z: f32) {
        for frag in self.fragments.iter_mut() {
            let delta_pz = beta_z * frag.mass();
            frag.shift_momentum_z(delta_pz);
        }
    }

    pub fn add_fragment(&mut self, fragment: Fragment) {
        self.fragments.push(fragment)
    }
}

/// Represents one AMD nuclear collision. Most use cases will only care about the data stored in
/// the last TimeStep, but some applications care about fragment tracking through all of the
/// TimeSteps
#[derive(Archive, Debug, Serialize, Deserialize)]
#[archive(check_bytes)]
pub struct Event {
    /// Offset of the nuclei initially in fm
    impact_parameter: f32,
    /// List of timesteps from the AMD.
    time_steps: Vec<TimeStep>,
}

impl Event {
    /// Event requires all information to be set up prior to construction
    pub fn new(impact_parameter: f32, time_steps: Vec<TimeStep>) -> Self {
        Self {
            time_steps,
            impact_parameter,
        }
    }

    /// Returns the impact parameter in fm.
    pub fn impact_parameter(&self) -> f32 {
        self.impact_parameter
    }

    /// Returns a reference to the TimeStep list
    pub fn time_steps(&self) -> &Vec<TimeStep> {
        &self.time_steps
    }
    
    /// Returns the number of timesteps in the events  
    pub fn num_time_steps(&self) -> usize {
        self.time_steps.len()
    }

    /// Returns a time_step reference by index. Will panic! if `idx < self.time_steps.len()`
    pub fn time_step(&self, idx: usize) -> &TimeStep {
        &self.time_steps[idx]
    }

    /// Returns a reference to the last time_step. This will panic if there is no timesteps, but
    /// that would mean one is trying to analyze no data.
    pub fn last_time_step(&self) -> &TimeStep {
        self.time_steps.last().unwrap()
    }

    /// Returns a direct iterator over the timesteps
    pub fn time_step_iter(&self) -> std::slice::Iter<'_, TimeStep> {
        self.time_steps.iter()
    }

    /// Adjust the frame of the event in all time steps by a passed velocity in units of c.  
    pub fn shift_frame(&mut self, beta_z: f32) {
        for ts in self.time_steps.iter_mut() {
            ts.shift_frame(beta_z);
        }
    }
}
