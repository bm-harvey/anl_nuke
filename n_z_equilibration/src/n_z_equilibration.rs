use std::path::Path;

//use event::*;
use amd::event::{Event, Fragment, TimeStep};

use nukers::anl_module::*;

use rmp::encode::RmpWrite;

#[derive(serde::Serialize, Debug, Default)]
pub struct NZEquilibrationData {
    first_collision_time: Option<f32>,
    plf_tlf_separation_time: Option<f32>,
    plf_breakup_time: Option<f32>,
    impact_parameter: f32,
    final_multiplicity: usize,
    plf_frag_z: Option<u8>,
    plf_frag_a: Option<u8>,
    heavier_frag_z: Option<u8>,
    heavier_frag_a: Option<u8>,
    lighter_frag_z: Option<u8>,
    lighter_frag_a: Option<u8>,
    angular_alignment: Option<f32>,
    delta_lf: Option<f32>,
    delta_hf: Option<f32>,
    heavier_frag_exp_z: Option<u8>,
    heavier_frag_exp_a: Option<u8>,
    lighter_frag_exp_z: Option<u8>,
    lighter_frag_exp_a: Option<u8>,
    angular_alignment_exp: Option<f32>,
}

impl NZEquilibrationData {
    pub fn set_first_collision_time(&mut self, value: f32) -> &mut Self {
        self.first_collision_time = Some(value);
        self
    }

    pub fn set_plf_tlf_separation_time(&mut self, value: f32) -> &mut Self {
        self.plf_tlf_separation_time = Some(value);
        self
    }
    pub fn set_plf_breakup_time(&mut self, value: f32) -> &mut Self {
        self.plf_breakup_time = Some(value);
        self
    }

    pub fn set_alpha(&mut self, value: f32) -> &mut Self {
        self.angular_alignment = Some(value);
        self
    }

    pub fn set_delta_lf(&mut self, value: f32) -> &mut Self {
        self.delta_lf = Some(value);
        self
    }
    pub fn set_delta_hf(&mut self, value: f32) -> &mut Self {
        self.delta_hf = Some(value);
        self
    }

    pub fn set_impact_parameter(&mut self, value: f32) -> &mut Self {
        self.impact_parameter = value;
        self
    }

    pub fn set_final_multiplicity(&mut self, value: usize) -> &mut Self {
        self.final_multiplicity = value;
        self
    }

    pub fn set_plf(&mut self, fragment: &Fragment) -> &mut Self {
        self.plf_frag_z = Some(fragment.charge_num());
        self.plf_frag_a = Some(fragment.mass_num());
        self
    }

    pub fn set_lighter_fragment(&mut self, fragment: &Fragment) -> &mut Self {
        self.lighter_frag_z = Some(fragment.charge_num());
        self.lighter_frag_a = Some(fragment.mass_num());
        self
    }
    pub fn set_heavier_fragment(&mut self, fragment: &Fragment) -> &mut Self {
        self.heavier_frag_z = Some(fragment.charge_num());
        self.heavier_frag_a = Some(fragment.mass_num());
        self
    }
}

#[derive(Default, Debug, serde::Serialize)]
pub struct NZEquilibration {
    light_frag_min_z: Option<u8>,
    heavy_frag_min_z: Option<u8>,
    light_frag_max_z: Option<u8>,
    heavy_frag_max_z: Option<u8>,
    combined_nucleus_min_mass: Option<u8>,
    plf_min_charge: Option<u8>,
    tlf_min_charge: Option<u8>,
    eventwise_data: Vec<NZEquilibrationData>,
    time_buffer: f32,
    extend_fragments_to_end_of_event: bool,
}
impl<'a> NZEquilibration {
    pub fn new() -> Self {
        NZEquilibration {
            light_frag_min_z: Some(4),
            heavy_frag_min_z: Some(11),
            light_frag_max_z: Some(20),
            heavy_frag_max_z: Some(35),
            combined_nucleus_min_mass: Some(100),
            plf_min_charge: Some(15),
            tlf_min_charge: Some(15),
            eventwise_data: Vec::new(),
            time_buffer: 0.,
            extend_fragments_to_end_of_event: false,
        }
    }

    pub fn boxed(self) -> Box<Self> {
        Box::new(self)
    }

    pub fn track_final_fragments_to_end_of_event(mut self, track: bool) -> Self {
        self.extend_fragments_to_end_of_event = track;
        self
    }
    pub fn with_time_buffer(mut self, time_buffer: f32) -> Self {
        self.time_buffer = time_buffer;
        self
    }
    pub fn with_light_frag_max_charge(mut self, max_z: u8) -> Self {
        self.light_frag_max_z = Some(max_z);
        self
    }

    pub fn with_heavy_frag_max_charge(mut self, max_z: u8) -> Self {
        self.heavy_frag_max_z = Some(max_z);
        self
    }
    pub fn with_light_frag_min_charge(mut self, min_z: u8) -> Self {
        self.light_frag_min_z = Some(min_z);
        self
    }

    pub fn with_heavy_frag_min_charge(mut self, min_z: u8) -> Self {
        self.heavy_frag_min_z = Some(min_z);
        self
    }

    pub fn with_combined_nucleus_min_mass(mut self, min_a: u8) -> Self {
        self.combined_nucleus_min_mass = Some(min_a);
        self
    }

    pub fn with_plf_min_charge(mut self, min_a: u8) -> Self {
        self.plf_min_charge = Some(min_a);
        self
    }

    pub fn with_tlf_min_charge(mut self, min_a: u8) -> Self {
        self.tlf_min_charge = Some(min_a);
        self
    }

    pub fn get_mass_sorted_fragment_indices(time_step: &TimeStep) -> Vec<usize> {
        let mut result = (0..time_step.fragments().len()).collect::<Vec<_>>();
        result.sort_by(|idx_1, idx_2| {
            let f1 = time_step.fragment(*idx_1);
            let f2 = time_step.fragment(*idx_2);
            // explicitly reverse the order so that the returned indices go from large -> small
            // fragments
            Fragment::mass_cmp(f2, f1)
        });
        result
    }

    pub fn get_commonality_sorted_fragment_indices(
        time_step: &TimeStep,
        fragment: &Fragment,
    ) -> Vec<usize> {
        let mut result = (0..time_step.fragments().len()).collect::<Vec<_>>();
        result.sort_by(|&idx_1, &idx_2| {
            let f1 = time_step.fragment(idx_1);
            let f2 = time_step.fragment(idx_2);
            let count_1 = NZEquilibration::num_shared_nucleons(f1, fragment);
            let count_2 = NZEquilibration::num_shared_nucleons(f2, fragment);
            // explicitly reverse the order so the most similar fragments are presented first
            count_2.cmp(&count_1)
        });
        result
    }

    pub fn fragments_share_a_nucleon(fragment_1: &Fragment, fragment_2: &Fragment) -> bool {
        for nucleon in fragment_1.nucleons().iter() {
            if fragment_2.contains_nucleon_by_id(nucleon) {
                return true;
            }
        }
        false
    }

    pub fn get_heavy_fragment_indices_in_time_step(
        fragment: &Fragment,
        other_time_step: &TimeStep,
    ) -> Vec<usize> {
        NZEquilibration::get_mass_sorted_fragment_indices(other_time_step)
            .iter()
            .filter(|frag_idx| {
                NZEquilibration::fragments_share_a_nucleon(
                    other_time_step.fragment(**frag_idx),
                    fragment,
                )
            })
            .map(|x| *x)
            .collect()
    }

    pub fn get_similar_fragment_indices_in_time_step(
        fragment: &Fragment,
        other_time_step: &TimeStep,
    ) -> Vec<usize> {
        NZEquilibration::get_mass_sorted_fragment_indices(other_time_step)
            .iter()
            .filter(|frag_idx| {
                NZEquilibration::fragments_share_a_nucleon(
                    other_time_step.fragment(**frag_idx),
                    fragment,
                )
            })
            .map(|x| *x)
            .collect()
    }

    pub fn get_largest_parent(
        fragment: &Fragment,
        time_idx: usize,
        event: &'a Event,
    ) -> Option<&'a Fragment> {
        if time_idx == 0 {
            return None;
        }

        let indices = NZEquilibration::get_heavy_fragment_indices_in_time_step(
            fragment,
            event.time_step(time_idx - 1),
        );

        Some(
            event
                .time_step(time_idx - 1)
                .fragment(*indices.first().unwrap()),
        )

        //NZEquilibration::get_children_indices(
    }
    pub fn num_shared_nucleons(frag1: &Fragment, frag2: &Fragment) -> usize {
        let mut count = 0;

        frag1.nucleons().iter().for_each(|nucleon| {
            if frag2.contains_nucleon_by_id(nucleon) {
                count += 1;
            }
        });
        count
    }

    pub fn get_most_similar_parent(
        fragment: &Fragment,
        time_idx: usize,
        event: &'a Event,
    ) -> Option<&'a Fragment> {
        if time_idx == 0 {
            return None;
        }

        let indices = NZEquilibration::get_similar_fragment_indices_in_time_step(
            fragment,
            event.time_step(time_idx - 1),
        );

        let result = indices.iter().max_by(|idx_1, idx_2| {
            let num_1 = NZEquilibration::num_shared_nucleons(
                fragment,
                event.time_step(time_idx - 1).fragment(**idx_1),
            );
            let num_2 = NZEquilibration::num_shared_nucleons(
                fragment,
                event.time_step(time_idx - 1).fragment(**idx_2),
            );
            num_1.cmp(&num_2)
        });

        Some(event.time_step(time_idx - 1).fragment(*result.unwrap()))
    }
    pub fn get_most_similar_child(
        fragment: &Fragment,
        time_idx: usize,
        event: &'a Event,
    ) -> Option<&'a Fragment> {
        if time_idx + 1 >= event.num_time_steps() {
            return None;
        }

        let indices = NZEquilibration::get_similar_fragment_indices_in_time_step(
            fragment,
            event.time_step(time_idx + 1),
        );

        Some(event.time_step(time_idx + 1).fragment(indices[0]))
    }
    pub fn get_largest_child(
        fragment: &Fragment,
        time_idx: usize,
        event: &'a Event,
    ) -> Option<&'a Fragment> {
        if time_idx + 1 >= event.num_time_steps() {
            return None;
        }

        let indices = NZEquilibration::get_heavy_fragment_indices_in_time_step(
            fragment,
            event.time_step(time_idx + 1),
        )
        .iter()
        .map(|&x| x)
        .collect::<Vec<usize>>();

        Some(
            event
                .time_step(time_idx + 1)
                .fragment(*indices.first().unwrap()),
        )
    }

    pub fn get_descendent(
        fragment: &'a Fragment,
        time_idx: usize,
        event: &'a Event,
    ) -> &'a Fragment {
        let largest_child = NZEquilibration::get_most_similar_child(fragment, time_idx, event);
        match largest_child {
            None => fragment,
            Some(child) => NZEquilibration::get_descendent(child, time_idx + 1, event),
        }
    }

    pub fn get_common_ancestor(
        fragment_1: &'a Fragment,
        fragment_2: &'a Fragment,
        time_idx: usize,
        event: &'a Event,
    ) -> (Option<&'a Fragment>, Option<usize>) {
        let fragment_1_parent =
            NZEquilibration::get_most_similar_parent(fragment_1, time_idx, event);
        let fragment_2_parent =
            NZEquilibration::get_most_similar_parent(fragment_2, time_idx, event);

        if fragment_1_parent.is_none() || fragment_2_parent.is_none() {
            (None, None)
        } else {
            let fragment_1_parent = fragment_1_parent.unwrap();
            let fragment_2_parent = fragment_2_parent.unwrap();

            if std::ptr::eq(fragment_1_parent, fragment_2_parent) {
                (Some(fragment_1_parent), Some(time_idx))
            } else {
                NZEquilibration::get_common_ancestor(
                    fragment_1_parent,
                    fragment_2_parent,
                    time_idx - 1,
                    event,
                )
            }
        }
    }

    pub fn angular_alignment(fragment_1: &Fragment, fragment_2: &Fragment) -> f32 {
        let velocity_hf = fragment_1.v_vec();
        let velocity_lf = fragment_2.v_vec();

        let velocity_rel = &velocity_hf - &velocity_lf;

        let momentum_cm = fragment_1.p_vec() + fragment_2.p_vec();
        let velocity_cm = momentum_cm.as_normalized_by(fragment_1.mass() + fragment_2.mass());

        velocity_rel.inner_angle_deg(&velocity_cm)
    }

    pub fn delta(fragment: &Fragment) -> f32 {
        ((fragment.neutron_num() as f32) - (fragment.charge_num() as f32))
            / (fragment.mass_num() as f32)
    }
}

impl AnlModule<Event> for NZEquilibration {
    fn name(&self) -> String {
        "n_z_equilibration".into()
    }

    fn initialize(&mut self, _output_directory: &Path) {
        if self.combined_nucleus_min_mass.is_none() {
            println!("Missing a setting for initial collision required mass");
            panic!();
        }
        if self.light_frag_min_z.is_none() {
            println!("Missing a setting for minimum Z of lighter fragment");
            panic!();
        }
        if self.heavy_frag_min_z.is_none() {
            println!("Missing a setting for minimum Z of heavier fragment");
            panic!();
        }
        if self.light_frag_max_z.is_none() {
            println!("Missing a setting for maximum Z of lighter fragment");
            panic!();
        }
        if self.heavy_frag_max_z.is_none() {
            println!("Missing a setting for maximum Z of heavier fragment");
            panic!();
        }
        if self.plf_min_charge.is_none() {
            println!("Missing a setting for PLF minimum mass");
            panic!();
        }
        if self.tlf_min_charge.is_none() {
            println!("Missing a setting for TLF minimum mass");
            panic!();
        }
    }

    fn analyze_event(&mut self, event: &Event, _event_idx: usize) {
        let mut anl_result = NZEquilibrationData::default();
        anl_result.set_impact_parameter(event.impact_parameter());
        anl_result.set_final_multiplicity(event.last_time_step().multiplicity());

        // the analysis from an experimental perspective
        let final_time_step = event.last_time_step();
        let final_particle_idices =
            NZEquilibration::get_mass_sorted_fragment_indices(final_time_step);

        let forward_particles = final_particle_idices
            .iter()
            .filter(|idx| final_time_step.fragment(**idx).p_vec().z() > 0.)
            .collect::<Vec<_>>();

        if forward_particles.len() >= 2 {
            let experimental_hf = final_time_step.fragment(*forward_particles[0]);
            let experimental_lf = final_time_step.fragment(*forward_particles[1]);
            if experimental_lf.charge_num() >= self.light_frag_min_z.unwrap() {
                anl_result.angular_alignment_exp = Some(NZEquilibration::angular_alignment(
                    experimental_hf,
                    experimental_lf,
                ));
                anl_result.heavier_frag_exp_z = Some(experimental_hf.charge_num());
                anl_result.heavier_frag_exp_a = Some(experimental_hf.mass_num());
                anl_result.lighter_frag_exp_z = Some(experimental_lf.charge_num());
                anl_result.lighter_frag_exp_a = Some(experimental_lf.mass_num());
            }
            //anl_result.angular_alignment_exp = NZEquilibration::angular_alignment(experimental_hf, experimental_lf);
        }

        // find the initial collision time by finding when the largest fragment in the event has a
        // sufficient mass
        let combined_nucleus_min_mass = self.combined_nucleus_min_mass.unwrap();
        let mut first_collision_time_idx: Option<usize> = None;
        for (idx, ts) in event.time_step_iter().enumerate() {
            if ts.largest_fragment().mass_num() >= combined_nucleus_min_mass {
                first_collision_time_idx = Some(idx);
                anl_result.set_first_collision_time(ts.time());
                break;
            }
        }
        // If there is not a collision (the nuclei just simply missed each other), then return the
        // data we have, and don't proceed further... there is no sense looking for a breakup of a
        // combined nucleus that doesn't exist
        if first_collision_time_idx.is_none() {
            self.eventwise_data.push(anl_result);
            return;
        }

        // Now we want to know when the combined nucleus breaks up into two large fragments. Care
        // is taken to ensure that the large fragments don't recombine into the single large
        // fragment again.
        let first_collision_time_idx = first_collision_time_idx.unwrap();
        let mut last_collision_time_idx = first_collision_time_idx;
        for time_idx in first_collision_time_idx..event.num_time_steps() {
            let time_step = event.time_step(time_idx);
            if time_step.largest_fragment().mass_num() >= combined_nucleus_min_mass {
                last_collision_time_idx = time_idx;
            }
        }
        let plf_min_charge = self.plf_min_charge.unwrap();
        let tlf_min_charge = self.tlf_min_charge.unwrap();
        let mut plf_tlf_seperation_time: Option<f32> = None;
        let mut plf_tlf_seperation_time_idx: usize = last_collision_time_idx;
        let mut plf: Option<&Fragment> = None;
        let mut _tlf: Option<&Fragment> = None;
        for time_idx in (last_collision_time_idx + 1)..event.time_steps().len() {
            let time_step = event.time_step(time_idx);
            if time_step.multiplicity() >= 2 {
                let frag_indicies = NZEquilibration::get_mass_sorted_fragment_indices(time_step);
                plf = Some(time_step.fragment(frag_indicies[0]));
                _tlf = Some(time_step.fragment(frag_indicies[1]));

                if plf.unwrap().p_vec().z() * _tlf.unwrap().p_vec().z() > 0. {
                    // the two largest fragments are going in the same direction along the beamline...not what
                    // are looking for
                    continue;
                }

                let plf_descendent =
                    NZEquilibration::get_descendent(&plf.unwrap(), time_idx, event);
                let tlf_descendent =
                    NZEquilibration::get_descendent(&_tlf.unwrap(), time_idx, event);

                if std::ptr::eq(plf_descendent, tlf_descendent) {
                    continue;
                }
                let merge_time = NZEquilibration::get_common_ancestor(
                    plf_descendent,
                    tlf_descendent,
                    time_idx,
                    event,
                )
                .1;
                if merge_time.is_none() {
                    continue;
                }
                if merge_time.unwrap() != time_idx - 1 {
                    continue;
                }
                //println!("merge time {}", merge_time.unwrap());
                //println!("current time {}", time_idx);

                if plf.unwrap().p_vec().z() < 0. {
                    (_tlf, plf) = (plf, _tlf);
                }

                if plf.unwrap().charge_num() >= plf_min_charge
                    && _tlf.unwrap().charge_num() >= tlf_min_charge
                {
                    plf_tlf_seperation_time = Some(time_step.time());
                    plf_tlf_seperation_time_idx = time_idx;
                    anl_result.set_plf_tlf_separation_time(time_step.time());
                    anl_result.set_plf(plf.unwrap());
                    break;
                } else {
                    plf = None;
                    _tlf = None;
                }
            }
        }
        // If there as not a PLF and TLF formed, then it does not make sense to look for when they
        // broke apart.
        if plf_tlf_seperation_time.is_none() {
            self.eventwise_data.push(anl_result);
            return;
        }

        // Look for the PLF breaking apart
        let mut plf = plf.unwrap();
        for time_idx in plf_tlf_seperation_time_idx..(event.num_time_steps() - 1) {
            if event.last_time_step().time() - event.time_step(time_idx).time() < self.time_buffer {
                break;
            }
            let _current_time_step = event.time_step(time_idx);
            let next_time_step = event.time_step(time_idx + 1);

            let children_indices =
                NZEquilibration::get_similar_fragment_indices_in_time_step(plf, next_time_step);
            if children_indices.len() == 1 {
                // there is only one child, so look again in the next timestep.
                plf = next_time_step.fragment(children_indices[0]);
                continue;
            }

            //there must be at least 2 children here.
            let heavier_frag = next_time_step.fragment(children_indices[0]);
            for child_idx in 1..children_indices.len() {
                let lighter_frag = next_time_step.fragment(children_indices[child_idx]);
                if lighter_frag.charge_num() >= self.light_frag_min_z.unwrap()
                    && lighter_frag.charge_num() <= self.light_frag_max_z.unwrap()
                    && heavier_frag.charge_num() >= self.heavy_frag_min_z.unwrap()
                    && heavier_frag.charge_num() <= self.heavy_frag_max_z.unwrap()
                {
                    let hf_descendent =
                        NZEquilibration::get_descendent(&heavier_frag, time_idx + 1, event);
                    let lf_descendent =
                        NZEquilibration::get_descendent(&lighter_frag, time_idx + 1, event);

                    if std::ptr::eq(hf_descendent, lf_descendent) {
                        continue;
                    }

                    let common_ancestor = NZEquilibration::get_common_ancestor(
                        hf_descendent,
                        lf_descendent,
                        event.num_time_steps(),
                        event,
                    )
                    .0
                    .unwrap();

                    if !std::ptr::eq(plf, common_ancestor) {
                        continue;
                    }
                    if heavier_frag.p_vec().z() * lighter_frag.p_vec().z() < 0. {
                        continue;
                    }

                    // A PLF breakup was found.
                    if self.extend_fragments_to_end_of_event {
                        anl_result.set_heavier_fragment(hf_descendent);
                        anl_result.set_lighter_fragment(lf_descendent);
                        anl_result.set_plf_breakup_time(next_time_step.time());
                        anl_result.set_alpha(NZEquilibration::angular_alignment(
                            hf_descendent,
                            lf_descendent,
                        ));
                        anl_result.set_delta_lf(NZEquilibration::delta(lf_descendent));
                        anl_result.set_delta_hf(NZEquilibration::delta(hf_descendent));
                    } else {
                        anl_result.set_heavier_fragment(heavier_frag);
                        anl_result.set_lighter_fragment(lighter_frag);
                        anl_result.set_plf_breakup_time(next_time_step.time());
                        anl_result.set_alpha(NZEquilibration::angular_alignment(
                            heavier_frag,
                            lighter_frag,
                        ));
                        anl_result.set_delta_lf(NZEquilibration::delta(lighter_frag));
                        anl_result.set_delta_hf(NZEquilibration::delta(heavier_frag));
                    }
                }
            }
            plf = heavier_frag;
        }

        // Everything was found and we are done with analysis
        self.eventwise_data.push(anl_result);
    } // analyze_event

    fn generate_output(&mut self, output_directory: &Path) {
        let output_file_name: String =
            format!("{}{}.json", output_directory.to_str().unwrap(), self.name());
        let output_file = std::fs::File::create(output_file_name).unwrap();
        let mut out_buf = std::io::BufWriter::with_capacity(100_000, output_file);
        let mut buf = serde_json::to_vec_pretty(&self).unwrap();
        out_buf.write_bytes(&mut buf).unwrap();
    }
}
