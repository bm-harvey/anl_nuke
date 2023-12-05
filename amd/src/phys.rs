use std::f32::consts::PI;

use std::ops::{Add, Sub};

use rkyv::{Archive, Deserialize, Serialize};

#[derive(Archive, Debug, Serialize, Deserialize)]
#[archive(check_bytes)]
pub struct PhysVec {
    coordinate: [f32; 3],
}

impl PhysVec {
    pub fn new(coordinate: [f32; 3]) -> Self {
        PhysVec { coordinate }
    }

    pub fn num_dims(&self) -> usize {
        self.coordinate.len()
    }

    pub fn at(&self, idx: usize) -> f32 {
        self.coordinate[idx]
    }

    pub fn at_or_zero(&self, idx: usize) -> f32 {
        if !self.coordinate.is_empty() {
            self.at(idx)
        } else {
            0.
        }
    }
    pub fn x(&self) -> f32 {
        self.at_or_zero(0)
    }
    pub fn y(&self) -> f32 {
        self.at_or_zero(1)
    }
    pub fn z(&self) -> f32 {
        self.at_or_zero(2)
    }

    pub fn set_by_idx(&mut self, idx: usize, value: f32) {
        self.coordinate[idx] = value;
    }
    pub fn set_x(&mut self, value: f32) {
        self.coordinate[0] = value;
    }
    pub fn set_y(&mut self, value: f32) {
        self.coordinate[1] = value;
    }
    pub fn set_z(&mut self, value: f32) {
        self.coordinate[2] = value;
    }

    pub fn transverse(&self) -> f32 {
        (self.x().powi(2) + self.y().powi(2)).sqrt()
    }

    pub fn parallel(&self) -> f32 {
        self.z().abs()
    }

    pub fn dot(&self, other: &Self) -> f32 {
        (0..self.num_dims())
            .map(|idx| self.coordinate[idx] * other.coordinate[idx])
            .sum()
    }

    pub fn mag_sqr(&self) -> f32 {
        self.dot(self)
    }

    pub fn mag(&self) -> f32 {
        self.mag_sqr().sqrt()
    }

    pub fn as_normalized(&self) -> Self {
        let magnitude = self.mag();

        Self::new(
            (0..self.num_dims())
                .map(|idx| self.at(idx) / magnitude)
                .collect::<Vec<f32>>()
                .try_into()
                .unwrap(),
        )
    }

    pub fn as_normalized_by(&self, magnitude: f32) -> Self {
        Self::new(
            (0..self.num_dims())
                .map(|idx| self.at(idx) / magnitude)
                .collect::<Vec<f32>>()
                .try_into()
                .unwrap(),
        )
    }

    pub fn as_normalized_to(&self, new_magnitude: f32) -> Self {
        let old_magnitude = self.mag();
        let scaling_factor = new_magnitude / old_magnitude;
        Self::new(
            (0..self.num_dims())
                .map(|idx| self.at(idx) * scaling_factor)
                .collect::<Vec<f32>>()
                .try_into()
                .unwrap(),
        )
    }

    pub fn normalize(&mut self) -> &mut Self {
        let norm = self.mag();

        for idx in 0..self.num_dims() {
            self.coordinate[idx] /= norm;
        }
        self
    }

    pub fn theta_rad(&self) -> f32 {
        (self.transverse() / self.z()).atan()
    }

    pub fn phi_rad(&self) -> f32 {
        let mut result = (self.y() / self.x()).atan();
        if self.x().is_sign_negative() {
            result += PI;
        }
        result
    }

    pub fn inner_angle_rad(&self, other: &Self) -> f32 {
        (self.dot(other) / (self.mag() * other.mag())).acos()
    }

    pub fn theta_deg(&self) -> f32 {
        self.theta_rad() * 180. / PI
    }

    pub fn phi_deg(&self) -> f32 {
        self.phi_rad() * 180. / PI
    }

    pub fn inner_angle_deg(&self, other: &Self) -> f32 {
        self.inner_angle_rad(other) * 180. / PI
    }
}

impl<'a, 'b> Add<&'b PhysVec> for &'a PhysVec {
    type Output = PhysVec;
    fn add(self, rhs: &'b PhysVec) -> Self::Output {
        let mut result = PhysVec::new([0., 0., 0.]);

        (0..3).for_each(|idx| result.set_by_idx(idx, self.at(idx) + rhs.at(idx)));

        result
    }
}
impl<'a, 'b> Sub<&'b PhysVec> for &'a PhysVec {
    type Output = PhysVec;
    fn sub(self, rhs: &'b PhysVec) -> Self::Output {
        let mut result = PhysVec::new([0., 0., 0.]);

        (0..3).for_each(|idx| result.set_by_idx(idx, self.at(idx) - rhs.at(idx)));

        result
    }
}
