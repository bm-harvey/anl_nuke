use std::f64::consts::PI;
use std::ops::{Add, AddAssign, Div, Mul, Sub, SubAssign};

#[derive(
    Clone,
    Debug,
    Default,
    rkyv::Archive,
    rkyv::Serialize,
    rkyv::Deserialize,
    serde::Serialize,
    serde::Deserialize,
)]
#[archive(check_bytes)]
/// Represent the coordinates of a vector in 3D space to be used in physics calculations. The
/// coordinates represent laboratory coordinates. In the cartesian representations, the
/// beam-axis is the z-axis, and y-axis is away from gravity; the x-axis is given by the
/// coordinate system being right handed. In the spherical representation, the angles are given
/// with theta in [-pi, pi] and phi in [0, 2 pi). Theta is the angle off of the z-axis, and phi is
/// the azimuthal angle in the x-y plane. The data is stored in a 3 element array containing the
/// cartesian representation.
pub struct PhysVec {
    coordinate: [f64; PhysVec::DIMS],
}

impl PhysVec {
    const DIMS: usize = 3;

    /// Creates a new vector from a given array of coordinates.
    /// The array must have a length of 3.
    /// The given coordinates must be cartesian
    /// ```
    /// use faust::phys_vec::PhysVec;
    /// let vec = PhysVec::new([1., 2., 3.]);
    /// ```
    pub fn new(coordinate: [f64; PhysVec::DIMS]) -> Self {
        Self { coordinate }
    }

    /// Create a new vector from given cartesian coordinates.
    /// ```
    /// # use faust::phys_vec::PhysVec;
    /// let v_1 = PhysVec::from_cartesian(1., 2., 3.);
    /// let v_2 = PhysVec::new([1., 2., 3.]);
    /// assert_eq!(v_1, v_2);
    /// ```
    pub fn from_cartesian(x: f64, y: f64, z: f64) -> Self {
        Self::new([x, y, z])
    }

    /// Create a new vector from given spherical coordinates. Angles should be given in radians.
    /// ```
    /// use faust::phys_vec::PhysVec;
    /// use std::f64::consts::PI;
    ///
    /// let v_1 = PhysVec::from_spherical(1., PI / 2., PI / 4.);
    /// let v_2 = PhysVec::from_cartesian(1. / 2_f64.sqrt(), 1. / 2_f64.sqrt(), 0.);
    ///
    /// assert_eq!(v_1, v_2);
    /// ```
    pub fn from_spherical(mag: f64, theta_rad: f64, phi_rad: f64) -> Self {
        Self::new([
            mag * theta_rad.sin() * phi_rad.cos(),
            mag * theta_rad.sin() * phi_rad.sin(),
            mag * theta_rad.cos(),
        ])
    }

    /// Number of dimensions that this physvec has. Left in for backwards compatibility and in case
    /// `PhysVec` is ever expanded to support 4-Vectors.
    pub fn dims(&self) -> usize {
        self.coordinate.len()
    }

    /// Get the value of the cartesian coordinate at the given index.
    /// ```
    /// # use faust::phys_vec::PhysVec;
    /// let vec = PhysVec::new([1., 2., 3.]);
    /// assert_eq!(vec.at(0), 1.);
    /// assert_eq!(vec.at(1), 2.);
    /// assert_eq!(vec.at(2), 3.);
    /// // Indexes greater than the number of dimensions return 0.
    /// assert_eq!(vec.at(3), 0.);
    /// ```
    pub fn at(&self, idx: usize) -> f64 {
        if idx < Self::DIMS {
            self.coordinate[idx]
        } else {
            0.
        }
    }

    /// Get the x-value of the coordinate.
    /// ```
    /// # use faust::phys_vec::PhysVec;
    /// let vec = PhysVec::new([1., 2., 3.]);
    /// assert_eq!(vec.x(), 1.);
    /// ```
    pub fn x(&self) -> f64 {
        self.at(0)
    }

    /// Get the y-value of the coordinate.
    /// ```
    /// # use faust::phys_vec::PhysVec;
    /// let vec = PhysVec::new([1., 2., 3.]);
    /// assert_eq!(vec.y(), 2.);
    /// ```
    pub fn y(&self) -> f64 {
        self.at(1)
    }

    /// Get the z-value of the coordinate.
    /// ```
    /// # use faust::phys_vec::PhysVec;
    /// let vec = PhysVec::new([1., 2., 3.]);
    /// assert_eq!(vec.z(), 3.);
    /// ```
    pub fn z(&self) -> f64 {
        self.at(2)
    }

    /// Set the value of the coordinate at the given index. This functions has no effect if the
    /// given index is greater than or equal to the number of dimensions.
    /// ```
    /// # use faust::phys_vec::PhysVec;
    /// let mut vec = PhysVec::new([1., 2., 3.]);
    /// vec.set_by_idx(0, 4.);
    /// assert_eq!(vec.x(), 4.);
    /// ```
    pub fn set_by_idx(&mut self, idx: usize, value: f64) -> &mut Self {
        if idx < Self::DIMS {
            self.coordinate[idx] = value;
        }
        self
    }

    /// Set the x-value of the coordinate without affecting y or z.
    /// ```
    /// # use faust::phys_vec::PhysVec;
    /// let mut vec = PhysVec::from_cartesian(1., 2., 3.);
    /// vec.set_x(4.);
    /// assert_eq!(vec.x(), 4.);
    /// ```
    pub fn set_x(&mut self, value: f64) -> &mut Self {
        self.coordinate[0] = value;
        self
    }
    /// Set the y-value of the coordinate without affecting x or z.
    /// ```
    /// # use faust::phys_vec::PhysVec;
    /// let mut vec = PhysVec::from_cartesian(1., 2., 3.);
    /// vec.set_y(4.);
    /// assert_eq!(vec.y(), 4.);
    /// ```
    pub fn set_y(&mut self, value: f64) -> &mut Self {
        self.coordinate[1] = value;
        self
    }

    /// Set the z-value of the coordinate without affecting x or y.
    /// ```
    /// # use faust::phys_vec::PhysVec;
    /// let mut vec = PhysVec::from_cartesian(1., 2., 3.);
    /// vec.set_z(4.);
    /// assert_eq!(vec.z(), 4.);
    /// ```
    pub fn set_z(&mut self, value: f64) -> &mut Self {
        self.coordinate[2] = value;
        self
    }

    /// Get the phi angle in radians.
    /// ```
    /// # use faust::phys_vec::PhysVec;
    /// # use std::f64::consts::PI;
    ///
    /// let vec = PhysVec::from_spherical(1., PI / 3., PI / 4.);
    /// assert!((vec.phi_rad() - PI / 4.).abs() < 1e-10);
    ///
    /// let vec = PhysVec::from_spherical(1.,  PI , 3. * PI / 4.);
    /// assert!((vec.phi_rad() - 3. * PI / 4.).abs() < 1e-10);
    ///```
    pub fn phi_rad(&self) -> f64 {
        let result = (self.x() / self.transverse()).acos() * self.y().signum();
        if self.y() < 0. {
            2. * PI + result
        } else {
            result
        }
    }

    /// Get the phi angle in degrees.
    /// ```
    /// # use faust::phys_vec::PhysVec;
    /// # use std::f64::consts::PI;
    ///
    /// let vec = PhysVec::from_spherical(1., PI / 3., PI / 4.);
    /// assert!((vec.phi_deg() - 45.).abs() < 1e-10);
    ///
    /// let vec = PhysVec::from_spherical(1., PI / 3., 3. * PI / 4.);
    /// assert!((vec.phi_deg() - 135.).abs() < 1e-10);
    /// ```
    pub fn phi_deg(&self) -> f64 {
        self.phi_rad() * 180. / PI
    }

    /// Get the theta angle in radians.
    /// ```
    /// # use faust::phys_vec::PhysVec;
    /// # use std::f64::consts::PI;
    /// let vec = PhysVec::from_spherical(1., PI / 3., PI / 4.);
    /// assert!((vec.theta_rad() - PI / 3.).abs() < 1e-10);
    /// ```
    pub fn theta_rad(&self) -> f64 {
        let mut result = (self.transverse() / self.z().abs()).atan();
        if self.z() < 0. {
            result = PI - result;
        }
        result
    }

    /// Get the theta angle in degrees.
    /// ```
    /// # use faust::phys_vec::PhysVec;
    /// # use std::f64::consts::PI;
    /// let vec = PhysVec::from_spherical(1., PI / 3., PI / 4.);
    /// assert!((vec.theta_deg() - 60.).abs() < 1e-10);
    /// ```
    pub fn theta_deg(&self) -> f64 {
        self.theta_rad() * 180. / PI
    }

    /// Dot-product of two vectors.
    /// ```
    /// # use faust::phys_vec::PhysVec;
    /// let vec_1 = PhysVec::from_cartesian(1., 2., 3.);
    /// let vec_2 = PhysVec::from_cartesian(4., 5., 6.);
    /// assert_eq!(vec_1.dot(&vec_2), 32.);
    /// ```
    pub fn dot(&self, other: &Self) -> f64 {
        (0..PhysVec::DIMS)
            .map(|idx| self.coordinate[idx] * other.coordinate[idx])
            .sum()
    }
    /// Cross-product of two vectors.
    /// ```
    /// # use faust::phys_vec::PhysVec;
    /// let vec_1 = PhysVec::from_cartesian(1., 2., 3.);
    /// let vec_2 = PhysVec::from_cartesian(4., 5., 6.);
    /// let vec_3 = PhysVec::from_cartesian(-3., 6., -3.);
    /// assert_eq!(vec_1.cross(&vec_2), vec_3);
    /// ```
    pub fn cross(&self, other: &Self) -> PhysVec {
        Self::from_cartesian(
            self.y() * other.z() - self.z() * other.y(),
            self.z() * other.x() - self.x() * other.z(),
            self.x() * other.y() - self.y() * other.x(),
        )
    }

    /// Magnitude squared of the vector. This is faster than `mag()` because it does not take the
    /// square root.
    /// ```
    /// # use faust::phys_vec::PhysVec;
    /// let vec = PhysVec::from_cartesian(1., 2., 3.);
    /// assert_eq!(vec.mag_sqr(), 14.);
    /// ```
    pub fn mag_sqr(&self) -> f64 {
        self.dot(self)
    }

    /// Magnitude of the vector.
    /// ```
    /// # use faust::phys_vec::PhysVec;
    /// let vec = PhysVec::from_cartesian(1., 2., 3.);
    /// assert_eq!(vec.mag(), 14_f64.sqrt());
    /// ```
    pub fn mag(&self) -> f64 {
        self.mag_sqr().sqrt()
    }

    /// Scale the vector by the given factor without changing the direction.
    /// ```
    /// # use faust::phys_vec::PhysVec;
    /// let mut vec = PhysVec::from_cartesian(1., 2., 3.);
    /// let mag = vec.mag();
    /// let phi = vec.phi_rad();
    /// let theta = vec.theta_rad();
    /// vec.scale(2.);
    ///
    /// assert!((vec.mag() - 2. * mag).abs() < 1e-10);
    /// assert!((vec.phi_rad() - phi).abs() < 1e-10);
    /// assert!((vec.theta_rad() - theta).abs() < 1e-10);
    /// ```
    pub fn scale(&mut self, factor: f64) -> &mut Self {
        (0..PhysVec::DIMS).for_each(|idx| {
            self.coordinate[idx] *= factor;
        });
        self
    }

    /// Assign the phi without changing the magnitude or theta.
    /// ```
    /// # use faust::phys_vec::PhysVec;
    /// # use std::f64::consts::PI;
    /// let mut vec = PhysVec::from_spherical(1., PI / 3., PI / 4.);
    /// let theta = vec.theta_rad();
    /// let mag = vec.mag();
    /// vec.set_phi_rad(PI / 2.);
    ///
    /// assert!((vec.phi_rad() - PI / 2.).abs() < 1e-10);
    /// assert!((vec.theta_rad() - theta).abs() < 1e-10);
    /// assert!((vec.mag() - mag).abs() < 1e-10);
    /// ```
    pub fn set_phi_rad(&mut self, phi: f64) -> &mut Self {
        let theta = self.theta_rad();
        let mag = self.mag();
        self.set_x(mag * theta.sin() * phi.cos());
        self.set_y(mag * theta.sin() * phi.sin());
        self.set_z(mag * theta.cos());
        self
    }

    /// Assign the theta without changing the magnitude or phi.
    /// ```
    /// # use faust::phys_vec::PhysVec;
    /// # use std::f64::consts::PI;
    /// let mut vec = PhysVec::from_spherical(1., PI / 3., PI / 4.);
    /// let phi = vec.phi_rad();
    /// let mag = vec.mag();
    /// vec.set_theta_rad(PI / 2.);
    ///
    /// assert!((vec.phi_rad() - phi).abs() < 1e-10);
    /// assert!((vec.theta_rad() - PI / 2.).abs() < 1e-10);
    /// assert!((vec.mag() - mag).abs() < 1e-10);
    /// ```
    pub fn set_theta_rad(&mut self, theta: f64) -> &mut Self {
        let phi = self.phi_rad();
        let mag = self.mag();

        self.set_x(mag * theta.sin() * phi.cos());
        self.set_y(mag * theta.sin() * phi.sin());
        self.set_z(mag * theta.cos());

        self
    }

    /// Assign the phi without changing the magnitude or theta.
    /// ```
    /// # use faust::phys_vec::PhysVec;
    /// # use std::f64::consts::PI;
    /// let mut vec = PhysVec::from_spherical(1., PI / 3., PI / 4.);
    /// let theta = vec.theta_deg();
    /// let mag = vec.mag();
    /// vec.set_phi_deg(90.);
    /// assert!((vec.phi_deg() - 90.).abs() < 1e-10);
    /// assert!((vec.theta_deg() - theta).abs() < 1e-10);
    /// assert!((vec.mag() - mag).abs() < 1e-10);
    /// ```
    pub fn set_phi_deg(&mut self, phi: f64) -> &mut Self {
        self.set_phi_rad(phi * PI / 180.)
    }

    /// Get the transverse component of the vector.
    /// ```
    /// # use faust::phys_vec::PhysVec;
    /// let vec = PhysVec::from_cartesian(1., 2., 3.);
    /// assert_eq!(vec.transverse(), 5_f64.sqrt());
    /// ```
    pub fn transverse(&self) -> f64 {
        (self.x().powi(2) + self.y().powi(2)).sqrt()
    }

    /// Get the beam-line component of the vector.
    /// ```
    /// # use faust::phys_vec::PhysVec;
    /// let vec = PhysVec::from_cartesian(1., 2., 3.);
    /// assert_eq!(vec.parallel(), 3.);
    /// ```
    pub fn parallel(&self) -> f64 {
        self.z()
    }

    /// Create a normalized version of the current vector
    /// ```
    /// # use faust::phys_vec::PhysVec;
    /// let vec = PhysVec::from_cartesian(1., 2., 3.);
    /// let norm = vec.as_normalized();
    /// assert!((norm.mag() - 1.).abs() < 1e-10);
    /// assert!((norm.phi_rad() - vec.phi_rad()).abs() < 1e-10);
    /// assert!((norm.theta_rad() - vec.theta_rad()).abs() < 1e-10);
    /// ```
    pub fn as_normalized(&self) -> Self {
        let magnitude = self.mag();

        Self::new(
            (0..self.dims())
                .map(|idx| self.at(idx) / magnitude)
                .collect::<Vec<f64>>()
                .try_into()
                .unwrap(),
        )
    }

    /// Create a scaled version of `self`.
    /// ```
    /// # use faust::phys_vec::PhysVec;
    /// let vec = PhysVec::from_cartesian(1., 2., 3.);
    /// let scaled = vec.as_scaled_by(2.);
    ///
    /// assert!((scaled.mag() - 2. * vec.mag()).abs() < 1e-10);
    /// assert!((scaled.phi_rad() - vec.phi_rad()).abs() < 1e-10);
    /// assert!((scaled.theta_rad() - vec.theta_rad()).abs() < 1e-10);
    /// ```
    pub fn as_scaled_by(&self, scaling_factor: f64) -> Self {
        Self::from_cartesian(
            self.x() * scaling_factor,
            self.y() * scaling_factor,
            self.z() * scaling_factor,
        )
    }

    /// Create a scaled version of `self` with the given magnitude.
    /// ```
    /// # use faust::phys_vec::PhysVec;
    /// let vec = PhysVec::from_cartesian(1., 2., 3.);
    /// let scaled = vec.as_normalized_to(2.);
    /// assert!((scaled.mag() - 2.).abs() < 1e-10);
    /// assert!((scaled.phi_rad() - vec.phi_rad()).abs() < 1e-10);
    /// assert!((scaled.theta_rad() - vec.theta_rad()).abs() < 1e-10);
    /// ```
    pub fn as_normalized_to(&self, new_magnitude: f64) -> Self {
        let old_magnitude = self.mag();
        let scaling_factor = new_magnitude / old_magnitude;
        Self::new(
            (0..self.dims())
                .map(|idx| self.at(idx) * scaling_factor)
                .collect::<Vec<f64>>()
                .try_into()
                .unwrap(),
        )
    }

    /// Normalize the vector in-place to unit magnitude.
    /// ```
    /// # use faust::phys_vec::PhysVec;
    /// let mut vec = PhysVec::from_cartesian(1., 2., 3.);
    /// vec.normalize();
    /// assert!((vec.mag() - 1.).abs() < 1e-10);
    /// ```
    pub fn normalize(&mut self) -> &mut Self {
        let norm = self.mag();

        for idx in 0..self.dims() {
            self.coordinate[idx] /= norm;
        }
        self
    }

    /// Calculate the inner angle between two vectors. Return value will be in radians between 0 and pi.
    /// ```
    /// # use faust::phys_vec::PhysVec;
    /// # use std::f64::consts::PI;
    /// let v_1 = PhysVec::from_spherical(1., 0., 0.);
    /// let v_2 = PhysVec::from_spherical(1., PI / 2., 0.);
    /// assert!((v_1.inner_angle_rad(&v_2) - PI / 2.).abs() < 1e-10);
    /// ```
    pub fn inner_angle_rad(&self, other: &Self) -> f64 {
        (self.dot(other) / (self.mag() * other.mag())).acos()
    }

    /// Calculate the inner angle between two vectors. Return value will be in degrees between 0 and 180.
    /// ```
    /// # use faust::phys_vec::PhysVec;
    /// # use std::f64::consts::PI;
    /// let v_1 = PhysVec::from_spherical(1., 0., 0.);
    /// let v_2 = PhysVec::from_spherical(1., PI / 2., 0.);
    /// assert!((v_1.inner_angle_deg(&v_2) - 90.).abs() < 1e-10);
    /// ```
    pub fn inner_angle_deg(&self, other: &Self) -> f64 {
        self.inner_angle_rad(other) * 180. / PI
    }

    pub fn rotated_around_vec(&self, other: &Self, angle_rad: f64) -> Self {
        let self_par = self.dot(other) / other.dot(other) * other;
        let self_perp = self.clone() - &self_par;

        let ortho = other.cross(&self_par);

        let x_1 = angle_rad.cos() / self_perp.mag(); 
        let x_2 = angle_rad.sin() / ortho.mag(); 

        let rotated_perp = self_perp.mag() * (x_1 * self_perp + x_2 * ortho);

        rotated_perp + self_par
    }
}

/// Add two physical vectors together and return the result as a new `PhysVec`
/// ```
/// # use faust::phys_vec::PhysVec;
/// # use std::f64::consts::PI;
/// let v_1 = PhysVec::from_cartesian(1., 2., 3.);
/// let v_2 = PhysVec::from_cartesian(4., 5., 6.);
/// let v_3 = PhysVec::from_cartesian(5., 7., 9.);
/// let v_sum = &v_1 + &v_2;
///
/// assert_eq!(v_sum, v_3);
///```
impl<'a, 'b> Add<&'b PhysVec> for &'a PhysVec {
    type Output = PhysVec;
    fn add(self, rhs: &'b PhysVec) -> Self::Output {
        let mut result = PhysVec::new([0., 0., 0.]);
        (0..3).for_each(|idx| {
            result.set_by_idx(idx, self.at(idx) + rhs.at(idx));
        });
        result
    }
}

/// Add two physical vectors together and return the result as a new `PhysVec`. Left hand side will
/// get dropped (rather it will get modified and returned, but it will feel like it's getting
/// dropped).
/// ```
/// # use faust::phys_vec::PhysVec;
/// # use std::f64::consts::PI;
/// let v_1 = PhysVec::from_cartesian(1., 2., 3.);
/// let v_2 = PhysVec::from_cartesian(4., 5., 6.);
/// let v_3 = PhysVec::from_cartesian(5., 7., 9.);
/// let v_sum = v_1 + &v_2;
///
/// assert_eq!(v_sum, v_3);
///```
impl<'a> Add<&'a PhysVec> for PhysVec {
    type Output = PhysVec;
    fn add(mut self, rhs: &'a PhysVec) -> Self::Output {
        self += rhs;
        self
    }
}

impl Add<PhysVec> for &PhysVec {
    type Output = PhysVec;
    fn add(self, rhs: PhysVec) -> Self::Output {
        let mut result = PhysVec::new([0., 0., 0.]);
        for idx in 0..3 {
            result.set_by_idx(idx, self.at(idx) + rhs.at(idx));
        }
        result
    }
}

impl Add<PhysVec> for PhysVec {
    type Output = PhysVec;
    fn add(self, rhs: PhysVec) -> Self::Output {
        let mut result = PhysVec::new([0., 0., 0.]);
        for idx in 0..3 {
            result.set_by_idx(idx, self.at(idx) + rhs.at(idx));
        }
        result
    }
}

impl<'a> AddAssign<&'a PhysVec> for PhysVec {
    fn add_assign(&mut self, rhs: &'a PhysVec) {
        for idx in 0..3 {
            self.set_by_idx(idx, self.at(idx) + rhs.at(idx));
        }
    }
}

impl<'a, 'b> Sub<&'b PhysVec> for &'a PhysVec {
    type Output = PhysVec;
    fn sub(self, rhs: &'b PhysVec) -> Self::Output {
        let mut result = PhysVec::new([0., 0., 0.]);
        for idx in 0..3 {
            result.set_by_idx(idx, self.at(idx) - rhs.at(idx));
        }
        result
    }
}
impl<'a, 'b> Sub<&'b PhysVec> for &'a mut PhysVec {
    type Output = PhysVec;
    fn sub(self, rhs: &'b PhysVec) -> Self::Output {
        let mut result = PhysVec::new([0., 0., 0.]);
        for idx in 0..3 {
            result.set_by_idx(idx, self.at(idx) - rhs.at(idx));
        }
        result
    }
}

impl<'a> Sub<&'a PhysVec> for PhysVec {
    type Output = PhysVec;
    fn sub(mut self, rhs: &'a PhysVec) -> Self::Output {
        self -= rhs;
        self
    }
}

impl Sub<PhysVec> for &PhysVec {
    type Output = PhysVec;
    fn sub(self, rhs: PhysVec) -> Self::Output {
        let mut result = PhysVec::new([0., 0., 0.]);
        (0..3).for_each(|idx| {
            result.set_by_idx(idx, self.at(idx) - rhs.at(idx));
        });
        result
    }
}

impl Sub<PhysVec> for &mut PhysVec {
    type Output = PhysVec;
    fn sub(self, rhs: PhysVec) -> Self::Output {
        let mut result = PhysVec::new([0., 0., 0.]);
        (0..3).for_each(|idx| {
            result.set_by_idx(idx, self.at(idx) - rhs.at(idx));
        });
        result
    }
}

impl Sub<PhysVec> for PhysVec {
    type Output = PhysVec;
    fn sub(self, rhs: PhysVec) -> Self::Output {
        let mut result = PhysVec::new([0., 0., 0.]);
        for idx in 0..3 {
            result.set_by_idx(idx, self.at(idx) - rhs.at(idx));
        }
        result
    }
}

impl<'a> SubAssign<&'a PhysVec> for PhysVec {
    fn sub_assign(&mut self, rhs: &'a PhysVec) {
        for idx in 0..3 {
            self.set_by_idx(idx, self.at(idx) - rhs.at(idx));
        }
    }
}

impl Mul<f64> for PhysVec {
    type Output = PhysVec;
    fn mul(self, rhs: f64) -> PhysVec {
        self.as_scaled_by(rhs)
    }
}

impl<'a> Mul<f64> for &'a PhysVec {
    type Output = PhysVec;
    fn mul(self, rhs: f64) -> PhysVec {
        self.as_scaled_by(rhs)
    }
}

impl<'a> Mul<&'a PhysVec> for f64 {
    type Output = PhysVec;
    fn mul(self, rhs: &PhysVec) -> PhysVec {
        rhs.as_scaled_by(self)
    }
}

impl Mul<PhysVec> for f64 {
    type Output = PhysVec;
    fn mul(self, rhs: PhysVec) -> PhysVec {
        rhs.as_scaled_by(self)
    }
}

impl Div<f64> for PhysVec {
    type Output = PhysVec;
    fn div(self, rhs: f64) -> PhysVec {
        self.as_scaled_by(1. / rhs)
    }
}

impl<'a> Div<f64> for &'a PhysVec {
    type Output = PhysVec;
    fn div(self, rhs: f64) -> PhysVec {
        self.as_scaled_by(1. / rhs)
    }
}

impl PartialEq for PhysVec {
    fn eq(&self, other: &Self) -> bool {
        (0..3).all(|idx| (self.at(idx) - other.at(idx)).abs() < 1e-10)
    }
}
