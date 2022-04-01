// [[file:../dimer.note::87727f93][87727f93]]
use super::*;
// 87727f93 ends here

// [[file:../dimer.note::b4042e0e][b4042e0e]]
/// Return coordinates of the dimer endpoint `R1` after rotation with a trial angle `phi`.
pub fn get_dimer_trial_rotation_endpoint(r0: &[f64], n_unit: &[f64], t_unit: &[f64], phi: f64, dr: f64) -> DVector {
    r0.as_vector_slice() + dr * phi.cos() * n_unit.as_vector_slice() + dr * phi.sin() * t_unit.as_vector_slice()
}
// b4042e0e ends here

// [[file:../dimer.note::fbf1f130][fbf1f130]]
fn estimate_rotational_angle(c0: f64, c0d: f64) -> f64 {
    0.5 * (-0.5 * c0d / c0.abs()).atan()
}
// fbf1f130 ends here

// [[file:../dimer.note::*algo][algo:2]]
/// # Return a tuple of Fourier constants: (a0, a1, b1)
fn get_fourier_series_constants(c0: f64, c0d: f64, c1: f64, phi1: f64) -> (f64, f64, f64) {
    let phi = 2.0 * phi1;
    let b1 = 0.5 * c0d;
    let a1 = (c0 - c1 + b1 * phi.sin()) / (1.0 - phi.cos());
    let a0 = 2.0 * (c0 - a1);

    (a0, a1, b1)
}
// algo:2 ends here

// [[file:../dimer.note::*algo][algo:3]]
fn get_phi_min_by_fourier_series(a1: f64, b1: f64) -> f64 {
    0.5 * (b1 / a1).atan()
}
// algo:3 ends here

// [[file:../dimer.note::*algo][algo:4]]
/// Eq 24 in Heyden2005JCP
///
/// * a0, a1: Fourier series constants
///
/// * phi: rotation angle
fn get_curvature_by_fourier_series(a0: f64, a1: f64, b1: f64, phi: f64) -> f64 {
    a0 / 2.0 + a1 * (2.0 * phi).cos() + b1 * (2.0 * phi).sin()
}
// algo:4 ends here

// [[file:../dimer.note::b2282332][b2282332]]
pub fn get_extrapolated_forces(phi1: f64, phi_min: f64, f1: &[f64], f1_prime: &[f64]) -> DVector {
    let f1 = f1.as_vector_slice();
    let f1_prime = f1_prime.as_vector_slice();
    (phi1 - phi_min).sin() / phi1.sin() * f1
        + phi_min.sin() / phi1.sin() * f1_prime
        + (1.0 - phi_min.cos() - phi_min.sin() * (0.5 * phi1).tan()) * f1
}
// b2282332 ends here

// [[file:../dimer.note::c701d372][c701d372]]
impl DimerState {
    /// Return estimated angle for trial rotation
    pub fn estimated_rotational_angle(&self) -> f64 {
        let c0 = self.curvature();
        let c0d = self.curvature_derivative();
        estimate_rotational_angle(c0, c0d)
    }
}

pub struct FourierRotation {
    // Fourier constants
    a0: f64,
    // Fourier constants
    a1: f64,
    // Fourier constants
    b1: f64,
    // Initial curvature
    c0: f64,
}

impl FourierRotation {
    /// # Parameters
    ///
    /// * c0: Curvature when phi = 0
    /// * c0d: Curvature derivative when phi = 0
    /// * phi1: trial rotation angle
    /// * c1: Curvature when phi = phi1
    pub fn new(c0: f64, c0d: f64, phi1: f64, c1: f64) -> Self {
        let (a0, a1, b1) = get_fourier_series_constants(c0, c0d, c1, phi1);
        Self { a0, a1, b1, c0 }
    }

    /// Interpolated curvature when rotation angle is `phi`
    pub fn curvature(&self, phi: f64) -> f64 {
        get_curvature_by_fourier_series(self.a0, self.a1, self.b1, phi)
    }

    /// Return optimal rotation angle and the minimum curvature.
    pub fn optimal_rotation(&self) -> (f64, f64) {
        let mut phi_min = get_phi_min_by_fourier_series(self.a1, self.b1);

        (phi_min, self.curvature(phi_min))
    }
}
// c701d372 ends here