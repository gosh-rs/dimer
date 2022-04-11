// [[file:../dimer.note::87727f93][87727f93]]
use super::*;
// 87727f93 ends here

// [[file:../dimer.note::b4042e0e][b4042e0e]]
/// Return coordinates of the dimer endpoint `R1` after rotation with a trial angle `phi`.
fn rotate_dimer_endpoint1(r0: &DVector, n_unit: &DVector, t_unit: &DVector, phi: f64, dr: f64) -> DVector {
    r0 + dr * phi.cos() * n_unit + dr * phi.sin() * t_unit
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
fn get_extrapolated_force(phi1: f64, phi_min: f64, f1: &DVector, f1_prime: &DVector) -> DVector {
    (phi1 - phi_min).sin() / phi1.sin() * f1
        + phi_min.sin() / phi1.sin() * f1_prime
        + (1.0 - phi_min.cos() - phi_min.sin() * (0.5 * phi1).tan()) * f1
}
// b2282332 ends here

// [[file:../dimer.note::4e1ae80e][4e1ae80e]]
struct FourierRotation {
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
    /// * c0: curvature when phi = 0
    /// * c0d: curvature derivative when phi = 0
    /// * phi1: trial rotation angle
    /// * c1: curvature when phi = phi1
    fn new(c0: f64, c0d: f64, phi1: f64, c1: f64) -> Self {
        let (a0, a1, b1) = get_fourier_series_constants(c0, c0d, c1, phi1);
        Self { a0, a1, b1, c0 }
    }

    /// Interpolated curvature when rotation angle is `phi`
    fn curvature(&self, phi: f64) -> f64 {
        get_curvature_by_fourier_series(self.a0, self.a1, self.b1, phi)
    }

    /// Return optimal rotation angle and the minimum curvature.
    fn optimal_rotation(&self) -> (f64, f64) {
        let mut phi_min = get_phi_min_by_fourier_series(self.a1, self.b1);

        (phi_min, self.curvature(phi_min))
    }
}
// 4e1ae80e ends here

// [[file:../dimer.note::c701d372][c701d372]]
/// State variables after Fourier rotation
#[derive(Debug, Clone)]
pub struct FourierState {
    /// The minium curvature estimated using Fourier series
    pub curvature_min: f64,
    /// The rotation angle corresponding to minimum curvature mode
    pub phi_min: f64,
    /// The position of endpoint `1` when optimal rotation applied
    pub r1_min: DVector,
    /// Extrapolated force of endpoint `1` when optimal rotation applied
    pub f1_min: DVector,
}

impl RotationState {
    /// Return estimated angle for trial rotation
    pub fn estimated_rotational_angle(&self) -> f64 {
        let c0 = self.curvature();
        let c0d = self.curvature_derivative();
        estimate_rotational_angle(c0, c0d)
    }
}

impl RawDimer {
    /// Return coordinates of the dimer endpoint `R1` after rotation in
    /// direction `theta` by angle `phi` in the plane spanned by `tau` and
    /// `theta`.
    pub fn get_endpoint1_after_rotation(&self, tau: &DVector, theta: &DVector, phi: f64) -> DVector {
        let dr = (&self.r1 - &self.r0).norm();
        rotate_dimer_endpoint1(&self.r0, &tau, theta, phi, dr)
    }

    /// Estimate optimal rotation using Fourier series
    ///
    /// # Parameters
    ///
    /// * r1_prime: position of endpoint `1` in trial rotation
    /// * f1_prime: force of endpoint `1` in trial rotation
    /// * phi1: trial rotation angle applied
    /// * theta: rotation direction
    /// * extrapolated_force: if f1 is extrapolated
    pub fn fourier_rotate(
        &mut self,
        r1_prime: DVector,
        f1_prime: DVector,
        phi1: f64,
        theta: &DVector,
        extrapolated_force: bool,
    ) -> FourierState {
        // get dimer state before trial rotation
        let n0 = self.dimer_axis();
        let state = self.extrapolate();
        let c0 = state.curvature();
        let c0d = state.curvature_derivative();
        let f1 = self.f1.clone();

        // trial rotation: update dimer with new endpoint1
        self.r1 = r1_prime;
        self.f1 = f1_prime;
        let c1 = self.extrapolate().curvature();

        let fourier_rot = FourierRotation::new(c0, c0d, phi1, c1);
        let (mut phi_min, mut curvature_min) = fourier_rot.optimal_rotation();

        // if find maximum curvature, the rotation angle has to be increased by 90°.
        // when F1 is extrapolated, this also could be triggered due to numberic noise
        if curvature_min > c0 {
            info!("found maximum curvature: {curvature_min} > {c0}");
            let phi = phi_min + PI / 2.0;
            if extrapolated_force {
                warn!("using extrapolated force is not reliable for testing curvature minimum/maximum");
                phi_min = phi;
                // update curvature in new position
                curvature_min = fourier_rot.curvature(phi_min);
            } else {
                // update rotation angle only when curvature can be really lowered
                let c_min_prime = fourier_rot.curvature(phi);
                log_dbg!(c0, curvature_min, c_min_prime);
                if c_min_prime < c0 && c_min_prime < curvature_min {
                    phi_min = phi;
                    curvature_min = c_min_prime;
                }
            }
        }

        // estimate force on new endpoint1
        let r1_min = self.get_endpoint1_after_rotation(&n0, &theta, phi_min);
        let f1_min = get_extrapolated_force(phi1, phi_min, &f1, &self.f1);

        FourierState {
            r1_min,
            f1_min,
            phi_min,
            curvature_min,
        }
    }
}
// c701d372 ends here
