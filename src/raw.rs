// [[file:../dimer.note::d5c73cde][d5c73cde]]
fn get_dimer_endpoints(r0: &[f64], dr: f64, n: &[f64]) -> [DVector; 2] {
    let r1 = r0.as_vector_slice() + dr * n.as_vector_slice();
    let r2 = r0.as_vector_slice() - dr * n.as_vector_slice();
    [r1, r2]
}
// d5c73cde ends here

// [[file:../dimer.note::bfde551f][bfde551f]]
fn compute_dimer_endpoint2_force(f0: &[f64], f1: &[f64]) -> DVector {
    2.0 * f0.as_vector_slice() - f1.as_vector_slice()
}
// bfde551f ends here

// [[file:../dimer.note::*docs][docs:3]]
fn compute_rotational_force(f0: &[f64], f1: &[f64], dr: f64) -> DVector {
    (f1.as_vector_slice() - f0.as_vector_slice()) / dr
}
// docs:3 ends here

// [[file:../dimer.note::746f3305][746f3305]]
fn compute_dimer_axis(r0: &[f64], r1: &[f64]) -> DVector {
    (r1.as_vector_slice() - r0.as_vector_slice()).normalize()
}
// 746f3305 ends here

// [[file:../dimer.note::31639bf8][31639bf8]]
fn compute_rotational_direction(f_rot: &DVector, n_unit: &DVector) -> DVector {
    let f_v = f_rot - f_rot.dot(n_unit) * n_unit;
    f_v.normalize()
}
// 31639bf8 ends here

// [[file:../dimer.note::b8bd7012][b8bd7012]]
fn compute_dimer_curvature(f_rot: &DVector, n_unit: &DVector) -> f64 {
    -f_rot.dot(n_unit)
}
// b8bd7012 ends here

// [[file:../dimer.note::add08c3f][add08c3f]]
fn compute_dimer_curvature_derivative(f_rot: &DVector, theta_unit: &DVector) -> f64 {
    -2.0 * f_rot.dot(theta_unit)
}
// add08c3f ends here

// [[file:../dimer.note::de6e084c][de6e084c]]
use super::*;
// de6e084c ends here

// [[file:../dimer.note::4648b13c][4648b13c]]
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct RawDimer {
    /// Distance between dimer central image 0 and endpoint image 1
    pub dr: f64,
    /// Postions of image 0 in dimer
    pub r0: Vec<f64>,
    /// Postions of image 1 in dimer
    pub r1: Vec<f64>,
    /// Forces of image 0 in dimer
    pub f0: Vec<f64>,
    /// Forces of image 1 in dimer
    pub f1: Vec<f64>,
}

impl RawDimer {
    /// Return position of dimer center
    pub fn center(&self) -> &[f64] {
        &self.r0
    }
}
// 4648b13c ends here

// [[file:../dimer.note::62d61ee4][62d61ee4]]
/// Represents internal state of a dimer
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct DimerState {
    /// dimer mode curvature
    cx: f64,
    /// Rotational force
    fr: DVector,
    /// Dimer axis direction
    n: DVector,
}

impl DimerState {
    /// Dimer curvature
    pub fn curvature(&self) -> f64 {
        self.cx
    }

    /// Dimer curvature derivative
    pub fn curvature_derivative(&self) -> f64 {
        let theta = compute_rotational_direction(&self.fr, &self.n);
        compute_dimer_curvature_derivative(&self.fr, &theta)
    }

    /// The rotational force felt at dimer center
    pub fn rotational_force(&self) -> &[f64] {
        self.fr.as_slice()
    }

    /// Return dimer axis direction
    pub fn dimer_axis(&self) -> &[f64] {
        self.n.as_slice()
    }
}
// 62d61ee4 ends here

// [[file:../dimer.note::02a5a92f][02a5a92f]]
impl RawDimer {
    /// Extrapolate other missing dimer properties
    pub fn extrapolate(&self) -> DimerState {
        let fr = compute_rotational_force(&self.f0, &self.f1, self.dr);
        let n = compute_dimer_axis(&self.r0, &self.r1);
        let cx = compute_dimer_curvature(&fr, &n);
        DimerState { fr, n, cx }
    }
}
// 02a5a92f ends here
