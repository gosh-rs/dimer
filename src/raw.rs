// [[file:../dimer.note::d5c73cde][d5c73cde]]
fn get_dimer_endpoints(r0: &DVector, dr: f64, n: &DVector) -> [DVector; 2] {
    let r1 = r0 + dr * n;
    let r2 = r0 - dr * n;
    [r1, r2]
}
// d5c73cde ends here

// [[file:../dimer.note::bfde551f][bfde551f]]
fn compute_dimer_endpoint2_force(f0: &DVector, f1: &DVector) -> DVector {
    2.0 * f0 - f1
}
// bfde551f ends here

// [[file:../dimer.note::68051c57][68051c57]]
fn compute_rotational_force(f0: &DVector, f1: &DVector, dr: f64) -> DVector {
    (f1 - f0) / dr
}
// 68051c57 ends here

// [[file:../dimer.note::746f3305][746f3305]]
fn compute_dimer_axis(r0: &DVector, r1: &DVector) -> DVector {
    (r1 - r0).normalize()
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
    pub r0: DVector,
    /// Postions of image 1 in dimer
    pub r1: DVector,
    /// Forces of image 0 in dimer
    pub f0: DVector,
    /// Forces of image 1 in dimer
    pub f1: DVector,
    /// Curvature of dimer to be determined
    pub c0: Option<f64>,
    /// Real energy of image R0
    pub e0: f64,
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
    pub fn rotational_force(&self) -> &DVector {
        &self.fr
    }

    /// Return dimer axis direction
    pub fn dimer_axis(&self) -> &DVector {
        &self.n
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
