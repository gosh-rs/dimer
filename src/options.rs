// [[file:../dimer.note::c38894e0][c38894e0]]
use super::*;

/// Options for tuning DIMER algorithm from user side
#[derive(Deserialize, Debug, Clone, Serialize)]
#[serde(default)]
pub struct UserOptions {
    /// force component criteria for convergence test
    pub fmax: f64,

    /// dimer distance between image 0 and image 1
    pub distance: f64,

    /// Trial angle for the finite difference estimate of the rotational angle
    /// in radians.
    pub trial_rot_angle: f64,

    /// Use a fixed angle for trial rotation step.
    pub use_fixed_rot_angle: bool,

    /// The minimum rotational angle for skipping trial rotation step.
    pub min_rot_angle: f64,

    /// Maximum number of rotation steps allowed in iteratons.
    pub max_num_rot: usize,

    /// Use extrapolated force on the new dimer endpoint R1 to reduce one
    /// evulation per dimer rotation. Kastner2008JCP
    pub use_extrapolated_force: bool,

    /// Use Conjugate gradient algorithm to determine the rotation plane,
    /// instead of simple steepest descent direction.
    pub use_cg_rot: bool,
}

impl Default for UserOptions {
    fn default() -> Self {
        Self {
            fmax: 0.1,
            distance: 1E-3,
            min_rot_angle: 5f64.to_radians(),
            trial_rot_angle: PI / 4.0,
            use_fixed_rot_angle: true,
            max_num_rot: 5,
            use_extrapolated_force: false,
            use_cg_rot: true,
        }
    }
}
// c38894e0 ends here
