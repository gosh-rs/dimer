// [[file:../dimer.note::c492aa76][c492aa76]]
//! Implementation of the Conjugate Gradient (CG) optimization algorithm
//!
//! # References
//! - <https://en.wikipedia.org/wiki/Nonlinear_conjugate_gradient_method>
//! - <https://github.com/siesta-project/flos/blob/master/flos/optima/cg.lua>

use super::*;
use vecfx::*;
// c492aa76 ends here

// [[file:../dimer.note::2f6d095d][2f6d095d]]
/// The beta value to determine the step of the steepest descent direction
#[derive(Debug, Clone)]
pub enum BetaKind {
    /// Polak-Ribiere
    PR,
    /// Fletcher-Reeves
    FR,
    /// Hestenes-Stiefel
    HS,
    /// Dai-Yuan
    DY,
}

impl Default for BetaKind {
    fn default() -> Self {
        BetaKind::PR
    }
}

/// Method for determining when to restart a CG optimization
#[derive(Debug, Clone)]
pub enum RestartMethod {
    /// When the scalar-projection of the two previous gradients is above 0.2
    Powell,
    /// When `beta < 0`, CG restarts the conjugate gradient
    Negative,
}

impl Default for RestartMethod {
    fn default() -> Self {
        RestartMethod::Powell
    }
}

/// History data during conjugate gradient optimization
#[derive(Debug, Clone)]
struct ConjugateGradientState {
    /// The forces
    forces: DVector,
    /// The conjugate direction
    conjct: DVector,
}

#[derive(Debug, Clone)]
pub struct ConjugateGradient {
    /// The state of previous step
    state: Option<ConjugateGradientState>,

    /// Method of calculating the beta constant
    beta: BetaKind,

    /// Damping factor for creating a smooth CG restart
    beta_damping: f64,

    /// The method to restart CG optimization
    restart: RestartMethod,
}

impl Default for ConjugateGradient {
    fn default() -> Self {
        ConjugateGradient {
            state: None,
            beta: BetaKind::default(),
            restart: RestartMethod::default(),
            beta_damping: 0.8,
        }
    }
}

pub type CG = ConjugateGradient;
// 2f6d095d ends here

// [[file:../dimer.note::1e041d59][1e041d59]]
impl ConjugateGradient {
    /// Return the new conjugate direction
    pub fn propagate(&mut self, forces: &DVector) -> DVector {
        self.propagate_dimer(forces, None)
    }

    /// Return the new conjugate direction
    pub fn propagate_dimer(&mut self, forces: &DVector, dimer_orientation: Option<&DVector>) -> DVector {
        // FIXME: ad hoc hacking for dimer rotation
        let forces = dimer_orientation.map_or(forces.clone(), |tau| forces.vector_rejection(tau));

        // take out previous state, or update it with current forces
        let state = self.state.get_or_insert(ConjugateGradientState {
            forces: forces.clone(),
            conjct: forces.clone(),
        });

        // udpate beta
        let beta = self.beta.update(&forces, &state);

        // restart
        let beta = self.beta_damping
            * match self.restart {
                RestartMethod::Negative => beta.max(0.0),
                // Here we check whether the gradient of the current iteration has
                // "lost" orthogonality to the previous iteration
                RestartMethod::Powell => {
                    let n = forces.norm_squared();
                    let m = forces.dot(&state.forces);
                    if n / m >= 0.2 {
                        0.0
                    } else {
                        beta
                    }
                }
                _ => unimplemented!(),
            };

        // Now calculate the new steepest descent direction
        let disp = &forces + beta * &state.conjct;

        // FIXME: ad hoc hacking for dimer rotation
        let disp = dimer_orientation.map_or(disp.clone(), |tau| disp.vector_rejection(tau));

        // save state
        state.forces = forces.clone();
        state.conjct = disp.clone();

        disp
    }
}

impl BetaKind {
    /// Return the beta value to determine the step of the steepest descent direction
    /// # Parameters
    /// - forces: current forces
    /// - state : stored data in previous step
    fn update(&self, forces: &DVector, state: &ConjugateGradientState) -> f64 {
        let forces_this = forces;
        let forces_prev = &state.forces;
        let conjct_prev = &state.conjct;

        match self {
            BetaKind::PR => forces_this.dot(&(forces_this - forces_prev)) / forces_prev.norm_squared(),
            BetaKind::FR => forces_this.norm_squared() / forces_prev.norm_squared(),
            BetaKind::HS => {
                let d = forces_this - forces_prev;
                -forces_this.dot(&d) / conjct_prev.dot(&d)
            }
            BetaKind::DY => {
                let d = forces_this - forces_prev;
                forces_this.norm_squared() / conjct_prev.dot(&d)
            }
            _ => {
                error!("unkown beta parameter scheme!");
                unimplemented!()
            }
        }
    }
}
// 1e041d59 ends here
