// [[file:../dimer.note::875f7ef9][875f7ef9]]
use super::*;

use crate::cg::CG;
// 875f7ef9 ends here

// [[file:../dimer.note::1b911cfd][1b911cfd]]
impl<'a> Dimer<'a> {
    /// Rebuild `RawDimer` from start (updating center and endpoint 1)
    fn reinitialize(&mut self) -> Result<(RawDimer, f64)> {
        let dr = self.vars.distance;
        let r0 = self.center.clone();
        let r1 = &r0 + dr * &self.orientation;

        self.dynamics.set_position(r0.as_slice());
        let f0 = self.dynamics.get_force()?.to_vector();
        let e0 = self.dynamics.get_energy()?;
        self.dynamics.set_position(r1.as_slice());
        let f1 = self.dynamics.get_force()?.to_vector();
        let raw_dimer = RawDimer { r0, r1, f0, f1 };
        Ok((raw_dimer, e0))
    }
}
// 1b911cfd ends here

// [[file:../dimer.note::5bff1ad1][5bff1ad1]]
impl<'a> Dimer<'a> {
    // update rotational direction perpendicular to dimer orientation.
    // use conjugate-gradient or steepest descent to determinte rotational direction
    fn get_rotational_direction(&mut self, f_rot: &DVector, cg: &mut crate::cg::CG) -> DVector {
        let tau = &self.orientation;
        if self.vars.use_cg_rot {
            cg.propagate_dimer(f_rot, Some(tau)).normalize()
        } else {
            f_rot.vector_rejection(tau).normalize()
        }
    }
}
// 5bff1ad1 ends here

// [[file:../dimer.note::04c155c6][04c155c6]]
fn check_dimer_rotation_convergence(phi_est: f64, phi_tol: f64) -> bool {
    if phi_est.abs() < phi_tol {
        info!(
            "rotational angle is small enough: |{:.2}|° < {:.2}°",
            phi_est.to_degrees(),
            phi_tol.to_degrees()
        );
        return true;
    }

    false
}
// 04c155c6 ends here

// [[file:../dimer.note::69cb7fbe][69cb7fbe]]
impl<'a> Dimer<'a> {
    /// Estimate optimal rotation by rotating the dimer in direction `theta`
    /// with trial angle `phi1` using Fourier transform.
    ///
    /// # Parameters
    ///
    /// * raw_dimer: RawDimer to be rotated
    /// * theta: rotation direction
    /// * phi1: trial rotation angle
    ///
    fn rotate_dimer_within(&mut self, raw_dimer: &mut RawDimer, theta: &DVector, phi1: f64, phi_est: f64) -> Result<f64> {
        // get endpoint 1 (R1, F1) after trial rotation
        let r1_prime = raw_dimer.get_endpoint1_after_rotation(&self.orientation, &theta, phi1);
        self.dynamics.set_position(r1_prime.as_slice());
        let f1_prime = self.dynamics.get_force()?.to_vector();
        let fourier_state = raw_dimer.fourier_rotate(r1_prime, f1_prime, phi1, theta, self.vars.use_extrapolated_force);
        let phi_min = fourier_state.phi_min;
        let curvature_min = fourier_state.curvature_min;
        info!(
            "{:^15}{:^15}{:^15}{:^15}",
            "phi_est/deg", "phi_trial/deg", "phi_min/deg", "c_min_est"
        );
        info!(
            "{:^15.2}{:^15.2}{:^15.2}{:^-15.3}",
            phi_est.to_degrees(),
            phi1.to_degrees(),
            phi_min.to_degrees(),
            curvature_min
        );

        raw_dimer.r1 = fourier_state.r1_min;
        raw_dimer.f1 = fourier_state.f1_min;

        Ok(curvature_min)
    }
}
// 69cb7fbe ends here

// [[file:../dimer.note::45c98025][45c98025]]
/// Represents the results obtained in rotation step
#[derive(Debug, Clone)]
pub struct RotationOutput {
    pub raw_dimer: RawDimer,
    /// The potential energy evaluated at dimer center
    pub energy: f64,
    /// The optimized curvature
    pub curvature_min: f64,
    /// The number of iterations used in rotation step
    pub n_iterations: usize,
}

/// The part for DIMER rotation
impl<'a> Dimer<'a> {
    /// Rotate the dimer axis into the lowest curvature mode of the potential
    /// energy at the dimer center estimated in Fourier series.
    pub(crate) fn next_rotation_step(&mut self, n_max_rot: usize) -> Result<RotationOutput> {
        let tau_ini = self.orientation.clone();
        let phi_tol = self.vars.min_rot_angle;

        let mut cg = CG::default();
        let (mut raw_dimer, e0) = self.reinitialize()?;
        // save the state before trial rotation
        let mut state = raw_dimer.extrapolate();
        let mut curvature_min = state.curvature();
        let mut niter = 0;
        loop {
            niter += 1;
            info!("dimer rotation iteration {niter}");
            // avoid trial rotation if estimated rotational angle `phi_est` is small
            // enough (Eq32 Heyden2005JCP)
            let phi_est = state.estimated_rotational_angle();
            let converged = check_dimer_rotation_convergence(phi_est, phi_tol);
            match (converged, niter >= n_max_rot) {
                (true, _) => {
                    info!("Optimal dimer rotation found within {niter} iterations.");
                    break;
                }
                (false, true) => {
                    warn!("Max allowed iterations {n_max_rot} reached, but dimer rotation not converged yet.");
                    break;
                }
                (false, false) => {}
            }

            // trial rotation: use a fixed angle or variable one. See p12 in Heyden2005JCP
            let phi1 = if self.vars.use_fixed_rot_angle {
                self.vars.trial_rot_angle
            } else {
                self.vars.trial_rot_angle.min(phi_est)
            };
            // rotate `raw_dimer` in optimal direction with a angle leading to lowest curvature
            let f_rot = state.rotational_force();
            assert!(f_rot.norm() > 0.0, "invalid rotational force: {:?}", &f_rot);
            let theta = self.get_rotational_direction(f_rot, &mut cg);
            let curvature_min_est = self.rotate_dimer_within(&mut raw_dimer, &theta, phi1, phi_est)?;
            // Update extrapolated force of endpint `1` if necessary
            if !self.vars.use_extrapolated_force {
                self.dynamics.set_position(raw_dimer.r1.as_slice());
                let f1 = self.dynamics.get_force()?.to_vector();
                let s = f1.cosine_similarity(&raw_dimer.f1);
                debug!("similarity between extrapolated force and real force of endpoint 1: {}", s);
                raw_dimer.f1 = f1;
            }
            // Update dimer state after rotation
            state = raw_dimer.extrapolate();
            // Update current dimer orientation, important for translation step
            self.orientation = state.curvature_mode().clone();
            // Recalculate curvature. If we do not use extrapolated f1, the
            // curvature_min should be updated with more accurate number
            curvature_min = state.curvature();
            debug!("real curvature vs estimated curvature: {curvature_min} vs. {curvature_min_est}");
        }
        // Total rotation angle during rotation steps
        let phi = self.orientation.cosine_similarity(&tau_ini).acos().to_degrees();
        info!("Total rotational angle = {phi:.2}°; optimized curvature = {curvature_min}");

        let out = RotationOutput {
            raw_dimer,
            curvature_min,
            energy: e0,
            n_iterations: niter,
        };
        Ok(out)
    }
}
// 45c98025 ends here
