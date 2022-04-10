// [[file:../dimer.note::875f7ef9][875f7ef9]]
use super::*;

use crate::cg::CG;
// 875f7ef9 ends here

// [[file:../dimer.note::1b911cfd][1b911cfd]]
impl<'a> Dimer<'a> {
    /// Rebuild `RawDimer` from start (updating center and endpoint 1)
    fn reinitialize(&mut self) -> Result<RawDimer> {
        let dr = self.vars.distance;
        let r0 = self.center.clone();
        let r1 = &r0 + dr * &self.orientation;

        self.dynamics.set_position(r0.as_slice());
        let f0 = self.dynamics.get_force()?.to_vector();
        // FIXME: about dimer center point energy
        // let e0 = self.dynamics.get_energy()?;
        self.dynamics.set_position(r1.as_slice());
        let f1 = self.dynamics.get_force()?.to_vector();
        let raw_dimer = RawDimer { r0, r1, f0, f1 };
        Ok(raw_dimer)
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
    /// * dimer_state: current state before trial rotation
    /// * theta: rotation direction
    /// * phi1: trial rotation angle
    ///
    fn next_rotation(&mut self, raw_dimer: &mut RawDimer, theta: &DVector, phi1: f64, phi_est: f64) -> Result<f64> {
        // get endpoint 1 (R1, F1) after trial rotation
        let r1_prime = raw_dimer.get_endpoint1_after_rotation(&self.orientation, &theta, phi1);
        self.dynamics.set_position(r1_prime.as_slice());
        let f1_prime = self.dynamics.get_force()?.to_vector();
        let fourier_state = raw_dimer.fourier_rotate(r1_prime, f1_prime, phi1, theta, self.vars.use_extrapolated_force);
        let phi_min = fourier_state.phi_min;
        let curvature_min = fourier_state.curvature_min;
        info!(
            "{:^15}{:^15}{:^15}{:^15}",
            "phi_est/deg", "phi_trial/deg", "phi_min/deg", "c_min"
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
/// The part for DIMER rotation
impl<'a> Dimer<'a> {
    /// Rotate the dimer axis into the lowest curvature mode of the potential
    /// energy at the dimer center. Return optimized curvature value on success.
    pub fn get_optimal_rotation(&mut self, n_max_rot: usize) -> Result<(RawDimer, f64)> {
        let phi_tol = self.vars.min_rot_angle;

        let mut cg = CG::default();
        let mut raw_dimer = self.reinitialize()?;
        // save the state before trial rotation
        let mut state = raw_dimer.extrapolate();
        let mut curvature_min = state.curvature();
        for niter in 1.. {
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
                    warn!("Max allowed iterations reached, but dimer rotation not converged yet.");
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
            self.next_rotation(&mut raw_dimer, &theta, phi1, phi_est)?;
            // update extrapolated force of endpint `1` if necessary
            if !self.vars.use_extrapolated_force {
                self.dynamics.set_position(raw_dimer.r1.as_slice());
                let f1 = self.dynamics.get_force()?.to_vector();
                let s = f1.cosine_similarity(&raw_dimer.f1);
                debug!("similarity between extrapolated force and real force at R1: {}", s);
                raw_dimer.f1 = f1;
            }
            // update dimer state after rotation
            state = raw_dimer.extrapolate();
        }
        curvature_min = state.curvature();
        let tau_min = state.curvature_mode();
        // Total rotation angle during rotation steps
        let phi = self.orientation.cosine_similarity(tau_min).acos();
        info!("Total rotational angle = {:.2}°; c_min = {curvature_min}", phi.to_degrees());
        // update current dimer orientation, important for translation step
        self.orientation = tau_min.clone();

        Ok((raw_dimer, curvature_min))
    }
}
// 45c98025 ends here
