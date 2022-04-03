// [[file:../dimer.note::875f7ef9][875f7ef9]]
use super::*;

use crate::cg::CG;
use crate::fourier::FourierRotation;
// 875f7ef9 ends here

// [[file:../dimer.note::1b911cfd][1b911cfd]]
impl<'a> Dimer<'a> {
    /// Build internal dimer data (raw_dimer) from dimer setup, and return
    /// dimer state on new position.
    fn re_orientate(&mut self, n_unit: &DVector) -> Result<RotationState> {
        let dr = self.vars.distance;
        let r0 = self.center.clone();
        let r1 = &r0 + dr * n_unit;

        // FIXME: convergence test before evaluation
        if let Some(mut raw_dimer) = self.inner.take() {
            if !self.vars.use_extrapolated_forces {
                self.dynamics.set_position(r1.as_slice());
                let f1 = self.dynamics.get_force()?.to_vector();
                let s = f1.cosine_similarity(&raw_dimer.f1);
                debug!("similarity between extrapolated force and real force at R1: {}", s);
                raw_dimer.f1 = f1;
            }
            raw_dimer.r1 = r1;
            self.inner = Some(raw_dimer);
        } else {
            self.dynamics.set_position(r0.as_slice());
            let e0 = self.dynamics.get_energy()?;
            let f0 = self.dynamics.get_force()?.to_vector();
            self.dynamics.set_position(r1.as_slice());
            let f1 = self.dynamics.get_force()?.to_vector();
            // set active atoms weight matrix
            let raw_dimer = RawDimer {
                r0,
                r1,
                f0,
                f1,
                dr,
                e0,
                c0: None,
            };
            self.inner = Some(raw_dimer);
        }

        // update rotational direction perpendicular to dimer orientation
        let state = self.inner.as_ref().unwrap().extrapolate();
        let f_rot = state.rotational_force();
        assert!(f_rot.norm() > 0.0, "invalid rotational force: {:?}", &f_rot);

        Ok(state)
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
    fn next_rotation(&mut self, dimer_state: &RotationState, theta: &DVector, phi1: f64, phi_est: f64) -> Result<()> {
        let mut raw_dimer = self.inner.take().expect("no raw dimer");
        // get endpoint 1 (R1, F1) after trial rotation
        let r1_prime = raw_dimer.get_endpoint1_after_rotation(&self.orientation, &theta, phi1);
        self.dynamics.set_position(r1_prime.as_slice());
        let f1_prime = self.dynamics.get_force()?.to_vector();
        raw_dimer.r1 = r1_prime;
        raw_dimer.f1 = f1_prime.clone();
        let c1 = raw_dimer.extrapolate().curvature();

        // estimate optimal rotation using Fourier series
        let c0 = dimer_state.curvature();
        let c0d = dimer_state.curvature_derivative();
        let fourier_rot = FourierRotation::new(c0, c0d, phi1, c1);
        let (mut phi_min, mut curvature_min) = fourier_rot.optimal_rotation();

        // If find maximum curvature, the rotation angle has to be increased by 90°
        // FIXME: 如果f1是插值出来的话, 由于数据误差, 也可能出现下面的情况
        if curvature_min > c0 {
            info!("dimer fourier: found maximum curvature: {} > {}", curvature_min, c0);
            if !self.vars.use_extrapolated_forces {
                phi_min = phi_min + PI / 2.0;
                // update curvature in new position
                curvature_min = fourier_rot.curvature(phi_min);
            } else {
                // FIXME: 如果真能使curvature变小才加90度
                warn!("using extrapolated forces is not reliable for testing curvature minimum/maximum");
                let phi = phi_min + PI / 2.0;
                let c_min_prime = fourier_rot.curvature(phi);
                dbg!(c0, curvature_min, c_min_prime);
                if c_min_prime < c0 && c_min_prime < curvature_min {
                    phi_min = phi;
                    curvature_min = c_min_prime;
                }
            }
        }

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

        let f1_min = raw_dimer.estimate_force_after_rotation(phi_min, phi1, &f1_prime);
        // FIXME: think about dimer axis
        let r1_min = raw_dimer.get_endpoint1_after_rotation(&self.orientation, &theta, phi_min);
        raw_dimer.r1 = r1_min;
        raw_dimer.f1 = f1_min;
        raw_dimer.c0 = curvature_min.into();
        self.orientation = raw_dimer.extrapolate().curvature_mode().clone();
        self.inner = raw_dimer.into();

        Ok(())
    }
}
// 69cb7fbe ends here

// [[file:../dimer.note::45c98025][45c98025]]
impl<'a> Dimer<'a> {
    /// Rotate the dimer axis into the lowest curvature mode of the potential
    /// energy at the dimer center.
    pub fn get_optimal_rotation(&mut self, n_max_rot: usize) -> Result<(f64, DVector)> {
        let mut cg = CG::default();
        let mut tau = self.orientation.clone();
        let phi_tol = self.vars.min_rot_angle;
        for niter in 1..n_max_rot {
            info!("dimer rotation iteration {niter}");
            let state = self.re_orientate(&tau)?;
            // avoid trial rotation if estimated rotational angle `phi_est` is small
            // enough (Eq32 Heyden2005JCP)
            let phi_est = state.estimated_rotational_angle();
            let converged = check_dimer_rotation_convergence(phi_est, phi_tol);
            if converged {
                info!("Optimal dimer rotation found within {niter} iterations.");
                break;
            }
            // trial rotation: use a fixed angle or variable one. See p12 in Heyden2005JCP
            let phi1 = if self.vars.use_fixed_rot_angle {
                self.vars.trial_rot_angle
            } else {
                self.vars.trial_rot_angle.min(phi_est)
            };
            // rotate the dimer in optimal direction
            let theta = self.get_rotational_direction(state.rotational_force(), &mut cg);
            self.next_rotation(&state, &theta, phi1, phi_est)?;
            if niter >= n_max_rot {
                info!("Max allowed iterations reached");
                break;
            }
        }
        // FIXME: rewrite
        let phi = self.orientation.cosine_similarity(&tau).acos();
        info!("Total rotational angle = {:.2}°", phi.to_degrees());

        // Total rotation angle during rotation steps
        // FIXME: rewrite
        let curvature_min = self.inner.as_ref().unwrap().c0.unwrap();
        let tau_min = self.orientation.clone();

        Ok((curvature_min, tau_min))
    }
}
// 45c98025 ends here
