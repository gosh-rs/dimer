// [[file:../dimer.note::875f7ef9][875f7ef9]]
use super::*;

use fourier::FourierRotation;
// 875f7ef9 ends here

// [[file:../dimer.note::1b911cfd][1b911cfd]]
impl<'a> Dimer<'a> {
    /// Build internal dimer data (raw_dimer) from dimer setup, return
    /// rotational forces on new position.
    pub fn re_orientate(&mut self) -> Result<DVector> {
        let n_unit = &self.orientation;
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
        // FIXME: avoid clone
        let f_rot = state.rotational_force().clone();
        assert!(f_rot.norm() > 0.0, "invalid rotational force: {:?}", &f_rot);

        Ok(f_rot)
    }
}
// 1b911cfd ends here

// [[file:../dimer.note::5bff1ad1][5bff1ad1]]
impl<'a> Dimer<'a> {
    // update rotational direction perpendicular to dimer orientation.
    // use conjugate-gradient or steepest descent to determinte rotational direction
    pub fn get_rotational_direction(&mut self, f_rot: &DVector, cg: &mut crate::cg::CG) -> DVector {
        let tau = &self.orientation;
        if self.vars.use_cg_rot {
            cg.propagate_dimer(f_rot, Some(tau)).normalize()
        } else {
            f_rot.vector_rejection(tau).normalize()
        }
    }
}
// 5bff1ad1 ends here

// [[file:../dimer.note::69cb7fbe][69cb7fbe]]
impl<'a> Dimer<'a> {
    /// Return Ok(true) if converged
    pub fn next_rotation_step(&mut self, theta: &DVector) -> Result<bool> {
        let mut raw_dimer = self.inner.take().ok_or(format_err!("raw dimer is empty"))?;

        // avoid trial rotation if estimated rotational angle `phi_est` is small
        // enough (Eq32 Heyden2005JCP)
        let phi_tol = self.vars.min_rot_angle;
        let dimer_state = raw_dimer.extrapolate();
        let phi_est = dimer_state.estimated_rotational_angle();
        let converged = check_dimer_rotation_convergence(phi_est, phi_tol);
        if converged {
            return Ok(true);
        }

        // trail rotation: use a fixed angle or variable one. See p12 in Heyden2005JCP
        let phi1 = if self.vars.use_fixed_rot_angle {
            self.vars.trial_rot_angle
        } else {
            self.vars.trial_rot_angle.min(phi_est)
        };
        let c0 = dimer_state.curvature();
        let c0d = dimer_state.curvature_derivative();

        // FIXME: refactor
        let r1_prime = raw_dimer.get_dimer_trial_rotation_endpoint(&self.orientation, &theta, phi1);
        self.dynamics.set_position(r1_prime.as_slice());
        let f1_prime = self.dynamics.get_force()?.to_vector();
        let mut trial_dimer = raw_dimer.clone();
        trial_dimer.r1 = r1_prime;
        trial_dimer.f1 = f1_prime;
        let state = trial_dimer.extrapolate();
        let c1 = state.curvature();

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

        // FIXME: refactor
        let f1_min = RawDimer::get_extrapolated_forces(phi1, phi_min, &raw_dimer.f1, &trial_dimer.f1);
        let r1_min = raw_dimer.get_dimer_trial_rotation_endpoint(&self.orientation, &theta, phi_min);
        raw_dimer.r1 = r1_min;
        raw_dimer.f1 = f1_min;
        raw_dimer.c0 = curvature_min.into();
        let state = raw_dimer.extrapolate();
        // FIXME: rewrite
        self.orientation = state.dimer_axis().clone();
        self.inner = raw_dimer.clone().into();

        let converged = check_dimer_rotation_convergence(phi_min, phi_tol);
        Ok(converged)
    }
}
// 69cb7fbe ends here

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
