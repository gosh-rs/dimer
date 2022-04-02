// [[file:../dimer.note::782181ce][782181ce]]
use super::*;

use crate::cg::CG;
// 782181ce ends here

// [[file:../dimer.note::df98a463][df98a463]]
// 主要用于优化收敛控制
pub struct DimerOutput {
    /// DIMER energy, which is also equal to real energy
    pub total_energy: f64,
    /// DIMER forces on projection
    pub total_forces: Vec<f64>,
    /// DIMER force norm criterion
    pub fmax: f64,
    /// Real forces without DIMER projection
    pub fmax_real: f64,
    /// Lowest curvature
    pub curvature: f64,
    /// Lowest curvature mode
    pub curvature_mode: Vec<f64>,
}

impl<'a> Dimer<'a> {
    ///求得R0对应的最优curvature及curvature mode(或tangent)
    pub(crate) fn get_optimal_rotation(&mut self) -> Result<(f64, DVector)> {
        let tau_ini = self.orientation.clone();

        let mut niter = 0;
        let mut curvature_min = 0.0;

        let mut cg = crate::cg::CG::default();
        loop {
            niter += 1;
            info!("dimer rotation iteration {}", niter);
            let f_rot = self.re_orientate()?;
            let theta = self.get_rotational_direction(&f_rot, &mut cg);
            if self.next_rotation_step(&theta)? {
                info!("Optimal dimer rotation found within {} iterations.", niter);
                break;
            }
            if niter >= self.vars.max_num_rot {
                info!("Max allowed iterations reached");
                break;
            }
        }
        let phi = self.orientation.cosine_similarity(&tau_ini).acos();
        info!("Total rotational angle = {:.2}°", phi.to_degrees());

        // Total rotation angle during rotation steps
        // FIXME: rewrite
        let curvature_min = self.inner.as_ref().unwrap().c0.unwrap();
        let tau_min = self.orientation.clone();

        Ok((curvature_min, tau_min))
    }

    /// Return the total energy and forces for DIMER optimization.
    pub fn evaluate(&mut self) -> Result<DimerOutput> {
        let (c_min, t_min) = self.get_optimal_rotation()?;
        let (total_energy, fmax_real, fall) = self.next_translation_step(c_min, &t_min)?;
        let total_forces = fall.as_slice().to_vec();
        // let fmax = total_forces.iter().map(|x| x.vec2norm()).float_max();
        let fmax = todo!();

        Ok(DimerOutput {
            total_energy,
            total_forces,
            fmax,
            fmax_real,
            curvature: c_min,
            curvature_mode: t_min.as_slice().to_vec(),
        })
    }
}
// df98a463 ends here
