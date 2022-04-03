// [[file:../dimer.note::782181ce][782181ce]]
use super::*;
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
    /// Carry out optimization in Dimer algorithm, and return the total energy and forces.
    pub fn optimize(&mut self) -> Result<DimerOutput> {
        let (c_min, t_min) = self.get_optimal_rotation(self.vars.max_num_rot)?;
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
