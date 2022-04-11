// [[file:../dimer.note::782181ce][782181ce]]
use super::*;
// 782181ce ends here

// [[file:../dimer.note::df98a463][df98a463]]
/// Optimized results in DIMER algorithm
pub struct DimerOutput {
    /// DIMER energy, which is equal to potential energy when at dimer center.
    pub total_energy: f64,
    /// The effective force felt at dimer center (modified force for DIMER translation)
    pub effective_force: Vec<f64>,
    /// The optimized lowest curvature
    pub curvature: f64,
    /// The optimized lowest curvature mode
    pub curvature_mode: Vec<f64>,
}

/// Main entry point for DIMER algorithm.
impl<'a> Dimer<'a> {
    /// Carry out optimization in Dimer algorithm, and return the total energy and forces.
    pub fn evaluate(&mut self) -> Result<DimerOutput> {
        let rotation = self.next_rotation_step(self.vars.max_num_rot)?;
        let mut raw_dimer = rotation.raw_dimer;
        let c_min = rotation.curvature_min;
        let effective_force = self.next_translation_step(&mut raw_dimer, c_min).as_slice().to_vec();

        Ok(DimerOutput {
            effective_force,
            curvature: c_min,
            total_energy: rotation.energy,
            curvature_mode: self.orientation.as_slice().to_vec(),
        })
    }
}
// df98a463 ends here
