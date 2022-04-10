// [[file:../dimer.note::b46a5a4d][b46a5a4d]]
use super::*;
// b46a5a4d ends here

// [[file:../dimer.note::5205fe0e][5205fe0e]]
/// The part for DIMER translation
impl<'a> Dimer<'a> {
    /// Evalute DIMER total energy and forces for translation step. Return total
    /// energy, fmax of the real forces, and total force vector.
    ///
    /// # Parameters
    ///
    /// * c_min: optimized curvature value in rotation step
    ///
    pub fn next_translation_step(&mut self, raw_dimer: &mut RawDimer, c_min: f64) -> Result<(f64, f64, DVector)> {
        // re-use the energy and forces evaluated at rotation step
        // let e0 = raw_dimer.e0;
        let e0 = todo!();
        let f0 = &raw_dimer.f0;
        let t_min = &self.orientation;

        // update gradient for dimer translation
        let fall = if c_min.is_sign_positive() {
            info!("drag up directly");
            -f0.vector_projection(t_min)
        } else {
            f0 - 2.0 * f0.vector_projection(t_min)
        };

        // calculate fmax of real forces for statistics
        // FIXME: rewrite
        // let fmax_real = f0.as_slice().iter().map(|x| x.vec2norm()).float_max();
        let fmax_real = todo!();

        Ok((e0, fmax_real, fall))
    }
}
// 5205fe0e ends here
