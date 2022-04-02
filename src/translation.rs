// [[file:../dimer.note::b46a5a4d][b46a5a4d]]
use super::*;
// b46a5a4d ends here

// [[file:../dimer.note::5205fe0e][5205fe0e]]
impl<'a> Dimer<'a> {
    /// Evalute DIMER total energy and forces for translation step. Return total
    /// energy, fmax of the real forces, and total force vector.
    ///
    /// # Parameters
    ///
    /// * c_min: optimized curvature value in rotation step
    /// * t_min: optimized dimer mode in rotation step
    ///
    pub fn next_translation_step(&mut self, c_min: f64, t_min: &DVector) -> Result<(f64, f64, DVector)> {
        // reset internal dimer state
        let raw_dimer = self.inner.take().expect("raw dimer for trans");

        // re-use the energy and forces evaluated at rotation step
        let e0 = raw_dimer.e0;
        // FIXME: rewrite
        let f0 = &raw_dimer.f0;

        // update gradient for dimer translation
        let fall = if c_min.is_sign_positive() {
            info!("drag up directly");
            -f0.vector_projection(&t_min)
        } else {
            f0 - 2.0 * f0.vector_projection(&t_min)
        };

        // calculate fmax of real forces for statistics
        // FIXME: rewrite
        // let fmax_real = f0.as_slice().iter().map(|x| x.vec2norm()).float_max();
        let fmax_real = todo!();

        Ok((e0, fmax_real, fall))
    }
}
// 5205fe0e ends here
