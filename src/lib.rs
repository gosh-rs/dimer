// [[file:../dimer.note::c6f8257d][c6f8257d]]
mod cg;
mod dimer;
mod fourier;
mod options;
mod raw;
mod rotation;
mod translation;
// c6f8257d ends here

// [[file:../dimer.note::1e3853ed][1e3853ed]]
// #![deny(warnings)]

use gut::prelude::*;
use std::f64::consts::PI;
use vecfx::*;
use gosh::optim::Dynamics;

type DVector = nalgebra::DVector<f64>;

pub use options::UserOptions;
use raw::*;
// 1e3853ed ends here

// [[file:../dimer.note::a7df26ce][a7df26ce]]
pub struct Dimer<'a> {
    /// Potential for evaluation energy and forces
    dynamics: Dynamics<'a>,

    /// Dimer algorithm parameters
    vars: UserOptions,
    /// dimer orientation unit vector
    orientation: DVector,
    /// position vector of dimer center
    center: DVector,
    /// raw dimer struct
    inner: Option<RawDimer>,
}
// a7df26ce ends here
