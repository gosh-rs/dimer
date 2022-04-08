// [[file:../dimer.note::c6f8257d][c6f8257d]]
mod cg;
mod dimer;
mod fourier;
mod options;
mod raw;
mod rotation;
mod translation;

#[cfg(test)]
mod test;
// c6f8257d ends here

// [[file:../dimer.note::1e3853ed][1e3853ed]]
// #![deny(warnings)]

use gosh::optim::{Dynamics, EvaluateEnergyForce};
use gut::prelude::*;
use std::f64::consts::PI;
use vecfx::*;

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

impl<'a> Dimer<'a> {
    /// Construct a dimer from center and axis orientation.
    pub fn new(center: &[f64], orientation: &[f64], pot: impl EvaluateEnergyForce + 'a) -> Self {
        assert_eq!(center.len(), orientation.len(), "invalid data: {center:?}, {orientation:?}");
        let dynamics = Dynamics::new(center, pot);
        let orientation = orientation.to_vector().normalize();
        Self {
            center: center.to_vector(),
            dynamics,
            orientation,
            vars: UserOptions::default(),
            inner: None,
        }
    }
}
// a7df26ce ends here

// [[file:../dimer.note::cfd3ba0e][cfd3ba0e]]
#[cfg(feature = "adhoc")]
/// Docs for local mods
pub mod docs {
    macro_rules! export_doc {
        ($l:ident) => {
            pub mod $l {
                pub use crate::$l::*;
            }
        };
    }

    export_doc!(raw);
    export_doc!(fourier);
    export_doc!(options);
    export_doc!(rotation);
    export_doc!(translation);
    export_doc!(cg);
}
// cfd3ba0e ends here
