// [[file:../dimer.note::c6f8257d][c6f8257d]]
mod cg;
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

type DVector = nalgebra::DVector<f64>;

pub use options::UserOptions;
use raw::*;
// 1e3853ed ends here

// [[file:../dimer.note::a7df26ce][a7df26ce]]
pub struct Dimer {
    vars: UserOptions,
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
