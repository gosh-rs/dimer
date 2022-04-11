// [[file:../dimer.note::6b0d4c47][6b0d4c47]]
use dimer::Dimer;

use approx::*;
use gchemol::prelude::*;
use gchemol::Molecule;
use gosh::gchemol;
use gosh::model::BlackBoxModel;
use gosh::prelude::ChemicalModel;
use gut::prelude::*;
use vecfx::*;

#[test]
fn test_dimer_opt() -> Result<()> {
    gut::cli::setup_logger_for_test();

    let path = "./tests/files/HCN.xyz";
    let mut mol = Molecule::from_file(path)?;
    let center = mol.positions().collect_vec();
    let orientation = mol.velocities().collect_vec();
    let mut bbm = BlackBoxModel::from_dir("/share/apps/xtb/sp")?;
    let pot = |position: &[f64], force: &mut [f64]| {
        mol.set_positions(position.as_3d().to_vec());
        let mp = bbm.compute(&mol)?;
        let energy = mp.get_energy().unwrap();
        let forces = mp.get_forces().unwrap();
        force.clone_from_slice(forces.as_flat());
        Ok(energy)
    };

    let mut dimer = Dimer::new(center.as_flat(), orientation.as_flat(), pot);
    dimer.vars.min_rot_angle = 1f64.to_radians();
    dimer.vars.max_num_rot = 10;
    dimer.vars.max_num_trans = 50;
    dimer.vars.use_extrapolated_force = false;
    dimer.vars.use_fixed_rot_angle = true;
    dimer.vars.use_cg_rot = false;

    let o = dimer.evaluate()?;
    #[rustfmt::skip]
    let effective_force_expected = [  2.66399668,   2.98157015,   0.00000000,
                                      1.22798397,  -2.85601768,   0.00000000,
                                     -3.89198253,  -0.12555248,   0.00000000].to_vector();
    assert_relative_eq!(o.effective_force.to_vector(), effective_force_expected, epsilon = 1e-6);
    let c_min_expected = -32.89828875128214;
    assert_relative_eq!(o.curvature, c_min_expected, epsilon = 1e-6);
    #[rustfmt::skip]
    let c_mode_expected = [  0.52279381,  -0.38265624,   0.00000000,
                            -0.31599899,  -0.23553730,   0.00000000,
                            -0.20679556,   0.61819354,   0.00000000].to_vector();
    assert_relative_eq!(o.curvature_mode.to_vector(), c_mode_expected, epsilon = 1e-6);

    Ok(())
}
// 6b0d4c47 ends here
