// [[file:../dimer.note::3de9d57e][3de9d57e]]
use super::*;
// 3de9d57e ends here

// [[file:../dimer.note::6899fb9b][6899fb9b]]
#[derive(Debug, serde::Deserialize, serde::Serialize)]
struct RawData {
    r0: Vec<f64>,
    r1: Vec<f64>,
    f0: Vec<f64>,
    f1: Vec<f64>,
    dr: f64,
}

#[test]
fn test_dimer_rotation() -> Result<()> {
    let s = gut::fs::read_file("./tests/files/HCN.json")?;
    let raw_data: RawData = serde_json::from_str(&s)?;
    let r0 = raw_data.r0.to_vector();
    let r1 = raw_data.r1.to_vector();
    let f0 = raw_data.f0.to_vector();
    let f1 = raw_data.f1.to_vector();
    let tau = &r1 - &r0;
    let dr = raw_data.dr;
    let raw_dimer = RawDimer { r0, r1, f0, f1 };
    let state = raw_dimer.extrapolate();
    let c0 = state.curvature();
    approx::assert_relative_eq!(c0, -27.990186185223447, epsilon = 1e-6);
    let c0d = state.curvature_derivative();
    approx::assert_relative_eq!(c0d, -24.810830783954334, epsilon = 1e-6);
    let angle = state.estimated_rotational_angle();
    approx::assert_relative_eq!(angle, 0.20859479894110086, epsilon = 1e-6);

    let fr = state.rotational_force();
    #[rustfmt::skip]
    let fr_expected = [20.60444040,  -8.03260917,  0.00000000,
                       -8.25168033,  -5.86648454,  0.00000000,
                      -12.35271671,  13.89905034,  0.00000000];
    approx::assert_relative_eq!(fr, &fr_expected.to_vector(), epsilon = 1e-6);
    let cx = state.curvature_mode();
    #[rustfmt::skip]
    let cx_expected = [  0.504815,  -0.380250,   0.000000,
                        -0.394706,  -0.234350,   0.000000,
                        -0.110112,   0.614600,   0.000000];
    approx::assert_relative_eq!(cx, &cx_expected.to_vector(), epsilon = 1e-6);

    Ok(())
}
// 6899fb9b ends here
