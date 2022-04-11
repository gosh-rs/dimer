// [[file:../dimer.note::3de9d57e][3de9d57e]]
use super::*;
// 3de9d57e ends here

// [[file:../dimer.note::917f277b][917f277b]]
#[derive(Debug, serde::Deserialize, serde::Serialize)]
struct RawData {
    r0: Vec<f64>,
    r1: Vec<f64>,
    f0: Vec<f64>,
    f1: Vec<f64>,
    dr: f64,
}

fn get_raw_dimer() -> Result<RawDimer> {
    let s = gut::fs::read_file("./tests/files/HCN.json")?;
    let raw_data: RawData = serde_json::from_str(&s)?;
    let r0 = raw_data.r0.to_vector();
    let r1 = raw_data.r1.to_vector();
    let f0 = raw_data.f0.to_vector();
    let f1 = raw_data.f1.to_vector();
    let tau = &r1 - &r0;
    let dr = raw_data.dr;
    let raw_dimer = RawDimer { r0, r1, f0, f1 };
    Ok(raw_dimer)
}
// 917f277b ends here

// [[file:../dimer.note::6899fb9b][6899fb9b]]
#[test]
fn test_raw_dimer() -> Result<()> {
    let raw_dimer = get_raw_dimer()?;
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

// [[file:../dimer.note::170e45af][170e45af]]
#[test]
fn test_fourier_rotation() -> Result<()> {
    let mut raw_dimer = get_raw_dimer()?;

    let phi1 = 0.7853981633974483;

    #[rustfmt::skip]
    let r1_prime = [ 0.69860200757930350,  0.31180993004655666, 0.0,
                    -0.40309271543630126,  1.37241379142601970, 0.0,
                    -0.29640328699113040, -0.28631672394427166, 0.0].to_vector();

    #[rustfmt::skip]
    let f1_prime = [1.688936429470,  3.4258501799580,  0.0,
                    1.608288755908, -2.5162030267550,  0.0,
                   -3.297225185378, -0.9096471532024,  0.0].to_vector();

    #[rustfmt::skip]
    let theta = [ 0.52191443789529680,  0.210445700894246520, 0.0,
                  0.22540289240012049,  0.055864255148702295, 0.0,
                 -0.74730702733915100, -0.266313451548034930, 0.0].to_vector();

    let state = raw_dimer.fourier_rotate(r1_prime, f1_prime, phi1, &theta, false);

    #[rustfmt::skip]
    let r1_min = [  0.698408,   0.311562,   0.000000,
                   -0.403355,   1.372309,   0.000000,
                   -0.295947,  -0.285964,   0.000000].to_vector();
    #[rustfmt::skip]
    let f1_min = [  1.775292,   3.477056,   0.000000,
                    1.608411,  -2.513140,   0.000000,
                   -3.383703,  -0.963917,   0.000000].to_vector();

    let phi_min = 0.053786278813014274;
    let curvature_min = -28.65807149654434;
    approx::assert_relative_eq!(state.curvature_min, curvature_min, epsilon = 1e-6);
    approx::assert_relative_eq!(state.phi_min, phi_min, epsilon = 1e-6);
    approx::assert_relative_eq!(state.r1_min, r1_min, epsilon = 1e-6);
    approx::assert_relative_eq!(state.f1_min, f1_min, epsilon = 1e-6);

    Ok(())
}
// 170e45af ends here
