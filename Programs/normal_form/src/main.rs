use num;
use num::Rational32;
use anyhow::{bail};

type Rational = Rational32;
type Complex = num::Complex<i32>;

// type Point = num::Complex<f64>;

type Point = (Complex, Complex, Complex);
type ProjTrans = [[Complex; 3]; 3];
type Cubic = [[i32; 3]; 3]; // a_ij x^i y^j z^(3-i-j)

// type Cubic<T> = [[T; 3]; 3]; // a_ij x^i y^j z^(3-i-j)
// type Cubic = Cubic<i32>;

// 6 a_30 x + 2 a_21 y + 2 a_20 z

// 3 a_30 x^2 + a_10 z^2 + a_11 y z + 2 a_21 x y + a_12 y^2 + 2 a_20 x z

// a_30 x^3 + a_03 y^3 + a_00 z^3 + a_01 y z^2 + a_10 x z^2 + a_11 x y z + a_21 x^2 y + a_12 x y^2 + a_20 x^2 z + a_02 y^2 z

// fn partial_diff_x(cubic: Cubic) -> [[i32; 2]; 2] {
//     // i a_ij x^(i-1) y^j z^(3-i-j)
//     // 3 a_30 x^2 + a_10 z^2 + a_11 y z + 2 a_21 x y + a_12 y^2 + 2 a_20 x z
// }

// fn partial_diff_y(cubic: Cubic) -> [[i32; 2]; 2] {
//     // j a_ij x^i y^(j-1) z^(3-i-j)
//     // 3 a_03 y^2 + a_01 y z^2 + a_10 x z^2 + a_11 x y z + a_21 x^2 y + a_12 x y^2 + a_20 x^2 z + a_02 y^2 z
// }

// fn is_cubic_singular(cubic: Cubic) -> bool {}

// Is rational point on cubic always inflection point?

fn gessian(cubic: Cubic) -> Cubic {
    let a = cubic;

    let mut gessian = [[0; 3]; 3];

    gessian[3][0] = 8 * a[1][1] * a[2][0] * a[2][1]
        - 8 * a[1][0] * a[2][1].pow(2)
        - 6 * a[1][1].pow(2) * a[3][0]
        - 8 * a[1][2] * (a[2][0].pow(2) - 3 * a[1][0] * a[3][0]);

    gessian[0][0] = 8 * a[0][1] * a[1][0] * a[1][1]
        - 6 * a[0][0] * a[1][1].pow(2)
        - 8 * a[0][1].pow(2) * a[2][0]
        - 8 * a[0][2] * (a[1][0].pow(2) - 3 * a[0][0] * a[2][0]);

    gessian[0][3] = -6 * a[0][3] * (a[1][1].pow(2) - 4 * a[0][1] * a[2][1])
        - 8 * ((-a[0][2]) * a[1][1] * a[1][2]
            + a[0][1] * a[1][2].pow(2)
            + a[0][2].pow(2) * a[2][1]);

    gessian[0][1] = -24 * a[0][3] * (a[1][0].pow(2) - 3 * a[0][0] * a[2][0])
        + 2 * a[0][1]
            * (a[1][1].pow(2) + 8 * a[1][0] * a[1][2]
                - 4 * a[0][2] * a[2][0])
        - 8 * a[0][1].pow(2) * a[2][1]
        + 24 * a[0][0] * ((-a[1][1]) * a[1][2] + a[0][2] * a[2][1]);

    //     // a_30 x^3 + a_03 y^3 + a_00 z^3 + a_01 y z^2 + a_10 x z^2 + a_11 x y z + a_21 x^2 y + a_12 x y^2 + a_20 x^2 z + a_02 y^2 z

    //     // x = x1, y = x2, z = x3
    //     // gessian[3][0] = 2
    //     //     * (-4 * a[1][2] * a[2][0].pow(2) + 4 * a[1][1] * a[2][0] * a[2][1]
    //     //         - 4 * a[1][0] * a[2][1].pow(2)
    //     //         - 3 * a[1][1].pow(2) * a[3][0]
    //     //         + 12 * a[1][0] * a[1][2] * a[3][0]);

    //     // gessian[0][0] = 2
    //     //     * (-4 * a[0][2] * a[1][0].pow(2) + 4 * a[0][1] * a[1][0] * a[1][1]
    //     //         - 3 * a[0][0] * a[1][1].pow(2)
    //     //         - 4 * a[0][1].pow(2) * a[2][0]
    //     //         + 12 * a[0][0] * a[0][2] * a[2][0]);

    //     // gessian[0][3] = 2
    //     //     * (-3 * a[0][3] * a[1][1].pow(2) + 4 * a[0][2] * a[1][1] * a[1][2]
    //     //         - 4 * a[0][1] * a[1][2].pow(2)
    //     //         - 4 * a[0][2].pow(2) * a[2][1]
    //     //         + 12 * a[0][1] * a[0][3] * a[2][1]);

    //     // [[], []];
    //     todo!();
}

// fn all_intersection_points(cubic1: Cubic, cubic2: Cubic) -> [Point; 9] {
//     todo!();
// }

// fn find_all_inflection_points(cubic: Cubic) -> [Point; 9] {
//     todo!();
// }

// fn find_rational_point(cubic: Cubic) -> Rational {
//     todo!();
// }

// Not needed for our problem, but might be useful for more general problems
fn is_point_of_inflection(cubic: Cubic, point: Point) -> bool {
    return false;
}

// fn step1(
// import subprocess
// code = 'ToString[' + match.group(1) + ', TeXForm]'
// snip.rv = subprocess.check_output(['wolframscript', '-code', code])

fn normal_form(cubic: Cubic, point: Point) -> (Complex, Complex, ProjTrans) {
    let A = [
        [
            Rational::from_integer(1),
            Rational::from_integer(1),
            Rational::from_integer(0),
        ],
        [
            Rational::new(1, 2),
            Rational::new(-1, 2),
            Rational::from_integer(0),
        ],
        [
            Rational::from_integer(0),
            Rational::from_integer(0),
            Rational::from_integer(1),
        ],
    ];

    //// Not needeed for our problem
    // if !is_cubic_singular(cubic) {
    //     return bail!("Cubic is not regular");
    // }

    if is_point_of_inflection(cubic, point) {
        unimplemented!();
    }

    todo!();
}

fn main() {
    // a^3 + b^3 + c^3 + (1 - N) (a^2 b + a^2 c + b^2 a + b^2 c + c^2 a + c^2 b) + (3 - 2N) a b c = 0
    // rational point = (1 : -1 : 0), it is inflection point

    // a_30 x^3 + a_03 y^3 + a_00 z^3 + a_01 y z^2 + a_10 x z^2 + a_11 x y z + a_21 x^2 y + a_12 x y^2 + a_20 x^2 z + a_02 y^2 z

    const N: i32 = 4;

    #[rustfmt::skip]
    let cubic = [
        [1, 1 - N, 1],
        [1 - N, 3 - 2 * N, 1 - N],
        [1, 1 - N, 1 - N]
    ];

    let point = (Complex::new(1, 0), Complex::new(-1, 0), Complex::new(0, 0));

    normal_form(cubic, point);
}
