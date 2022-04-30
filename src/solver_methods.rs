use crate::prelude::*;

pub(crate) struct LaxFriedrichsMarker;

impl MethodMarker for LaxFriedrichsMarker {
    type Method = LaxFriedrichs;

    fn create_method(_n: usize, _s: &State) -> Self::Method {
        LaxFriedrichs
    }
}
pub(crate) struct LaxFriedrichs;

impl SolverMethod for LaxFriedrichs {
    fn solver_step(&mut self, q: &mut Matrix, f: &Matrix, s: &mut State) {
        let q_copy = q.clone();
        for j in 1..s.n-1 {
            for v in 0..3 {
                let left = 0.5 * (q_copy[[v, j+1]] + q_copy[[v, j-1]] );
                let right = (s.dt / (2. * s.h)) * (f[[v, j+1]] - f[[v, j-1]]);

                q[[v,j]] = left - right;
            }
        }
    }
}

pub(crate) struct LaxWendroffMarker;

impl MethodMarker for LaxWendroffMarker {
    type Method = LaxWendroff;

    fn create_method(n: usize, s: &State) -> Self::Method {
        let mut q_star = Matrix::zeros((3, n));
        let mut f_star = Matrix::zeros((3, n));

        recalculate_q(&mut q_star, s);
        recalculate_f(&mut f_star, s);

        LaxWendroff {
            q_star, 
            f_star
        }
    }
}

pub(crate) struct LaxWendroff {
    q_star: Matrix,
    f_star: Matrix
}

impl SolverMethod for LaxWendroff {
    fn solver_step(&mut self, q: &mut Matrix, f: &Matrix, s: &mut State) {
        //perform first step
        for j in 1..s.n-1 {
            for v in 0..3 {
                let left = 0.5 * (q[[v, j]] + q[[v,j+1]]);
                let right = (s.dt / (2. * s.h)) * (f[[v,j+1]] - f[[v, j]]);

                self.q_star[[v,j]] = left - right;
            }
        }

        // pull u/ rho / e/ E ... etc out of q* so we can use them
        // to calculate F*
        recalculate_variables(&self.q_star, s);
        // use the new variables to calculate F*
        recalculate_f(&mut self.f_star, s);

        // perform second step
        for j in 1..s.n-1 {
            for v in 0..3 {
                let left = q[[v,j]];
                let right = (s.dt / s.h) * (self.f_star[[v, j]] - self.f_star[[v, j-1]]);

                q[[v,j]] = left - right;
            }
        }
    }
}
