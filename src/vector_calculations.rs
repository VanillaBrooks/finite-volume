use crate::prelude::*;


pub(crate) fn recalculate_q(q: &mut Matrix, s: &State) {
    for i in 0..s.n {
        q[[0, i]] = s.rho[i]
    }

    for i in 0..s.n {
        q[[1, i]] = s.rho[i] * s.u[i]
    }

    for i in 0..s.n {
        q[[2, i]] = s.rho[i] * s.E[i]
    }
}

pub(crate) fn recalculate_f(f: &mut Matrix, s: &State) {
    for i in 0..s.n {
        f[[0, i]] = s.rho[i] * s.u[i]
    }

    for i in 0..s.n {
        f[[1, i]] = (s.rho[i] * s.u[i].powf(2.)) + s.p[i]
    }

    for i in 0..s.n {
        f[[2, i]] = s.u[i] * (s.p[i] + (s.rho[i] * s.E[i]))
    }
}

pub(crate) fn recalculate_variables(q: &Matrix, s: &mut State) {
    for i in 0..s.n {
        s.rho[i] = q[[0, i]]
    }

    for i in 0..s.n {
        s.u[i] = q[[1, i]] / s.rho[i]
    }

    for i in 0..s.n {
        s.E[i] = q[[2, i]] / s.rho[i]
    }

    for i in 0..s.n {
        s.e[i] = s.E[i] - (0.5 * s.u[i].powf(2.));

        if s.E[i] < 0.5 * s.u[i].powf(2.) {
            println!("e will have negative at i = {i}");
        }
    }

    for i in 0..s.n {
        s.p[i] = (s.gamma - 1.) * s.rho[i] * s.e[i]
    }

    for i in 0..s.n {
        s.c[i] = (s.gamma * s.p[i] / s.rho[i]).sqrt();
    }
}
