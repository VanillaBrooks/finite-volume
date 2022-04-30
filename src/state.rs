use crate::prelude::*;

pub(crate) struct State {
    pub(crate) ic: InitialConditions,
    pub(crate) n: usize,
    pub(crate) h: f64,
    pub(crate) dt: f64,
    pub(crate) alpha: f64,
    pub(crate) gamma: f64,
    pub(crate) rho: Vector,
    pub(crate) p: Vector,
    pub(crate) u: Vector,
    pub(crate) e: Vector,
    pub(crate) E: Vector,
    pub(crate) c: Vector,
}

impl State{
    pub(crate) fn new(n: usize, h: f64, alpha: f64, dt:f64) -> Self {
        
		let gamma = 1.4;

		let ic = InitialConditions::new(gamma);

		let mut rho = Vector::zeros(n);
		let mut p = Vector::zeros(n);
		let mut u = Vector::zeros(n);
		let mut e = Vector::zeros(n);
		let mut E = Vector::zeros(n); // e + 1/2 u^2
		let mut c = Vector::zeros(n);

		// the discontinuity should be exactly in the middle
		let half_n = (n-1) / 2;

        rho.slice_mut(s!(0..half_n)).fill(ic.rhoL);
        rho.slice_mut(s!(half_n..n)).fill(ic.rhoR);

        p.slice_mut(s!(0..half_n)).fill(ic.pL);
        p.slice_mut(s!(half_n..n)).fill(ic.pR);

        u.slice_mut(s!(0..half_n)).fill(ic.uL);
        u.slice_mut(s!(half_n..n)).fill(ic.uR);

        e.slice_mut(s!(0..half_n)).fill(ic.eL);
        e.slice_mut(s!(half_n..n)).fill(ic.eR);

        c.slice_mut(s!(0..half_n)).fill(ic.cL);
        c.slice_mut(s!(half_n..n)).fill(ic.cR);

        for i in 0..half_n {
            E[i] = ic.eL + 0.5 * ic.uL.powf(2.);
        }
        for i in half_n..n{
            E[i] = ic.eR + 0.5 * ic.uR.powf(2.);
        }

        Self { rho, p, u, e, E, c, alpha, h, n, ic, dt, gamma }
    }

    pub(crate) fn check_state(&self, location: &str) {
        if self.rho.has_negs() {
            panic!("rho has negative values at {location}")
        }
        if self.rho.has_negs() {
            panic!("rho has nans at {location}")
        }
        if self.rho.has_zero() {
            panic!("rho has zero at {location}")
        }

        if self.e.has_nans() {
            dbg!(&self.e);
            panic!("entropy has nans at {location}")
        }
        if self.e.has_negs() {
            dbg!(&self.e);
            panic!("entropy has negatives at {location}")
        }
    }
}

pub(crate) struct InitialConditions{
    pub(crate) uR: f64,
    pub(crate) uL: f64,
    //
    pub(crate) pR: f64,
    pub(crate) pL: f64,
    //
    pub(crate) rhoR: f64,
    pub(crate) rhoL: f64,
    //
    pub(crate) eR: f64,
    pub(crate) eL: f64,
    //
    pub(crate) cR: f64,
    pub(crate) cL: f64,
    //
}

impl InitialConditions {
    fn new(gamma: f64) -> Self {
		let uL = 0.;
		let uR = 0.;
		let pL = 10f64.powf(5.);
		let pR = 10f64.powf(3.);
		let rhoL  = 1.0;
		let rhoR = 0.01;
		// from the paper p = (gamma - 1) rho e
		// so you can calculate these
		let eL = pL / ((gamma - 1.) * rhoL);
		let eR = pR / ((gamma - 1.) * rhoR);

		// given in the paper again
		let cL = (gamma * pL/ rhoL).sqrt();
		let cR = (gamma * pR / rhoR).sqrt();

        Self {uR, uL, pR, pL, rhoR, rhoL, eR, eL, cR, cL }
    }
}
