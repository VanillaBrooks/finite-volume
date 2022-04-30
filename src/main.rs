#![allow(non_snake_case)]

mod io_utils;
mod prelude;
mod vector_calculations;
mod traits;
mod state;
mod solver_methods;

use prelude::*;

use io_utils::{Stepper, FlowfieldWriter};



fn run_solver<MM, M>(method_marker: MM) 
    where MM: MethodMarker<Method = M>,
    M : SolverMethod
{
    
    // length in meters
    let L = 10.;
    // this should be odd
    let n = 1501;
    let h = L / (n-1) as f64;
    // sptial coordinates
    let x = Vector::range(0.0, L, h);

    // TODO: experiement with this
    let CFL = 0.1;

    //
    // Initial conditions
    //
    
    let mut f = Matrix::zeros((3, n));
    let mut q = Matrix::zeros((3, n));

    // our state
    let alpha = 0.01;
    let mut s = State::new(n, h, alpha, 0.01);

    // set the method we are using
    let mut method = MM::create_method(n, &s);

    let umax = s.ic.cL;
    let div : usize= 100usize.pow(0);
    let dt = (CFL * h / umax) / (div as f64);
    s.dt = dt;

    let mut step_number = 0;
    let t_end = 3.9 / 1000.;
    let total_steps = (t_end / dt).ceil() as usize;

    // create something to write flowfields
    let stepper = Stepper::new(20 * div, total_steps);
    let mut writer = FlowfieldWriter::new(stepper, n);

    println!("dt is {} and total steps are {}", dt, total_steps);

    // calculate q according to the state variables `s`
    // after this, q does not need to be calculated explicitly
    // with this function because updates to q are handled in
    // `solver_step!`
    recalculate_q(&mut q, &s);

    while step_number < total_steps {
        let t = step_number as f64 * dt;
        step_number += 1;
        
        //println!("step {} / {} ", step_number, total_steps);

        writer.write_flowfield(step_number, &s, t);

        // adjust F to the new values of s that were pulled from u 
        // in the last time step
        recalculate_f(&mut f, &s);

        //solve using the chosen method
        //solver_step(method, q, F, dt, h, s)
        method.solver_step(&mut q, &f, &mut s);

        // update the state `s` according to `q` that was calculated
        // in solver_step!()
        recalculate_variables(&q, &mut s);
        s.check_state("recalculate variables main loop");
    }

    writer.close()
}

fn main() {
    //let method = LaxFriedrichsMarker;
    let method = LaxWendroffMarker;
    run_solver(method)
}
