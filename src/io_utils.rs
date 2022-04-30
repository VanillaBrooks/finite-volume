use crate::prelude::*;

pub(crate) struct Stepper {
    io_steps: usize,
    total_writes: usize
}

impl Stepper {
    pub(crate) fn new(io_steps: usize, total_steps: usize) -> Self {
		let total_writes = if io_steps == 1 {
			total_steps / io_steps
        } else {
			total_steps / io_steps + 1
        };

        Self { total_writes , io_steps }
    }

    fn should_step(&self, step: usize) -> bool {
        (step % self.io_steps == 0) || step == 1
    }
}

pub(crate) struct FlowfieldWriter {
    n: usize,
    stepper: Stepper,
    curr_write_number: usize,
    h5_file: hdf5::File,
    velocity: hdf5::dataset::Dataset,
    p: hdf5::dataset::Dataset,
    rho: hdf5::dataset::Dataset,
    e: hdf5::dataset::Dataset,
    c: hdf5::dataset::Dataset,
    E: hdf5::dataset::Dataset,
    t: hdf5::dataset::Dataset,
}

impl FlowfieldWriter {
    pub(crate) fn new(stepper: Stepper, n: usize) -> Self {
        let filename = "./results/flowfield.h5";
        std::fs::remove_file(filename).ok();

        let h5_file = hdf5::File::create(filename).unwrap();

        let creator = |name| {
            h5_file.new_dataset::<f64>()
                .shape((stepper.total_writes, n))
                .create(name).unwrap()
        };

        let velocity = creator("velocity");
        let p = creator("pressure");
        let rho = creator("rho");
        let e = creator("entropy");
        let c = creator("mach");
        let E = creator("energy");

        let t = h5_file.new_dataset::<f64>()
            .shape((stepper.total_writes, 1))
            .create("time").unwrap();

        FlowfieldWriter { n, stepper, curr_write_number: 0, h5_file, velocity, p, rho, e, c, E, t }
    }

    pub(crate) fn write_flowfield(&mut self, step: usize, s: &State, t: f64) {
        if !self.stepper.should_step(step) {
            return
        }

        let n = self.n;
        let write = self.curr_write_number;

        let slice = s!(write, 0..n);

        let mut mach = Vector::zeros(s.n);
        for i in 0..s.n {
            mach[i] = s.u[i] / s.c[i]
        }
        let mut time = Vector::zeros(1);
        time.fill(t);


        self.velocity.write_slice(&s.u, slice).expect("failed to write to hdf5 velocity");
        self.p.write_slice(&s.p, slice).expect("failed to write to hdf5 pressure");
        self.rho.write_slice(&s.rho, slice).expect("failed to write to hdf5 density");
        self.e.write_slice(&s.e, slice).expect("failed to write to hdf5 entropy");
        self.c.write_slice(&mach, slice).expect("failed to write to hdf5 mach");
        self.E.write_slice(&s.E, slice).expect("failed to write to hdf5 energy");
        self.t.write_slice(&time, s!(write, ..)).expect("failed to write to hdf5 time");


        self.curr_write_number += 1;
    }

    pub(crate) fn close(self) {
        self.h5_file.close().unwrap()
    }
}
