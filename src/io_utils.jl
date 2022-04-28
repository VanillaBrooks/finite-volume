module IoUtils
	using HDF5
	export Stepper, create_stepper, create_writer, write_flowfield, close_flowfield

	struct Stepper
		io_steps::Int
		total_writes::Int
	end

	function create_stepper(io_steps::Int, total_steps::Int)::Stepper
		if io_steps == 1
			total_writes = Int(floor(total_steps / io_steps))
		else
			total_writes = Int(floor(total_steps / io_steps)) + 1
		end

		return Stepper(io_steps, total_writes)
	end

	# check if we should write output now
	function should_step(stepper::Stepper, step::Int)
		return (step % stepper.io_steps == 0) || step == 1
	end

	mutable struct FlowfieldWriter
		n::Int
		h5file::HDF5.File
		u_dset::HDF5.Dataset
		p_dset::HDF5.Dataset
		t_dset::HDF5.Dataset
		ρ_dset::HDF5.Dataset
		e_dset::HDF5.Dataset
		c_dset::HDF5.Dataset
		E_dset::HDF5.Dataset
		stepper::Stepper
		write_number::Int
	end

	function create_writer(stepper::Stepper, out::T, n::Int) where T <: AbstractFloat
		filename = "./results/flowfield.h5"
		Base.Filesystem.rm(filename, force=true)

		flowfield_file_id = h5open(filename, "w")

		u_dset = create_dataset(flowfield_file_id, "velocity", datatype(T), dataspace(stepper.total_writes, n))
		p_dset= create_dataset(flowfield_file_id, "pressure", datatype(T), dataspace(stepper.total_writes, n))
		ρ_dset = create_dataset(flowfield_file_id, "rho", datatype(T), dataspace(stepper.total_writes, n))
		t_dset = create_dataset(flowfield_file_id, "time", datatype(T), dataspace(stepper.total_writes, 1))
		e_dset = create_dataset(flowfield_file_id, "entropy", datatype(T), dataspace(stepper.total_writes, n))
		c_dset = create_dataset(flowfield_file_id, "mach", datatype(T), dataspace(stepper.total_writes, n))
		E_dset = create_dataset(flowfield_file_id, "energy", datatype(T), dataspace(stepper.total_writes, n))
		
		return FlowfieldWriter(
			n,
			flowfield_file_id,
			u_dset,
			p_dset,
			t_dset,
			ρ_dset,
			e_dset,
			c_dset,
			E_dset,
			stepper,
			1
		)
	end

	function write_flowfield(
		writer::FlowfieldWriter, 
		step:: Int, 
		u::Vector{T}, 
		p::Vector{T}, 
		t::T,
		ρ::Vector{T},
		e::Vector{T},
		c::Vector{T},
		E::Vector{T}
	) where T <: AbstractFloat
		# break early if we are not supposed to step now
		if !should_step(writer.stepper, step)
			return
		end

		println("writing flowfield for step ", step, " with io steps as ", writer.stepper.io_steps, " and total writes ", writer.stepper.total_writes, "write number is ", writer.write_number)

		writer.u_dset[writer.write_number, :] = u
		writer.p_dset[writer.write_number, :] = p
		writer.t_dset[writer.write_number, 1] = t
		writer.ρ_dset[writer.write_number, :] = ρ
		writer.e_dset[writer.write_number, :] = e
		writer.c_dset[writer.write_number, :] = c
		writer.E_dset[writer.write_number, :] = E


		writer.write_number += 1
	end

	function close_flowfield(
		writer::FlowfieldWriter
	)
		close(writer.h5file)
	end
end
