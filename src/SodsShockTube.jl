module SodsShockTube
	include("./state.jl")
	include("./vector_calculations.jl")
	include("./solver_methods.jl")
	include("./io_utils.jl")

	using .IoUtils: Stepper, create_stepper, create_writer, write_flowfield, close_flowfield
	using .StateContainers: State, create_state, InitialConditions, check_state
	using .Methods: LaxFriedrichs, LaxFriedrichsMarker, LaxWendroff, LaxWendroffMarker, LaxWendroffViscousFluxes, LaxWendroffViscousFluxesMarker, create_method, MethodMarker
	using .Methods: solver_step!
	using .Calculations: recalculate_q!, recalculate_F!, recalculate_variables!
	
	export main
	export LaxFriedrichsMarker, LaxWendroffMarker, LaxWendroffViscousFluxesMarker

	function main(solver_method_marker::M) where M <: MethodMarker
		# length in meters
		L = 10
		# this should be odd
		n = 1501
		h = L / (n-1)
		# sptial coordinates
		x = 0:h:L

		# TODO: experiement with this
		CFL = 0.1

		#
		# Initial conditions
		#

		F = zeros(3, n)
		q = zeros(3, n)

		# our state
		alpha = 0.01
		s = create_state(n, h, alpha)

		# set the method we are using
		method = create_method(solver_method_marker, n, s)

		umax = max(s.ic.cL, s.ic.cR, s.ic.uL, s.ic.uR)
		dt = CFL * h / umax
		#dt /= 100

		t = 0.0
		step_number = 0
		t_end = 3.9 / 1000
		total_steps = Int(ceil(t_end / dt))

		# create something to write flowfields
		stepper = create_stepper(10, total_steps)
		writer = create_writer(stepper, s.ic.cL, n)

		println("dt is ", dt, " and total steps ", total_steps)

		recalculate_q!(q, s)
		while step_number < total_steps
			t = step_number * dt
			step_number += 1
			println("step ", step_number, " / ", total_steps)

			write_flowfield(writer, step_number, s, t)

			# construct q vector according to paper
			# row is the variable, the column is the 
			# spatial values
			recalculate_F!(F, s)

			# solve using the chosen method
			solver_step!(method, q, F, dt, h, s)

			println("recalculating final variables")
			recalculate_variables!(q, s)
			println("finished recalc final variables")
			
			#println("finished step", step_number)
		end

		close_flowfield(writer)
	end


end # module

using .SodsShockTube
#main(LaxFriedrichsMarker())
#main(LaxWendroffMarker())
main(LaxWendroffViscousFluxesMarker())
