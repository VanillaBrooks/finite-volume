module SodsShockTube
	include("./state.jl")
	include("./vector_calculations.jl")
	include("./solver_methods.jl")
	include("./io_utils.jl")

	using .IoUtils: Stepper, create_stepper, create_writer, write_flowfield, close_flowfield, create_stepper_num_writes
	using .StateContainers: State, create_state, InitialConditions, check_state
	using .Methods: LaxFriedrichs, LaxFriedrichsMarker, LaxWendroff, LaxWendroffMarker, LaxWendroffViscousFluxes, LaxWendroffViscousFluxesMarker, create_method, MethodMarker
	using .Methods: solver_step!, unique_name
	using .Calculations: recalculate_q!, recalculate_F!, recalculate_variables!
	
	export main
	export LaxFriedrichsMarker, LaxWendroffMarker, LaxWendroffViscousFluxesMarker

	function main(solver_method_marker::M, alpha::T, n::Int) where M <: MethodMarker where T <: AbstractFloat
		# length in meters
		L = 10
		# this should be odd
		#n = 1501
		#n = 601
		h = L / (n-1)
		# sptial coordinates
		x = 0:h:L

		# TODO: experiement with this
		CFL = 0.3

		#
		# Initial conditions
		#

		F = zeros(3, n)
		q = zeros(3, n)

		# our state
		s = create_state(n, h, alpha)

		# set the method we are using
		method = create_method(solver_method_marker, n, s)

		umax = max(s.ic.cL, s.ic.cR, s.ic.uL, s.ic.uR)
		dt = CFL * h / umax

		t = 0.0
		step_number = 0
		t_end = 3.9 / 1000
		total_steps = Int(ceil(t_end / dt))

		# create something to write flowfields
		#stepper = create_stepper(20 * div, total_steps)
		stepper = create_stepper_num_writes(1, 2)
		flowfield_name = unique_name(solver_method_marker, alpha, n)
		writer = create_writer(stepper, s.ic.cL, n, flowfield_name)

		println("dt is ", dt, " and total steps ", total_steps, " alpha is ", alpha)

		# calculate q according to the state variables `s`
		# after this, q does not need to be calculated explicitly
		# with this function because updates to q are handled in
		# `solver_step!`
		recalculate_q!(q, s)

		write_flowfield(writer, step_number, s, t)
		while step_number < total_steps
			t = step_number * dt
			step_number += 1
			#println("step ", step_number, " / ", total_steps)


			# adjust F to the new values of s that were pulled from u 
			# in the last time step
			recalculate_F!(F, s)

			# # solve using the chosen method
			solver_step!(method, q, F, dt, h, s)

			# update the state `s` according to `q` that was calculated
			# in solver_step!()
			#
			# For Lax Wendroff II Method, this is where the crash occurs
			recalculate_variables!(q, s)
		end

		recalculate_variables!(q, s, update_speed_of_sound=true)
		write_flowfield(writer, step_number, s, t)
		close_flowfield(writer)
	end
end # module

using .SodsShockTube
n = 601
main(LaxFriedrichsMarker(), 0.0, n)
main(LaxWendroffMarker(), 0.0, n)

alpha_values = [0.0001, 0.01, 0.2, 0.5]

for alpha in alpha_values 
	main(LaxWendroffViscousFluxesMarker(), alpha, n)
end


# grid resolution study
for n in [101, 201, 401, 801]
	main(LaxFriedrichsMarker(), 0.0, n)
end
