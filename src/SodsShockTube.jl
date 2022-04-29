module SodsShockTube
	include("./state.jl")
	include("./solver_methods.jl")
	include("./io_utils.jl")

	using .IoUtils: Stepper, create_stepper, create_writer, write_flowfield, close_flowfield
	using .StateContainers: State, create_state, InitialConditions, check_state
	using .Methods: LaxFriedrichs, LaxFriedrichsMarker, LaxWendroff, LaxWendroffMarker, LaxWendroffViscousFluxes, LaxWendroffViscousFluxesMarker, create_method, MethodMarker
	
	export main
	export LaxFriedrichsMarker, LaxWendroffMarker, LaxWendroffViscousFluxes

	function main(solver_method_marker::M) where M <: MethodMarker
		# length in meters
		L = 10
		# this should be odd
		n = 121
		h = L / (n-1)
		# sptial coordinates
		x = 0:h:L

		# TODO: experiement with this
		CFL = 1.0

		#
		# Initial conditions
		#

		F = zeros(3, n)
		q = zeros(3, n)

		# our state
		s = create_state(n, h)

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
		stepper = create_stepper(30, total_steps)
		writer = create_writer(stepper, s.ic.cL, n)

		println("dt is ", dt, " and total steps ", total_steps)

		while step_number < total_steps
			t = step_number * dt
			step_number += 1
			println("step ", step_number, " / ", total_steps)

			write_flowfield(writer, step_number, s, t)

			# construct q vector according to paper
			# row is the variable, the column is the 
			# spatial values
			recalculate_q!(q, s)
			recalculate_F!(F, s)

			# solve using the chosen method
			solver_step!(method, q, F, dt, h, s)

			#println("recalculating final variables")
			recalculate_variables!(q, s)
			#println("finished recalc final variables")
			
			#println("finished step", step_number)
		end

		close_flowfield(writer)
	end

	function solver_step!(method::LaxFriedrichs, q::Matrix{T}, F::Matrix{T}, dt::T, h::T, _...) where T <: AbstractFloat
		# TODO: Lax Freidrichs method
		# according to the paper
		n = size(q)[2]

		for j = 2:n-1
			for v = 1:3
				left = (1/2) * (q[v, j+1] + q[v, j-1])
				right = (dt / (2*h)) * (F[v, j+1] - F[v, j-1])
				q[v, j] = left - right
			end
		end
	end


	function recalculate_variables!(
		q::Matrix{T},
		s::State{T},
	) where T <: AbstractFloat
		s.ρ .= q[1, :]
		s.u .= q[2, :] ./ s.ρ[:]
		s.E .= q[3, :] ./ s.ρ[:]
		s.e .= s.E .- (1/2) .* s.u.^2
		s.p .= (s.gamma - 1) .* s.ρ .* s.e

		check_state(s)

		s.c .= sqrt.(s.gamma .* s.p ./ s.ρ)
	end

	function recalculate_q!(
		q::Matrix{T},
		s::State,
	) where T <: AbstractFloat
		q[1, :] = s.ρ
		q[2, :] = s.ρ .* s.u
		q[3, :] = s.ρ .* s.E
	end

	function recalculate_F!(
		F::Matrix{T}, 
		s::State,
	) where T <: AbstractFloat
		F[1, :] = s.ρ .* s.u
		F[2, :] = (s.ρ .* s.u.^2) .+ s.p
		F[3, :] = s.u .* ((s.ρ .* s.E) .+ s.p)
	end

end # module

using .SodsShockTube
main(LaxFriedrichsMarker())
#main(LaxWendroffMarker())
