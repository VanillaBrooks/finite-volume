module SodsShockTube
	include("./io_utils.jl")

	export main
	export LaxFriedrichsMarker, LaxWendroffMarker, LaxWendroffViscousFluxes
	using .IoUtils: Stepper, create_stepper, create_writer, write_flowfield, close_flowfield

	abstract type MethodMarker end

	struct LaxFriedrichsMarker <: MethodMarker end

	struct LaxWendroffMarker <: MethodMarker end

	struct LaxWendroffViscousFluxesMarker <: MethodMarker end

	function main(solver_method_marker::M) where M <: MethodMarker
		#
		# parameters
		#

		# specific heat ratio for air
		gamma = 1.4
		# 1 -> LF
		# 2 -> LW-II
		# 3 -> something + artificial viscousity
		flag = 1 

		#
		# spatial dimension
		#

		# length in meters
		L = 10
		# this should be odd
		n = 99
		# the number of delta steps is one minus
		# the total number of points you have
		h = L / (n-1)
		# sptial coordinates
		x = 0:h:L


		#
		# Time step
		#

		# TODO: experiement with this
		CFL = 1.0
		# 3.9 ms

		#
		# Initial conditions
		#

		# these are given in the problem statelemt
		# from the paper
		uL = 0.;
		uR = 0.;
		pL = 10e5;
		pR = 10^3;
		ρL  = 1;
		ρR = 0.010;
		# form the paper p = (gamma - 1) rho e
		# so you can calculate these
		eL = pL / ((gamma - 1) * ρL);
		eR = pR / ((gamma - 1) * ρR);

		# given in the paper again
		cL = sqrt(gamma * pL/ ρL)
		cR = sqrt(gamma * pR / ρR)	

		#
		# states / variables
		#
		ρ = zeros(n)
		p = zeros(n)
		u = zeros(n)
		e = zeros(n)
		E = zeros(n) # e + 1/2 u^2
		c = zeros(n)

		F = zeros(3, n)
		q = zeros(3, n)

		# the discontinuity should be exactly in the middle
		half_n = Int(floor((n-1) / 2));

		ρ[1:half_n] .= ρL;
		ρ[half_n:n] .= ρR;

		p[1:half_n] .= pL;
		p[half_n:n] .= pR;

		u[1:half_n] .= uL;
		u[half_n:n] .= uR;

		e[1:half_n] .= eL;
		e[half_n:n] .= eR;

		E[1:half_n] .= eL + (1/2) .* uL.^2;
		E[half_n:n] .= eR + (1/2) .* uR.^2;

		c[1:half_n] .= cL
		c[half_n:n] .= cR

		# set the method we are using
		method = create_method(solver_method_marker, n)


		umax = max(cL, cR, uL, uR)
		dt = CFL * h / umax

		t = 0.0
		step_number = 0
		t_end = 3.9 / 1000
		total_steps = Int(ceil(t_end / dt))

		# create something to write flowfields
		stepper = create_stepper(1, total_steps)
		writer = create_writer(stepper, uL, n)

		println("dt is ", dt, " and total steps ", total_steps)

		while step_number < total_steps
			t = step_number * dt
			step_number += 1
			write_flowfield(writer, step_number, u, p, t, ρ)

			# construct q vector according to paper
			# row is the variable, the column is the 
			# spatial values
			recalculate_q!(q, u, ρ, E)
			recalculate_F!(F, u, ρ, e, p)

			# solve using the chosen method
			solver_step!(method, q, F, dt, h, u, ρ, E, e, p, gamma)

			recalculate_variables!(q, u, ρ, E, e, p, gamma)
			
			t += dt
		end

		close_flowfield(writer)


		# then, do the plotting below
	end

	struct LaxFriedrichs end

	struct LaxWendroff{T}
		q_star::Matrix{T}
		F_star::Matrix{T}
	end

	struct LaxWendroffViscousFluxes{T}
	end

	function create_method(marker::LaxFriedrichsMarker, n::Int)::LaxFriedrichs 
		return LaxFriedrichs()
	end

	function create_method(marker::LaxWendroffMarker, n::Int)::LaxWendroff
		q_star = zeros(3, n)
		F_star = zeros(3, n)
		return LaxWendroff(q_star, F_star)
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

	function solver_step!(
		method::LaxWendroff, 
		q::Matrix{T}, 
		F::Matrix{T}, 
		dt::T, 
		h::T,
		u::Vector{T},
		ρ::Vector{T},
		E::Vector{T},
		e::Vector{T},
		p::Vector{T},
		gamma::T
	) where T <: AbstractFloat
		n = size(q)[2]

		# first, calculate q_star
		for j = 2:n-1
			for v = 1:3
				left = (1/2) * (q[v, j] + q[v, j+1])
				right = (dt / (2*h)) * (F[v, j+1] - F[v, j])
				method.q_star[v, j] = left - right
			end
		end

		# then calculate the F* that comes from that q*
		recalculate_variables!(q, u, ρ, E, e, p, gamma)
		recalculate_F!(method.F_star, u, ρ, e, p)

		for j = 2:n-1
			for v = 1:3
				left = q[v, j]
				right = (dt / h)  * (method.F_star[j+1] - method.F_star[j-1])
				q[v, j] = left - right
			end
		end

		# # since we need j+1 we can only go to n-1
		# qstar[1:3, 1:n-1] = todo
		# # here the input to the flux is qstar
		# # (instead of q)
		# Fstar(1:3, 1:n-1] = 
		# # then do the second step of the equations
		# q[1:3, 1:n-1] = 
	end

	function solver_step!(method::LaxWendroffViscousFluxes)
	end

	function recalculate_variables!(
		q::Matrix{T},
		u::Vector{T},
		ρ::Vector{T},
		E::Vector{T},
		e::Vector{T},
		p::Vector{T},
		gamma::T
	) where T <: AbstractFloat
		ρ .= q[1, :]
		u .= q[2, :] ./ ρ[:]
		E .= q[3, :] ./ ρ[:]
		e .= E - (1/2) .* u.^2
		p .= (gamma - 1) .* ρ .* e
	end

	function recalculate_q!(
		q::Matrix{T},
		u::Vector{T},
		ρ::Vector{T},
		E::Vector{T},
		#e::Vector{T},
		#p::Vector{T}
	) where T <: AbstractFloat
		q[1, :] = ρ
		q[2, :] = ρ .* u
		q[3, :] = ρ .* E
	end

	function recalculate_F!(
		F::Matrix{T}, 
		u::Vector{T},
		ρ::Vector{T},
		e::Vector{T},
		p::Vector{T}
		#E::Vector{T},
	) where T <: AbstractFloat
		F[1, :] = ρ .* u
		F[2, :] = ρ .* u.^2
		# TODO: can probably just use E to calculate this here
		F[3, :] = u .* (ρ .* (e .+ (1/2) .* u.^2) .+ p)
	end
end # module

using .SodsShockTube
main(LaxFriedrichsMarker())
#main(LaxWendroffMarker())
