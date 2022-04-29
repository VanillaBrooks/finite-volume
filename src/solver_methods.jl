module Methods
	using ..StateContainers: State

	export MethodMarker
	export LaxFriedrichs, LaxFriedrichsMarker
	export LaxWendroff, LaxWendroffMarker
	export LaxWendroffViscousFluxes, LaxWendroffViscousFluxesMarker
	export create_method

	abstract type MethodMarker end
	struct LaxFriedrichsMarker <: MethodMarker end
	struct LaxWendroffMarker <: MethodMarker end
	struct LaxWendroffViscousFluxesMarker <: MethodMarker end

	struct LaxFriedrichs end

	struct LaxWendroff{T}
		q_star::Matrix{T}
		F_star::Matrix{T}
	end

	struct LaxWendroffViscousFluxes{T}
	end

	function create_method(marker::LaxFriedrichsMarker, n::Int, _...)::LaxFriedrichs 
		return LaxFriedrichs()
	end

	function create_method(marker::LaxWendroffMarker, n::Int, 
		s::State{T}
	)::LaxWendroff where T <: AbstractFloat
		q_star = zeros(3, n)
		F_star = zeros(3, n)

		recalculate_q!(q_star, s)
		recalculate_F!(F_star, s)

		return LaxWendroff(q_star, F_star)
	end

	function solver_step!(
		method::LaxWendroff, 
		q::Matrix{T}, 
		F::Matrix{T}, 
		dt::T, 
		h::T,
		s::State
	) where T <: AbstractFloat
		n = size(q)[2]

		# first, calculate q_star
		for j = 1:n-1
			for v = 1:3
				left = (1/2) * (q[v, j] + q[v, j+1])
				right = (dt / (2*h)) * (F[v, j+1] - F[v, j])
				method.q_star[v, j] = left - right
			end
		end

		# then calculate the F* that comes from that q*
		println("recalculating variables after first q_star calculation")
		recalculate_variables!(method.q_star, s)

		if any(isnan, method.q_star)
			error("nan found in q_star")
		end

		#println("recalculating F")
		recalculate_F!(method.F_star, s)

		if any(isnan, method.F_star)
			matches = [
			   any(isnan, method.F_star[1, :]),
			   any(isnan, method.F_star[2, :]),
			   any(isnan, method.F_star[3, :]),
			]

			error("nan found in F_star", matches)
		end

		println("finished recalculating F")

		for j = 2:n-1
			for v = 1:3
				left = q[v, j]
				# TODO: This might be j, j-1
				#right = (dt / h)  * (method.F_star[j+1] - method.F_star[j-1])
				right = (dt / h)  * (method.F_star[j] - method.F_star[j-1])
				q[v, j] = left - right
			end
		end

		if any(isnan, q)
			error("nan found in q")
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
end
