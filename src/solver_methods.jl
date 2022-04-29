module Methods
	using ..StateContainers: State, check_state
	using ..Calculations: recalculate_q!, recalculate_F!, recalculate_variables!

	export MethodMarker
	export LaxFriedrichs, LaxFriedrichsMarker
	export LaxWendroff, LaxWendroffMarker
	export LaxWendroffViscousFluxes, LaxWendroffViscousFluxesMarker
	export create_method
	export solver_step!

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

	function create_method(marker::LaxWendroffViscousFluxesMarker, n::Int, 
		s::State{T}
	)::LaxWendroff where T <: AbstractFloat
		q_star = zeros(3, n)
		F_star = zeros(3, n)

		recalculate_q!(q_star, s)
		recalculate_F!(F_star, s)

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
		#println("recalculating variables after first q_star calculation")
		recalculate_variables!(method.q_star, s)
		recalculate_F!(method.F_star, s)
		#println("finished recalculating F")

		for j = 2:n-1
			for v = 1:3
				left = q[v, j]
				# TODO: This might be j, j-1 OR j+1 j-1
				right = (dt / h)  * (method.F_star[j] - method.F_star[j-1])
				q[v, j] = left - right
			end
		end
	end

	function solver_step!(method::LaxWendroffViscousFluxes, 
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
		#println("recalculating variables after first q_star calculation")
		recalculate_variables!(method.q_star, s)
		recalculate_F!(method.F_star, s)
		adjust_Fstar!(method.F_star, s)

		for j = 2:n-1
			for v = 1:3
				left = q[v, j]
				# TODO: This might be j, j-1 OR j+1 j-1
				right = (dt / h)  * (method.F_star[j] - method.F_star[j-1])
				q[v, j] = left - right
			end
		end
	end
end
