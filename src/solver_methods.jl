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

	mutable struct LaxWendroff{T}
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

		return LaxWendroff(q_star, F_star)#, copy(q_star), copy(F_star))
	end

	function create_method(marker::LaxWendroffViscousFluxesMarker, n::Int, 
		s::State{T}
	)::LaxWendroff where T <: AbstractFloat
		error!("not ready yet")
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

		q[:, 2:n-1] = (1/2) * (q[:, 3:n] + q[:, 1:n-2]) - (dt / (2h)) *( F[:, 3:n] - F[:, 1:n-2])
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

		for j = 2:n-1
			method.q_star[:, j]  = ((1/2) * (q[:, j]+ q[:, j+1])) -   ((dt / (2h)) * (F[:, j+1 ] - F[:, j]))
		end

		# recalculate state variables for j+1/2 terms, then calculate F* j+1/2
		# then calculate the F* for the new state `s`
		recalculate_variables!(method.q_star, s)
		recalculate_F!(method.F_star, s)

		for j = 2:n-1
			# q_next = q[j] - Δt / h (F*[j+1/2] - F*[j-1/2])
			# since F*[j-1/2] == F*[j+1/2 - 1] (barring j =1 which we hand wave) 
			# the indexing becomes just
			# q_next = q[j] - Δt / h (F*[j] - F*[j-1])
			q[:, j] = q[:, j] - ((dt / h) .* (method.F_star[:, j] -  method.F_star[:, j-1]))
		end
		println("finish u calc")
	end

	function solver_step!(method::LaxWendroffViscousFluxes, 
		q::Matrix{T}, 
		F::Matrix{T}, 
		dt::T, 
		h::T,
		s::State
	) where T <: AbstractFloat
		error("not implemented")
	end
end
