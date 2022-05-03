module Methods
	using ..StateContainers: State, check_state
	using ..Calculations: recalculate_q!, recalculate_F!, recalculate_variables!, adjust_Fstar!
	using Printf

	export MethodMarker
	export LaxFriedrichs, LaxFriedrichsMarker
	export LaxWendroff, LaxWendroffMarker
	export LaxWendroffViscousFluxes, LaxWendroffViscousFluxesMarker
	export create_method
	export solver_step!
	export unique_name

	abstract type MethodMarker end
	struct LaxFriedrichsMarker <: MethodMarker end
	struct LaxWendroffMarker <: MethodMarker end
	struct LaxWendroffViscousFluxesMarker <: MethodMarker end

	struct LaxFriedrichs end

	mutable struct LaxWendroff{T}
		q_star::Matrix{T}
		F_star::Matrix{T}
	end

	mutable struct LaxWendroffViscousFluxes{T}
		q_star::Matrix{T}
		F_star::Matrix{T}
	end

	function create_method(marker::LaxFriedrichsMarker, n::Int, _...)::LaxFriedrichs 
		return LaxFriedrichs()
	end

	function unique_name(marker::LaxFriedrichsMarker, alpha::Float64, n::Int)::String

		return "flowfield_" * "lax_friedrichs_n_" * string(n) * ".h5"
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

	function unique_name(marker::LaxWendroffMarker,  _...)::String
		return "flowfield_" * "lax_wendroff.h5"
	end

	function create_method(marker::LaxWendroffViscousFluxesMarker, n::Int, 
		s::State{T}
	)::LaxWendroffViscousFluxes where T <: AbstractFloat
		q_star = zeros(3, n)
		F_star = zeros(3, n)

		recalculate_q!(q_star, s)
		recalculate_F!(F_star, s)

		return LaxWendroffViscousFluxes(q_star, F_star)
	end

	function unique_name(marker::LaxWendroffViscousFluxesMarker,  alpha::Float64, n::Int)::String
		alpha = @sprintf "%2.2f" alpha
		return "flowfield_" * "lax_wendroff_fluxes_alpha_" *alpha * ".h5"
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

		for j = 1:n-1
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
	end

	function solver_step!(method::LaxWendroffViscousFluxes, 
		q::Matrix{T}, 
		F::Matrix{T}, 
		dt::T, 
		h::T,
		s::State
	) where T <: AbstractFloat
		n = size(q)[2]

		for j = 1:n-1
			method.q_star[:, j]  = ((1/2) * (q[:, j]+ q[:, j+1])) -   ((dt / (2h)) * (F[:, j+1 ] - F[:, j]))
		end

		# recalculate state variables for j+1/2 terms, then calculate F* j+1/2
		# then calculate the F* for the new state `s`
		recalculate_variables!(method.q_star, s)
		recalculate_F!(method.F_star, s)
		method.F_star = adjust_Fstar!(method.F_star, s)

		for j = 2:n-1
			# q_next = q[j] - Δt / h (F*[j+1/2] - F*[j-1/2])
			# since F*[j-1/2] == F*[j+1/2 - 1] (barring j =1 which we hand wave) 
			# the indexing becomes just
			# q_next = q[j] - Δt / h (F*[j] - F*[j-1])
			q[:, j] = q[:, j] - ((dt / h) .* (method.F_star[:, j] -  method.F_star[:, j-1]))
		end
	end
end
