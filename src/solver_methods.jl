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
		q_star_plus_half::Matrix{T}
		F_star_plus_half::Matrix{T}
		q_star_minus_half::Matrix{T}
		F_star_minus_half::Matrix{T}
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

		return LaxWendroff(q_star, F_star, copy(q_star), copy(F_star))
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

		#for j = 2:n-1
		#	for v = 1:3
		#		left = (1/2) * (q[v, j+1] + q[v, j-1])
		#		right = (dt / (2*h)) * (F[v, j+1] - F[v, j-1])
		#		q[v, j] = left - right
		#	end
		#end
		
		#for j = 2:n-1
		#	left = (1/2) * (q[:, j+1] + q[:, j-1])
		#	right = (dt / (2*h)) * (F[:, j+1] - F[:, j-1])
		#	q[:, j] = left[:] - right[:]
		#end

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

		# # first, calculate q_star
		# for j = 2:n-1
		# 	for v = 1:3
		# 		q_left = q[v, j-1]
		# 		q_mid = q[v, j]
		# 		q_right = q[v,j+1]

		# 		F_left = F[v, j-1]
		# 		F_mid = F[v, j]
		# 		F_right = F[v,j+1]

		# 		method.q_star_plus_half[v, j] = ((1/2) * (q_mid + q_right)) - ((dt / (2*h)) * (F_right- F_mid))
		# 		method.q_star_minus_half[v, j] = ((1/2) * (q_left + q_mid)) - ((dt / (2*h)) * (F_mid - F_left))
		# 	end
		# end

		# (1/2) ( q[j] + q[j+1] ) - (dt / 2h) * (F[j+1] - F[j])
		method.q_star_plus_half[:, 2:n-1]  = ((1/2) * (q[:, 2:n-1]+ q[:, 3:n])) -   ((dt / (2h)) * (F[:, 3:n  ] - F[:, 2:n-1]))
		# (1/2) ( q[j] + q[j-1] ) - (dt / 2h) * (F[j] - F[j-1])
		method.q_star_minus_half[:, 2:n-1] = ((1/2) * (q[:, 2:n-1]+ q[:, 1:n-2])) - ((dt / (2h)) * (F[:, 2:n-1] - F[:, 1:n-2]))

		# recalculate variables for +1/2 terms, then calculate F* j+1/2
		# then calculate the F* for the new state `s`
		recalculate_variables!(method.q_star_plus_half, s)
		recalculate_F!(method.F_star_plus_half, s)

		# recalculate variables for -1/2 terms, then calculate F* j-1/2
		# then calculate the F* for the new state `s`
		recalculate_variables!(method.q_star_minus_half, s)
		recalculate_F!(method.F_star_minus_half, s)

		q[:, 2:n-1] = q[:, 2:n-1] - ((dt / h) .* (method.F_star_plus_half[:, 2:n-1] -  method.F_star_minus_half[:, 2:n-1]))

		#q[:, 2:n-1] = q[:, 2:n-1] - ((dt / h) * (method.F_star_plus_half[:, 2:n-1] -  method.F_star_plus_half[:, 1:n-2]))

		#for j = 2:n-1
		#	for v = 1:3
		#		left = q[v, j]
		#		right = (dt / h)  * (method.F_star_plus_half[v, j+1] - method.F_star_minus_half[v, j-1])
		#		q[v, j] = left - right
		#	end
		#end

		println("finish u calc")
	end

	function solver_step!(method::LaxWendroffViscousFluxes, 
		q::Matrix{T}, 
		F::Matrix{T}, 
		dt::T, 
		h::T,
		s::State
	) where T <: AbstractFloat

	end
end
