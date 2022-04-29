module Calculations
	using ..StateContainers: State, check_state
	export recalculate_variables!, recalculate_q!, recalculate_F!

	function recalculate_variables!(
		q::Matrix{T},
		s::State{T},
	) where T <: AbstractFloat
		s.ρ .= q[1, :]
		s.u .= q[2, :] ./ s.ρ
		s.E .= q[3, :] ./ s.ρ
		s.e .= s.E - ((1/2) .* s.u.^2)
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

	function adjust_Fstar!(
		F::Matrix{T}, 
		s::State,
	) where T <: AbstractFloat
		n = size(s.u)[2]

		du_dx_buffer = copy(s.u)
		forward_difference(s.u, du_dx_buffer)
		error!("we should not be here, we should have errored on the function call above")

		vec = zeros(3,n)
		vec[2] .= 1
		vec[3] = s.u

		term = s.alpha * s.h^2 .* s.rho .* abs(du_dx_buffer) .* du_dx_buffer .* vec 

		F -= term
	end

	function forward_difference!(
		input::Vector{T},
		output::Vector{T},
		h::T
	) where T <: AbstractFloat
		for i = 1:length(input)-1
			output[i] = (input[i+1] - input[i]) / h
		end
	end
end
