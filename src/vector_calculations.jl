module Calculations
	using ..StateContainers: State, check_state
	export recalculate_variables!, recalculate_q!, recalculate_F!

	function recalculate_variables!(
		q::Matrix{T},
		s::State{T};
		update_speed_of_sound::Bool=false
	) where T <: AbstractFloat
		s.ρ .= q[1, :]
		s.u .= q[2, :] ./ s.ρ
		s.E .= q[3, :] ./ s.ρ
		s.e .= s.E - ((1/2) .* s.u.^2)
		s.p .= (s.gamma - 1) .* s.ρ .* s.e

		if update_speed_of_sound 
			s.c .= sqrt.(s.gamma .* s.p ./ s.ρ)
		end
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
		n = length(s.u)

		du_dx_buffer = copy(s.u)
		forward_difference!(s.u, du_dx_buffer, s.h)

		vec = zeros(3,n)
		vec[2, :] .= 1
		vec[3, :] = s.u

		term = s.alpha * s.h^2 .* s.ρ .* abs.(du_dx_buffer) .* du_dx_buffer .* transpose(vec[:, 1:n])

		F -= transpose(term)
		return F
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
