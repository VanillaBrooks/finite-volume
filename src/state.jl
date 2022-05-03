module StateContainers
	export State, create_state, InitialConditions, check_state

	struct InitialConditions{T}
		uR::T
		uL::T
		#
		pR::T
		pL::T
		#
		ρR::T
		ρL::T
		#
		eR::T
		eL::T
		#
		cR::T
		cL::T
	end

	function create_ic(gamma::T)::InitialConditions{T} where T <: AbstractFloat
		# these are given in the problem statelemt
		# from the paper
		uL = 0.;
		uR = 0.;
		pL = 10.0^5;
		pR = 10.0^3;
		ρL  = 1.0;
		ρR = 0.01;
		# form the paper p = (gamma - 1) rho e
		# so you can calculate these
		eL = pL / ((gamma - 1) * ρL);
		eR = pR / ((gamma - 1) * ρR);

		# given in the paper again
		cL = sqrt(gamma * pL/ ρL)
		cR = sqrt(gamma * pR / ρR)

		return InitialConditions(uR, uL, pR, pL, ρR, ρL, eR, eL, cR, cL)
	end

	mutable struct State{T}
		gamma::T
		alpha::T
		h::T
		ic::InitialConditions
		u::Vector{T}
		ρ::Vector{T}
		e::Vector{T}
		c::Vector{T}
		E::Vector{T}
		p::Vector{T}
	end

	function create_state(n::Int, h::T, alpha::T)::State{T} where T <: AbstractFloat
		gamma = 1.4

		ic = create_ic(gamma)

		ρ = zeros(n)
		p = zeros(n)
		u = zeros(n)
		e = zeros(n)
		E = zeros(n) # e + 1/2 u^2
		c = zeros(n)

		# the discontinuity should be exactly in the middle
		half_n = Int(floor((n+1) / 2));

		ρ[1:half_n] .= ic.ρL;
		ρ[half_n+1:n] .= ic.ρR;

		p[1:half_n] .= ic.pL;
		p[half_n+1:n] .= ic.pR;

		u[1:half_n] .= ic.uL;
		u[half_n+1:n] .= ic.uR;

		e[1:half_n] .= ic.eL;
		e[half_n+1:n] .= ic.eR;

		E[1:half_n] .= ic.eL + (1/2) .* ic.uL.^2;
		E[half_n+1:n] .= ic.eR + (1/2) .* ic.uR.^2;

		c[1:half_n] .= ic.cL
		c[half_n+1:n] .= ic.cR

		return State(
			gamma,
			alpha,
			h,
			ic,
			u,
			ρ,
			e,
			c,
			E,
			p
		)
	end

	function check_state(s::State)
		if any(isnan, s.ρ)
			error("Nan exists in rho")
		end
		if any(iszero, s.ρ)
			error("zero exists in rho")
		end
		if length(filter(x -> x < 0, s.ρ)) > 0
			display(transpose(s.ρ))
			error("less than zero exists in rho")
		end
		if length(filter(x -> x < 0, s.E)) > 0
			display(transpose(s.E))
			error("less than zero exists in E")
		end
		if length(filter(x -> x < 0, s.e)) > 0
			display(transpose(s.p))
			error("less than zero exists in e")
		end
		if length(filter(x -> x < 0, s.p)) > 0
			display(transpose(s.p))
			error("less than zero exists in p")
		end
		if any(isnan, s.u)
			error("Nan exists in u")
		end
		if any(isnan, s.E)
			error("Nan exists in E")
		end
		if any(isnan, s.e)
			error("Nan exists in e")
		end
		if any(isnan, s.p)
			error("Nan exists in p")
		end
	end
end
