module SodsShockTube

	greet() = print("Hello World!")

	function main()
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
		t_end = 3.9 / 1000

		#
		# Initial conditions
		#

		# these are given in the problem statelemt
		# from the paper
		uL = ;
		uR = ;
		pL = ;
		pR = ;
		rhoL = ;
		rhoR = ;
		# form the paper p = (gamma - 1) rho e
		# so you can calculate these
		eL = ;
		eR = ;

		# given in the paper again
		cL ~= sqrt(gamma p / rho)
		cL ~= sqrt(gamma p / rho)
		

		#
		# states / variables
		#
		rho = zeros(n)
		p = zeros(n)
		u = zeros(n)
		e = zeros(n)
		E = zeros(n) # e + 1/2 u^2
		c = zeros(n)

		# the discontinuity should be exactly in the middle
		half_n = Int(floor((n-1) / 2));
		rho[1:half_n] = rhoL;
		rho[half_n:n] = rhoR;

		p[1:half_n] = pL;
		p[half_n:n] = pR;

		# TODO: repeat this exercise for u, e, E, c, etc
		
		#
		#
		#

		# TODO: calculate what umax is
		# it could be the speed of sound (c) or speed of the 
		# wave (u?)
		dt = CFL * h / umax

		while t <= t_end
			# construct q vector according to paper
			# row is the variable, the column is the 
			# spatial values
			q = [
				rho;
				rho .* u;
				rho .* E
			]

			# construct F vector
			F = [
				rho .* u;
				rho .* u.^2;
				u .* (rho .* (e .+ (1/2) .* u.^2) .+ p)
			]

			if flag == 1
				# TODO: Lax Freidrichs method
				# according to the paper
				q[1:3, 2:n-1] = todo
			elseif flag == 2
				# since we need j+1 we can only go to n-1
				qstar[1:3, 1:n-1] = todo
				# here the input to the flux is qstar
				# (instead of q)
				Fstar(1:3, 1:n-1] = 

				# then do the second step of the equations
				q[1:3, 1:n-1] = 
			elseif flag == 3
			end

			# then we go back to calculate what the following
			# variables are for the next time step
			rho = q[1, 1:n];
			u = ;
			E = ;
			e = ;
			p = ;
		end


		# then, do the plotting below
	end

end # module
