function simsys(sys,nt,seed)
	###
	# simulated systems
	###
	
	# random seed for the simulations
	srand(seed)
    
	if sys==1
		#
		# system 1
	
		# ref: "Subspace Identification for Linear Systems"
		# P. Van Overschee and B. De Moor, page 81
		#
		A = [0.6 0.6 0.;-0.6 0.6 0.;0. 0. 0.4]
		K = [0.1706 -0.1507 0.2772]'
		C = [0.7831 0.5351 0.9701]
		R = 6.3663
		Q = K * R * K'
        
		et = sqrt(R)*randn(1,nt)
		w = K*et
        
		ns = 3
 		x = zeros(ns,nt)
		x[:,1] = w[:,1]
		for t in 1:nt-1
			x[:,t+1] = A * x[:,t] + w[:,t]
		end
		y = C * x + et
        
		salida = y, A, K, C, Q, R
		#
	elseif sys==2
		###
		# ref: "Subspace Identification for Linear Systems"
		# P. Van Overschee and B. De Moor, page 115
		###
		A = [0.6 0.6 0.;-0.6 0.6 0.;0. 0. 0.7]
		B = [1.6161 -0.3481 2.6319]'
		K = [-1.1472 -1.5204 -3.1993]'
		C = [-0.4373 -0.5046 -0.0936]
		D = -0.7759
		R = 0.0432
		Q = K * R * K'
        
		et = sqrt(R)*randn(1,nt)
		w = K*et
        
		# input = AR(1)
		u = zeros(1,nt)
		a = randn(nt)
		for t in 1:nt-1
			u[t+1] = 0.3*u[t] + a[t]
		end
        
		# states
		ns = size(A,1)
		x = zeros(ns,nt)
		x[:,1] = w[:,1]
		for t in 1:nt-1
			x[:,t+1:t+1] = A * x[:,t:t] + B * u[t] + w[:,t:t]
		end
		# output
		y = C * x + D * u + et
        
		salida = y, u, A, B, K, C, D, Q, R        
	end
    
	return salida
    
end
