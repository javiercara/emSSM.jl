function ACQR_bench(bench,nt,seed=1)
	#
	# benchmark models to test the functions
	#
	
	# random seed for the simulations
	srand(seed)
    
	if bench==1
		# 
		A = 0.6
		C = 1
		Q = 2
		R = 4
		
		v = sqrt(R)*randn(nt)
		w = sqrt(Q)*randn(nt)
		
		x = zeros(nt)
		for t in 1:nt-1
			x[t+1] = A*x[t] + w[t]
		end
		y = C*x + v
		
		salida = y,A,C,Q,R
		
	elseif bench==2
		#
		# Based on: 
		# "Subspace Identification for Linear Systems"
		# P. Van Overschee and B. De Moor, page 81
		#
		A = [0.6 0.6 0.;-0.6 0.6 0.;0. 0. 0.4]
		C = [0.7831 0.5351 0.9701]
		R = 6.0
		Q = [3 0.5 0.2;0.5 2 0.3;0.2 0.3 1]
		U = chol(Q)
        
		v = sqrt(R)*randn(1,nt)
		ns = size(A,1)
		w = U'*randn(ns,nt)
        
 		x = zeros(ns,nt)
		x[:,1] = w[:,1]
		for t in 1:nt-1
			x[:,t+1] = A*x[:,t] + w[:,t]
		end
		y = C*x + v
        
		salida = y,A,C,Q,R
		  
	end
    
	return salida
    
end
