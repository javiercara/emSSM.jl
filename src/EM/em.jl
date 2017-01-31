function em(y;nx::Int=2,u::Array=zeros(2),init=false,ssmi=SSM0(),max_iter::Int=100,tol=1e-6,txo=true,subid=false)
	"""
	estimate A, B, C, D, Q, R, S, m1, P1 using the STATIONARY EM algorithm for model
	
	x_{t+1} = A*x_{t} + B*u_{t} + w_{t}
	y_{t}   = C*x_{t} + D*u_{t} + v_{t}
	
	cov(w_{t},v_{t}) = [Q S;S' R]
	x1 -> N(m1,P1)	
	
	Input
	----------------------
		y          : data
		nx         : size of the estimated A matrix
		u != 0			 : input are considered
		init				: if init, ssmi is considered
		ssmi       : initial point for the EM algorithm
		subid=true : SUBspace IDentification algorithm is used, not the EM 
		S=true     : matrix S is estimated
	
	
	"""
	
	# output declaration
	ssm = SSM0()
	
	# initial point
	# ***************************************************************
	if init
		Ai = ssmi.A
		Bi = ssmi.B
		Ci = ssmi.C
		Di = ssmi.D
		Qi = ssmi.Q
		Ri = ssmi.R
		m1i = ssmi.m1
	else
		# initial point with subspace algorithm
		i = nx+1
		Ai,Ci,Qi,Ri,Si = ACQRS_ssi1(y,nx,i)	
		m1i = Ci\y[:,1] # x11 = C^{-1}(y_1 - v_1)
		
		if u != zeros(2)
			# initial values for B and D
			y,ny,nt = byrow(y)
			u,nu,nt = byrow(u)
			Bi = zeros(nx,nu)
			Di = zeros(ny,nu)
		end
	end
	
	# estimation of parameters
	# *************************************************************
	if subid # subspace algorithm
		
		ssm.A = Ai
		ssm.C = Ci
		ssm.Q = Qi
		ssm.R = Ri
		ssm.S = Si
		ssm.m1 = m1i
		
		# loglik
		xtt,Ptt,xtt1,Ptt1,et,St,Kt,loglik = ACQR_kfilter_s(y,Ai,Ci,Qi,Ri,m1i) # change to kfilter with S
		ssm.loglik = loglik
		# aic
		ny = size(Ci,1)
		P = ny*nx + nx*nx + nx*(nx+1)/2 + ny*(ny+1)/2
		ssm.aic = -2*loglik + 2*P
		
	else # em algorithm
		
		if u == zeros(2) # em without inputs
			
			# em 
			A,C,Q,R,m1,P1,loglikv,aic = ACQR_em_s(y,nx,Ai,Ci,Qi,Ri,m1i,max_iter,tol,txo)
				
			# output
			ssm.A = A
			ssm.C = C
			ssm.Q = Q
			ssm.R = R
			ssm.m1 = m1
			ssm.P1 = P1
			ssm.loglik = loglikv[end]
			ssm.aic = aic
				
		else # em with inputs
			
			# em
			A,B,C,D,Q,R,m1,P1,loglikv,aic = ABCDQR_em_s(y,u,nx,Ai,Bi,Ci,Di,Qi,Ri,m1i,max_iter,tol,txo)
		
			# output
			ssm.A = A
			ssm.B = B
			ssm.C = C
			ssm.D = D
			ssm.Q = Q
			ssm.R = R
			ssm.m1 = m1
			ssm.P1 = P1
			ssm.loglik = loglikv[end]	
			ssm.aic = aic
		
		end
	end
	
	# output
	return ssm

end
