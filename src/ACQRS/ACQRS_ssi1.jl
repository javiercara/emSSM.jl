function ACQRS_ssi1(y,n,i)
	#
	# ssi algorithm for the identification of state space models
	# -- CVA version --
	#
	# Reference:
	# Subspace Identification for Linear Systems
	# Theory - Implementation - Applications
	# Peter Van Overschee / Bart De Moor
	# Kluwer Academic Publishers, 1996
	#
	# javier.cara@upm.es, 2015-08
	# --------------------------------------------------
	
	# Data by rows
	y,m,N = byrow(y)
	
	# Hankel matrix
	# --------------------------------------------------
	# number of Hankel matrix columns
	j = N - 2*i + 1
	if (j < m*2*i) 
		# rows(H) has to be > columns(H)
		error("Not enough data for building the Hankel matrix")
	end
	Y = hankel_yij(y/sqrt(j),2*i,j)
	
	# LQ factorization 
	# --------------------------------------------------
	F = qrfact(Y')
	L = F[:R]'
	L21 = L[m*i+1:2*m*i,1:m*i]
	
	# singular values
	# --------------------------------------------------
	F = svdfact(L21)
	U1 = F[:U][:,1:n]
	S1 = F[:S][1:n]
	
	# Matrices gam and gam1
	# --------------------------------------------------
	gam  = U1*diagm(sqrt(S1))
	gam1 = U1[1:m*(i-1),:]*diagm(sqrt(S1))
	# and pseudo-inverses
	gam_inv  = pinv(gam)
	gam1_inv = pinv(gam1)
	
	# Determine the states Xi and Xi1
	# --------------------------------------------------
	Xi  = gam_inv  * L21
	Xi1 = gam1_inv * L[m*(i+1)+1:2*m*i,1:m*(i+1)]
	
	# Computing the state matrices A and C
	# --------------------------------------------------
	Rhs = [ Xi zeros(n,m) ]	# Right hand side
	Lhs = [ Xi1; L[m*i+1:m*(i+1),1:m*(i+1)] ] # Left hand side
	
	# Least squares
	sol = Lhs/Rhs

	A = sol[1:n,1:n]
	C = sol[n+1:n+m,1:n]
	
	# Computing the covariance matrices Q, R and S
	# -------------------------------------------------
	# Residuals
	res = Lhs - sol*Rhs
	cov = res*res'
	
	Q = cov[1:n,1:n]
	S = cov[1:n,n+1:n+m]
	R = cov[n+1:n+m,n+1:n+m]
	
	return A,C,Q,R,S
	
end

##########################################################
function hankel_yij(y,i,j)
	# Hankel matrix with i rows and j columns
	# from data [y_1 y_2 ... y_N],  y_t in R^{m x 1}

	y,m,N = byrow(y)

	# Hankel matrix
	H = zeros(m*i,j)
	for k in 1:i
		H[(k-1)*m+1:k*m,:] = y[:,k:k+j-1]
	end

	return H

end
