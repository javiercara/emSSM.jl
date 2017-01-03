function aicSSM(m)
	# Akaike Information Criterion for State Space Models
	nx = size(m.A,1)
	ny = size
	

	P = ny*nx + nx*nx + nx*(nx+1)/2 + ny*(ny+1)/2 + nx*nu + ny*nu
	aic[k] = -2*loglik[end] + 2*P



end
