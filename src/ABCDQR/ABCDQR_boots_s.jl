function ABCDQR_boots_s(y,u,A,B,C,D,Q,R,x10,nb)
	#
	# computes nb STATIONARY bootstrap replications for state space model
	#
	# x_{t+1} = A*x_{t} + B*u_{t} + w_{t}
	# y_{t}   = C*x_{t} + D*u_{t} + v_{t}
	#
	# javier.cara@upm.es, 2016-04
	
	ny,nt = size(y)
	nx = size(A,1)
	
	# step 0: compute the innovations
	# -----------------------------------------------------
	xtt,Ptt,xtt1,Ptt1,et,St,Kt,loglik = ABCDQR_kfilter_s(y,u,A,B,C,D,Q,R,x10)
	
	# step 1: standardized innovations
	# -----------------------------------------------------
	# St^(1/2)
	Std,V = eig(St)
	St_05 = V*diagm(sqrt(Std))*V' # St is symmetric = inv(V)=V'
	# St^(-1/2)*et
	et_st = St_05\et
   
	# bootstrap replications
	ybo = zeros(ny,nt,nb)
	xtt1_bo = zeros(nx,nt)   
	for i in 1:nb   	
		# step 2: sampling et_st WITH replacement
   	# ------------------------------------------------------
   	orden = rand(1:nt,nt)
   	et_st_bo = et_st[:,orden]
   	
   	# step 3: bootstrap data set
   	# ------------------------------------------------------
		for t in 1:nt-1
			xtt1_bo[:,t+1] = A*xtt1_bo[:,t] + B*u[:,t] + A*Kt*St_05*et_st_bo[:,t]
		end
		ybo[:,:,i] = C*xtt1_bo + D*u + St_05*et_st_bo
   	
	end # nb
   
	return ybo
   
end
   	

