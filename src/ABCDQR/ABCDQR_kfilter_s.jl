function ABCDQR_kfilter_s(y,u,A,B,C,D,Q,R,x10)
	# 
	# STATIONARY Kalman filter for model
	# 
	# x_{t+1} = A*x_{t} + B*u_{t} + w_{t}
	# y_{t}   = C*x_{t} + D*u_{t} + v_{t}
	# 
	# cov(w_{t},v_{t}) = [Q 0;0 R]
	#
	# javier.cara@upm.es, 2015-10 
	# 

	(ny,nt) = size(y)
	nx = size(A,1)	

	# allocation
	xtt = zeros(nx,nt)
	xtt1 = zeros(nx,nt+1)
	et = zeros(ny,nt)
	loglik = 0.0

	# Discrete-time Algebraic Riccati Equation (DARE) solution
	S = zeros(nx,ny)
	Ptt1 = dare(A',C',Q,R,S)

	# variance of innovations
	St = C*Ptt1*C' + R 
	Stinv = eye(ny,ny)/St # numerically preferible to Stinv=inv(St)
	
	# Kalman gain matrix
	Kt = Ptt1*C'*Stinv
	
	# Var[x_t|y_t]
	Ptt = (eye(nx)-Kt*C)*Ptt1

	# Filter 
	xtt1[:,1] = x10
	for t in 1:nt		
		#  innovations
		et[:,t] = y[:,t] - C*xtt1[:,t] - D*u[:,t]
	
		# filtered values
		xtt[:,t] = xtt1[:,t] + Kt*et[:,t]
		
		# one-step ahead prediction
		xtt1[:,t+1] = A*xtt[:,t] + B*u[:,t]
		
		# likelihood
		l0 = et[:,t]'*Stinv*et[:,t] # typeof(l0) = Array{Float64,1}
		loglik = loglik + l0[1]	
	end

	loglik = - ny*nt/2*log(2*pi) - nt/2*log(det(St)) - 1/2*loglik

	return xtt,Ptt,xtt1,Ptt1,et,St,Kt,loglik

end

