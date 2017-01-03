function ACQR_ksmoother_s(A,xtt,Ptt,xtt1,Ptt1)
	# 
	# STATIONARY Kalman smoother for model
	# 
	# x_{t+1} = Ax_{t} + w_{t}
	# y_{t}   = Cx_{t} + v_{t}
	# 
	# cov(w_{t},v_{t}) = [Q 0;0 R]
	#
	# javier.cara@upm.es, 2016-02
	# 
	
	nx,nt = size(xtt)
	
	# allocation
	xtN = zeros(nx,nt+1)

	# values for t=nt+1
	xtN[:,nt+1] = xtt1[:,nt+1]
	
	# values for t=nt
	xtN[:,nt] = xtt[:,nt]
	
	# Kalman Smoother matrix J and transpose
	Jt = (Ptt*A') / Ptt1
	JtT = Jt'
	
	# smoother
	
	# PtN = Ptt + Jt*(Pt1N - Pt1t)*JtT =>
	# Jt*PtN*JtT - PtN + (Ptt - Jt*Pt1t*JtT) = 0 =>
	# PtN is the solution of a discrete lyapunov equation
	PtN =  dlyap(Jt,Ptt - Jt*Ptt1*JtT)
	Pt1tN = PtN*JtT
	
	for t = nt-1:-1:1	
		xtN[:,t] = xtt[:,t] + Jt*( xtN[:,t+1] - xtt1[:,t+1] )
	end	
	
	return xtN,PtN,Pt1tN
	
end


