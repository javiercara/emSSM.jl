function dare(A,B,Q,R,S=0)
	# 
	# solves discrete-time algebraic Riccati equation
	# X = A'XA - (A'XB+S)(R+B'XB)^{-1}(B'XA+S') + Q
	#
	# R^{-1} is not computed
	#
	# javier.cara@upm.es, 2016-02 
	# 
	
	# Dimensions
	n,m = size(B)
	
	if S== 0
		S = zeros(n,m)
	end

	L = [A zeros(n,n) B;Q -eye(n) S;S' zeros(m,n) R]
	M = [eye(n) zeros(n,n) zeros(n,m); zeros(n,n) -A' zeros(n,m);zeros(m,n) -B' zeros(m,m)]

	# Compute the eigenvalue decomposition
	(d,v) = eig(L,M)
	
	# Sort the eigenvalues
	ew = abs(d)
	pos = sortperm(ew)

	# Compute X
	v1 = v[1:n,pos[1:n]]
	v2 = v[n+1:2*n,pos[1:n]]
	X = real(v2/v1)
	
	return X

end

######
function dare_test()
    #=
    BibText
    @TECHREPORT{Benner95acollection,
    author = {Peter Benner and Alan J. Laub and Volker Mehrmann},
    title = {A Collection of Benchmark Examples for the Numerical Solution of
    Algebraic Riccati Equations II: Discrete-Time Case},
    institution = {FAK. F. MATHEMATIK, TU CHEMNITZ--ZWICKAU},
    year = {1995}
    }
    
    nov-2016
    =#

	# ------------------------------------------------------------------
	# ejm1 collection of benchmark examples
	A1 = [4 3;-4.5 -3.5]
	B1 = [1 -1]'
	R1 = 1
	Q1 = [9 6;6 4]
	X1 = (1+sqrt(5))/2*[9 6;6 4]

	X11 = dare(A1,B1,Q1,R1)
    
	# --------------------------------------------------------------------
	# ejm computed with mathematica
	Am = [1 -1 1;0 1 1;0 0 1]
	Bm = [1 0;1 0;0 1]
	Rm = [10 0;0 0.1]
	Qm = [10 0 0;0 1 0;0 0 0.1]
	Xm = [42.2835 -68.5247 -3.94783;-68.5247 154.043 16.0017;-3.94783 16.0017 8.33197]
    
	Xmm = dare(Am,Bm,Qm,Rm)
    
	#####
    
	return X1,X11,Xm,Xmm
    
end
