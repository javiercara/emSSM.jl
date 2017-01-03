function dlyap(A,Q)

	# Solving Discrete LYAPunov equation using Silvester equation
	#
	# Sylvester equation A*X + X*B + C = 0
	# if B = -(A')^(-1), C = Q*(A')^(-1) 
	# => A*X - X*(A')^(-1) + Q*(A')^(-1) = 0 
	# => A*X*A' - X + Q = 0  -> Discrete Lyapunov equation

	Atinv = inv(A')
	B = -Atinv
	C = Q*Atinv
	X = sylvester(A, B, C)
	
	return X
	
end

# ======================================================
function dlyap_test()
	
	# example from 
	# Barraud, A.Y., "A numerical algorithm to solve A'XA - X = Q",
	# IEEE Trans. Auto. Contr., AC-22, pp. 883-885, 1977. 

	A = [1.5 2.5 3.5 4.5 5.5 6.5 7.5;
		1 -1.5 2.5 3.5 4.5 5.5 6.5;
		0 0 1.5 -2.5 3.5 4.5 5.5;
		0 0 0 1.5 2.5 -3.5 4.5;
		0 0 0 1 1.5 2.5 -3.5;
		0 0 0 0 0 1.5 -2.5;
		0 0 0 0 0 1 1.5]

	q = [5.25 1.5 19.25 22.5 55.75 60 49.25;
		0 0 3.25 -0.25 -6 -11 11.25;
		0 0 62.75 83 183.25 222.75 129.75;
		0 0 0 56.5 208.25 232.25 250.75;
		0 0 0 0 532.25 739.5 476.75;
		0 0 0 0 0 631 745.75;
		0 0 0 0 0 0 133.75]

	Q = (q + q') - diagm(diag(q))

	x = [1 1 1 1 1 1 1;
		0 1 2 3 4 5 6;
		0 0 1 3 6 10 15;
		0 0 0 1 4 10 20;
		0 0 0 0 1 5 15;
		0 0 0 0 0 1 6;
		0 0 0 0 0 0 1];

	X = (x + x') - diagm(diag(x))
	
	A = A'
	Q = -Q

	X1 = dlyap(A,Q)

	return X,X1

end
