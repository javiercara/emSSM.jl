function byrow(y)
	# data by rows
	# if vector, converts to a rowmatrix
	
	if ndims(y) == 1
		# vector
		nt = length(y)
		ny = 1
		y = reshape(y,1,nt)
	else
		# matrix
		(ny,nt) = size(y)	
		
		if nt < ny
			y=y'
			(ny,nt) = size(y)
		end
	end
	
	return y,ny,nt
end
