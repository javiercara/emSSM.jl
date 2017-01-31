__precompile__()


module emSSM

#######################################################
# types
#######################################################

typealias ScalarOrArray{T} Union{T, Array{T}}

type SSM
	# algorithms restricted to matrices
	A::Array{Float64,2}
	B::Array{Float64,2}
	C::Array{Float64,2}
	D::Array{Float64,2}
	Q::Array{Float64,2}
	R::Array{Float64,2}
	S::Array{Float64,2}
	m1::Array{Float64,1}
	P1::Array{Float64,2}
	loglik::Float64
	aic::Float64
end

function SSM0()
	# initialize SSM type with zeros
	m = zeros(2,2) # matrix of zeros
	v = zeros(2) # vector of zeros
	A = SSM(m,m,m,m,m,m,m,v,m,0.,0.)
	
	return A
	
end

export SSM, SSM0

#####################################################
# em
#####################################################

include("EM/em.jl")

export em

#####################################################
# auxfun
#####################################################

include("auxfun/byrow.jl")
include("auxfun/dare.jl")
include("auxfun/dlyap.jl")

export byrow, dare, dlyap
	
####################################################
# ACQR
####################################################

include("ACQR/ACQR_kfilter_s.jl")
include("ACQR/ACQR_ksmoother_s.jl")
include("ACQR/ACQR_em_s.jl")
include("ACQR/ACQR_bench.jl")
include("ACQR/ACQR_boots_s.jl")

export
	ACQR_kfilter_s,
	ACQR_ksmoother_s,
	ACQR_em_s,
	ACQR_bench,
	ACQR_boots_s
	
####################################################
# ACQRS
####################################################

include("ACQRS/ACQRS_ssi1.jl")

export ACQRS_ssi1
	
####################################################
# ABCDQR
####################################################
	
include("ABCDQR/ABCDQR_kfilter_s.jl")
include("ABCDQR/ABCDQR_em_s.jl")
include("ABCDQR/ABCDQR_bench.jl")
include("ABCDQR/ABCDQR_boots_s.jl")

export 
	ABCDQR_kfilter_s,
	ABCDQR_em_s,
	ABCDQR_bench,
	ABCDQR_boots_s
	

end # module
