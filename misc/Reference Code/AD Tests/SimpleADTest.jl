using ForwardDiff
using Zygote: @ignore

function assembleK!(K, k)
	for i = 1:4
		for j = 1:4
			if i != j
				K[i,j] = 1.0
				continue
			end
		K[i,i] = k
		end
	end
end
function assembleF!(f)
	for i = 1:4
		f[i] = 2.0
	end
end
function assembleK(K,k)
	A = Zygote.Buffer(K)
	assembleK!(A, k)
	return copy(A)
end
function assembleF(f)
	F = Zygote.Buffer(f)
	assembleF!(F)
	return copy(F)
end


function solveCompliance(k)
	K = Matrix{Float64}(undef, 4,4)
	f = Array{Float64}(undef, 4)
	A = assembleK(K,k)
	F = assembleF(f)
	u = A \ F
	return u' * F
end

g = (k) -> ForwardDiff.derivative(solveCompliance, k)
gz = (k) -> Zygote.gradient(solveCompliance, k)

function finiteDifference(k)
	esp = 1e-6
	l = solveCompliance(k-esp)
	r = solveCompliance(k+esp)
	return (r - l) / (2 * esp)
end
