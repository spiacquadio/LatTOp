include("shapeFunctions.jl")

function computeStrain(u, ncon, pos, quad::QuadratureRules, shape::ShapeFunctions)
	nElements = size(ncon)[1]

	strain = zeros(nElements, 3)
	uElement = Vector{Float64}(undef, 8)
	B = zeros(3,8)
	for e = 1:nElements
		element = pos[ncon[e,:],:]
		for i = 1:4
			uElement[2*i-1] = u[2*ncon[e,i]-1]
			uElement[2*i] = u[2*ncon[e,i]]
		end

		for q = 1:quad.nQuadPoints
			detJ, dNdX = getJacobian(q, shape, element)
			getBmatrix!(B,dNdX[1,:], dNdX[2,:])
			for i = 1:3
				for j = 1:8
					strain[e,i] += B[i,j] * uElement[j]
				end
			end
		end
	end

	return strain;
end;