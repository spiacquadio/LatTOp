struct QuadratureRules
	nQuadPoints::Int16
	weights::Vector
	quadPoints::Matrix
end

function getQuadratureRules()
	QuadratureRules(4, [1,1,1,1], 1/sqrt(3) * [-1;-1;;  1;-1;;  1;1;;  -1;1])
end

struct ShapeFunctions
	N::Matrix
	dNdXi1::Matrix
	dNdXi2::Matrix
end

function getShapeFunctions(quad::QuadratureRules)
	xi1 = quad.quadPoints[1,:]
	xi2 = quad.quadPoints[2,:]
	shapeFuncs = [	
		(1 .- xi1) .* (1 .- xi2)/4;; 
		(1 .+ xi1) .* (1 .- xi2)/4;;
		(1 .+ xi1) .* (1 .+ xi2)/4;;
		(1 .- xi1) .* (1 .+ xi2)/4]
	return shapeFuncs
end

function getShapeDXi1(quad::QuadratureRules)
	xi2 = quad.quadPoints[2,:]
	shapeDxi1 = [
		(xi2 .- 1)/4;;
		(1 .- xi2)/4;;
		(1 .+ xi2)/4;;
	   -(1 .+ xi2)/4]
	return shapeDxi1
end

function getShapeDXi2(quad::QuadratureRules)
	xi1 = quad.quadPoints[1,:]
	shapeDxi2 = [
		(xi1 .- 1)/4;;
	   -(1 .+ xi1)/4;;
		(1 .+ xi1)/4;;
		(1 .- xi1)/4]
	return shapeDxi2
end

function getShapeAndDerivatives(quad::QuadratureRules)
	ShapeFunctions(getShapeFunctions(quad), getShapeDXi1(quad), getShapeDXi2(quad))
end

function getBmatrix!(B::Matrix, dNdX::Vector, dNdY::Vector)
	for i = 1:4
		B[1,2*i-1]	= dNdX[i];
		B[2,2*i]	= dNdY[i];
		B[3,2*i-1]	= dNdY[i];
		B[3,2*i]	= dNdX[i];
	end
end

function getJacobian(iQuad, shapeFuncs::ShapeFunctions, elePoints)
	dNdXi1 = shapeFuncs.dNdXi1[iQuad,:]'
	dNdXi2 = shapeFuncs.dNdXi2[iQuad,:]'
	jacobian = [
		dNdXi1 * elePoints[:,1]; dNdXi1 * elePoints[:,2];;
		dNdXi2 * elePoints[:,1]; dNdXi2 * elePoints[:,2]]

	detJ = (jacobian[1,1] * jacobian[2,2]) - (jacobian[2,1] * jacobian[1,2])
	dNdX = jacobian \ [dNdXi1; dNdXi2]
	return detJ, dNdX
end