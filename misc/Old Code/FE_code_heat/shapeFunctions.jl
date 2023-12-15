struct QuadratureRules
	nQuadPoints::Int16
	weights::Vector
	quadPoints::Matrix
end;

@enum ElementType begin
	Quad1 = 1;
	Lin1 = 2;
end;

function getQuadratureRules(type::ElementType)
	if type == Quad1
		QuadratureRules(4, [1,1,1,1], 1/sqrt(3) * [-1;-1;;  1;-1;;  1;1;;  -1;1])
	elseif type == Lin1
		QuadratureRules(1, [2], [0;;])
	end
end;

struct ShapeFunctions
	N::Matrix
	dNdXi1::Matrix
	dNdXi2::Matrix
end;

struct ShapeFunctions1D
	N::Matrix
	dNdXi::Matrix
end;

function getShapeFunctions(type::ElementType, quad::QuadratureRules)
	if type == Quad1
		xi1 = quad.quadPoints[1,:]
		xi2 = quad.quadPoints[2,:]
		shapeFuncs = [	
			(1 .- xi1) .* (1 .- xi2)/4;; 
			(1 .+ xi1) .* (1 .- xi2)/4;;
			(1 .+ xi1) .* (1 .+ xi2)/4;;
			(1 .- xi1) .* (1 .+ xi2)/4]
		return shapeFuncs
	elseif type == Lin1
		xi = quad.quadPoints
		shapeFuncs = 0.5 * [1 .- xi;; xi .- 1]
		return shapeFuncs
	end
end;

function getShapeDXi1(quad::QuadratureRules)
	xi2 = quad.quadPoints[2,:]
	shapeDxi1 = [
		(xi2 .- 1)/4;;
		(1 .- xi2)/4;;
		(1 .+ xi2)/4;;
	   -(1 .+ xi2)/4]
	return shapeDxi1
end;

function getShapeDXi2(quad::QuadratureRules)
	xi1 = quad.quadPoints[1,:]
	shapeDxi1 = [
		(xi1 .- 1)/4;;
	   -(1 .+ xi1)/4;;
		(1 .+ xi1)/4;;
		(1 .- xi1)/4]
	return shapeDxi1
end;

function getShapeDxi(quad::QuadratureRules)
	dNdXi = [-0.5;; 0.5]
	return dNdXi
end;

function getShapeAndDerivatives(type::ElementType, quad::QuadratureRules)
	if type == Quad1
		ShapeFunctions(getShapeFunctions(type, quad), getShapeDXi1(quad), getShapeDXi2(quad))
	elseif type == Lin1
		ShapeFunctions1D(getShapeFunctions(type, quad), getShapeDxi(quad))
	end
end;

function getJacobian(iQuad, shapeFuncs::ShapeFunctions, elePoints)
	dNdXi1 = shapeFuncs.dNdXi1[iQuad,:]'
	dNdXi2 = shapeFuncs.dNdXi2[iQuad,:]'
	jacobian = [
		dNdXi1 * elePoints[:,1]; dNdXi1 * elePoints[:,2];;
		dNdXi2 * elePoints[:,1]; dNdXi2 * elePoints[:,2]]

	detJ = (jacobian[1,1] * jacobian[2,2]) - (jacobian[2,1] * jacobian[1,2])
	dNdX = jacobian \ [dNdXi1; dNdXi2]
	return detJ, dNdX
end;