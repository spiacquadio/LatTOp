function squareGrid(nNodes, width)
	nElements = (nNodes - 1)^2
	ncon = zeros(Int64, nElements, 4)
	pos = zeros(Float64, nNodes^2, 2)
	counter = 1

	posRange = LinRange(-width, width, nNodes)
	for i = 1:nNodes
		for j = 1:nNodes
			pos[counter,:] = [posRange[j], posRange[i]]
			counter += 1
		end
	end

	for i = 1:nNodes-1
		for j = 1:nNodes-1
			ncon[(i-1) * (nNodes-1) + j, 1] = (i - 1) * nNodes + j
			ncon[(i-1) * (nNodes-1) + j, 2] = (i - 1) * nNodes + j + 1
			ncon[(i-1) * (nNodes-1) + j, 3] = (i    ) * nNodes + j + 1
			ncon[(i-1) * (nNodes-1) + j, 4] = (i    ) * nNodes + j
		end
	end

	return ncon, pos
end

function getBoundaries(nNodes)
	boundaryNodes = zeros(Int64, 4*(nNodes-1))
	counter = 1

	maxNodes = nNodes^2
	finalRow = maxNodes - nNodes
	for i = 1:maxNodes
		if i <= nNodes
			boundaryNodes[counter] = i
			counter += 1
		elseif i > finalRow
			boundaryNodes[counter] = i
			counter += 1
		elseif mod(i, nNodes) == 0
			boundaryNodes[counter] = i
			counter += 1
		elseif mod(i-1, nNodes) == 0
			boundaryNodes[counter] = i
			counter += 1
		end
	end

	return boundaryNodes
end
