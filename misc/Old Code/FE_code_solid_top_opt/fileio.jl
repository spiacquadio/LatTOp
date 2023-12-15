using Printf

## ----- Export to VTK function ----- ##
function export2vtk(filepath, points, ncon, sol = 0.0, cellData = 0.0)
	file = open(filepath, "w");
	write(file, "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n");
	write(file, "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\" compressor=\"vtkZLibDataCompressor\">\n");
	write(file, "\t<UnstructuredGrid>\n");
	
	nNodes = size(points)[1];
	nElements = size(ncon)[1];
	write(file, @sprintf("\t\t<Piece NumberOfPoints=\"%i\" NumberOfCells=\"%i\">\n", nNodes, nElements));

	# Write node position data
	write(file, "\t\t\t<Points>\n");
	write(file, "\t\t\t\t<DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n");
	for i = 1:nNodes
		if size(points)[2] == 2
			write(file, @sprintf("\t\t\t\t\t%f %f %f\n", points[i,1], points[i,2], 0.0));
		elseif size(points)[2] == 3
			write(file, @sprintf("\t\t\t\t\t%f %f %f\n", points[i,1], points[i,2], points[i,3]));
		end
	end
	write(file, "\t\t\t\t</DataArray>\n");
	write(file, "\t\t\t</Points>\n");


	# Write connectivity data
	write(file, "\t\t\t<Cells>\n");
	write(file, "\t\t\t\t<DataArray type=\"Int64\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">\n");
	for i = 1:nElements
		if length(ncon[i]) == 2
			write(file, @sprintf("\t\t\t\t\t%i %i\n", ncon[i][1]-1, ncon[i][2]-1));
		else length(ncon[i]) == 4
			write(file, @sprintf("\t\t\t\t\t%i %i %i %i\n", ncon[i][1]-1, ncon[i][2]-1, ncon[i][3]-1, ncon[i][4]-1));
		end
	end
	write(file, "\t\t\t\t</DataArray>\n");

	# Write element offset information
	write(file, "\t\t\t\t<DataArray type=\"Int64\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">\n");
	counter = 0
	for i = 1:nElements
		if length(ncon[i]) == 2
			write(file, @sprintf("\t\t\t\t\t%i\n", counter+2));
			counter += 2
		elseif length(ncon[i]) == 4
			write(file, @sprintf("\t\t\t\t\t%i\n", counter+4));
			counter += 4
		end
	end
	write(file, "\t\t\t\t</DataArray>\n");

	# Write element type information
	write(file, "\t\t\t\t<DataArray type=\"UInt8\" Name=\"types\" NumberOfComponents=\"1\" format=\"ascii\">\n");
	for i = 1:nElements
		if length(ncon[i]) == 2
			write(file, @sprintf("\t\t\t\t\t%i\n", 3));
		elseif length(ncon[i]) == 4
			write(file, @sprintf("\t\t\t\t\t%i\n", 9));
		end
	end
	write(file, "\t\t\t\t</DataArray>\n");
	write(file, "\t\t\t</Cells>\n");

	# Write solution data
	if length(sol) > 1
		nDofs = Int(length(sol) / nNodes);
		write(file, "\t\t\t<PointData>\n");
		write(file, @sprintf("\t\t\t\t<DataArray type=\"Float64\" Name=\"u\" NumberOfComponents=\"%i\" format=\"ascii\">\n", nDofs));
		write(file, @sprintf("\t\t\t\t\t"));
		for i = 1:nNodes
			for j = 1:nDofs
				write(file, @sprintf("%f ", sol[nDofs * i - (nDofs - j)]))
			end
		end
		write(file, "\n")
		write(file, "\t\t\t\t</DataArray>\n");
		write(file, "\t\t\t</PointData>\n");
	end
	if length(cellData) > 1
		write(file, "\t\t\t<CellData>\n");
		write(file, "\t\t\t\t<DataArray type=\"Float64\" Name=\"Density\" NumberOfComponents=\"1\" format=\"ascii\">\n");
		counter = 1
		for i = 1:nElements
			if length(ncon[i]) == 2
				write(file, @sprintf("\t\t\t\t\t%f\n", 0.0));
			elseif length(ncon[i]) == 4
				write(file, @sprintf("\t\t\t\t\t%f\n", cellData[counter]));
				counter += 1
			end
		end
		write(file, "\t\t\t\t</DataArray>\n");
		write(file, "\t\t\t</CellData>\n");
	end

	write(file, "\t\t</Piece>\n");
	write(file, "\t</UnstructuredGrid>\n");
	write(file, "</VTKFile>\n");
	close(file);
end;

## ----- Mesh import functions ----- ## 
function importMsh(filepath)
	file = open(filepath, "r")

	# Preallocate arrays and maps
	ncon = Vector{Vector{Int64}}([])
	pos = []
	groups = Dict()
	entityMap = Vector{Vector{Int}}([])
	elementTags = Vector{Vector{Int}}([])

	# Read through file
	while true
		line = readline(file)

		# Read mesh format
		if contains(line, "MeshFormat")
			version, type, _ = parseMeshFormat(file)
			if version != 4.1
				@printf "[Error] Version not supported! : %.1f" version
				break
			end
			if type == true
				@printf "[Error] File type not supported! : %i" type
				break
			end

		# Get physical groups
		elseif contains(line, "PhysicalNames")
			groups = parsePhysicalNames(file)

		# Get mesh entities
		elseif contains(line, "Entities")
			entityMap = parseEntities(file)

		# Read node positions
		elseif contains(line, "Nodes")
			pos = parseNodes(file)

		# Read element connectivity
		elseif contains(line, "Elements")
			ncon, elementTags = parseElements(file)
		
		# If EOF reached
		elseif isempty(line)
			break
			
		# Ignore other blocks
		else
			@printf "[Error] Block not recognised! : %s" line
		end
	end

	# Process the groups to get element tags in one array
	boundaryElements = Dict()
	for key in keys(groups)
		group, dim = groups[key]
		id = []

		# Loop through entity map
		for i = 1:size(entityMap)[1]
			if dim == entityMap[i][1]
				# Loop through groups entity is assigned to
				for j = 3:size(entityMap[i])[1]
					if group == entityMap[i][j]
						push!(id, entityMap[i][2])
					end
				end
			end
		end
		
		# Loop through elementTags to find
		idx = []
		for i = 1:size(elementTags)[1]
			if dim == elementTags[i][1]
				for j = 1:size(id)[1]
					if id[j] == elementTags[i][2]
						push!(idx, i)
					end
				end
			end
		end

		# Concat all elements to one groups
		elements = []
		for i = 1:size(idx)[1]
			elements = [elements; elementTags[idx[i]][3:end]]
		end
		boundaryElements[key] = elements
	end

	boundaryNodes = Dict()
	for key in keys(boundaryElements)
		nElements = length(ncon[boundaryElements[key]])
		nconTmp = ncon[boundaryElements[key]]
		nodes = Vector{Int64}(undef, 2*nElements)
		for i = 1:nElements
			nodes[2*i-1] = nconTmp[i][1]
			nodes[2*i] = nconTmp[i][2]
		end
		boundaryNodes[key] = unique(nodes)
	end

	return ncon, pos, boundaryElements, boundaryNodes
end

function parseMeshFormat(file)
	# Preallocate return values
	version = 0.0
	type = false
	dataSize = 0

	# Parse block
	while true
		line = readline(file)
		if contains(line, '$')
			break;
		end

		parts = split(line)
		version = parse(Float64, parts[1])
		type = parse(Bool, parts[2])
		dataSize = parse(Int64, parts[3])
	end

	return version, type, dataSize
end;

function parsePhysicalNames(file)
	# Preallocate map
	map = Dict()

	# Parse block
	line = readline(file)
	n = parse(Int64, line)
	for i = 1:n
		parts = split(readline(file))
		dim = parse(Int64, parts[1])
		val = parse(Int64, parts[2])
		name = parts[3][2:end-1]
		map[name] = [val dim]
	end
	if !contains(readline(file), '$')
		@printf "[Error] End of block not reached!"
	end

	return map
end;

function parseEntities(file)
	parts = split(readline(file))
	nPoints = parse(Int64, parts[1])
	nCurves = parse(Int64, parts[2])
	nSurfaces = parse(Int64, parts[3])
	nVolumes = parse(Int64, parts[4])

	map = Vector{Vector{Int}}([])

	# Parse points
	for i = 1:nPoints
		parts = split(readline(file))
		tag = parse(Int64, parts[1])
		nPhysicalTags = parse(Int64, parts[4])
		if nPhysicalTags > 0
			tmpArr = zeros(nPhysicalTags + 2)
			tmpArr[1] = 0
			tmpArr[2] = tag
			for j = 1:nPhysicalTags
				tmpArr[j+2] = parse(Int64, parts[4+j])
			end
			push!(map, tmpArr)
		end
	end

	# Parse Curves
	for i = 1:nCurves
		parts = split(readline(file))
		tag = parse(Int64, parts[1])
		nPhysicalTags = parse(Int64, parts[8])
		if nPhysicalTags > 0
			tmpArr = zeros(nPhysicalTags + 2)
			tmpArr[1] = 1
			tmpArr[2] = tag
			for j = 1:nPhysicalTags
				tmpArr[j+2] = parse(Int64, parts[8+j])
			end
			push!(map, tmpArr)
		end
	end

	# Parse Surfaces
	for i = 1:nSurfaces
		parts = split(readline(file))
		tag = parse(Int64, parts[1])
		nPhysicalTags = parse(Int64, parts[8])
		if nPhysicalTags > 0
			tmpArr = zeros(nPhysicalTags + 2)
			tmpArr[1] = 2
			tmpArr[2] = tag
			for j = 1:nPhysicalTags
				tmpArr[j+2] = parse(Int64, parts[8+j])
			end
			push!(map, tmpArr)
		end
	end

	# Parse volumes
	for i = 1:nVolumes
		@printf "[Warning] Volumes not implimented!"
	end

	if !contains(readline(file), '$')
		@printf "[Error] End of block not reached!"
	end

	return map;
end;

function parseNodes(file)
	parts = split(readline(file))
	nEntities = parse(Int64, parts[1])
	nNodes = parse(Int64, parts[2])

	nodes = Matrix{Float64}(undef, nNodes, 3)
	for e = 1:nEntities
		parts = split(readline(file))
		nBlock = parse(Int64, parts[4])
		tmpArr = Vector{Int64}(undef, nBlock)
		for i = 1:nBlock
			tmpArr[i] = parse(Int64, readline(file))
		end
		for i = 1:nBlock
			pos = split(readline(file))
			for j = 1:3
				nodes[tmpArr[i],j] = parse(Float64, pos[j])
			end
		end
	end

	if !contains(readline(file), '$')
		@printf "[Error] End of block not reached!"
	end

	return nodes;
end;

function parseElements(file)
	parts = split(readline(file))
	nEntities = parse(Int64, parts[1])
	nElements = parse(Int64, parts[2])

	ncon = Vector{Vector{Int64}}(undef, nElements)
	map = Vector{Vector{Int64}}(undef, nEntities)
	for e = 1:nEntities
		parts = split(readline(file))
		dim = parse(Int64, parts[1])
		tag = parse(Int64, parts[2])
		type = parse(Int64, parts[3])
		nBlock = parse(Int64, parts[4])

		entitiesTemp = Vector{Int64}(undef, nBlock + 2)
		entitiesTemp[1] = dim;
		entitiesTemp[2] = tag;

		nNodePerElement = 0
		if type == 1
			nNodePerElement = 2
		elseif type == 3
			nNodePerElement = 4
		else
			@printf "[Error] Element not supported!"
		end

		nconTemp = zeros(nNodePerElement)
		for i = 1:nBlock
			parts = split(readline(file))
			idx = parse(Int64, parts[1])
			entitiesTemp[i+2] = idx
			for j = 1:nNodePerElement
				nconTemp[j] = parse(Int64, parts[j+1])
			end
			ncon[idx] = nconTemp;
		end

		map[e] = entitiesTemp
	end

	if !contains(readline(file), '$')
		@printf "[Error] End of block not reached!"
	end

	return ncon, map
end;
