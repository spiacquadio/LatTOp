using ForwardDiff

#=
	In order to ForwardDiff, the function needs to have only one input variable, and one output.
	This is an issue as most functions have to return more than one array or matrix.
=#

function doSomething(x)
	nVals = size(x)[1]
	f = Matrix{Float64}(undef, nVals, nVals);
	fill!(f, 0.0);

	for i = 1:nVals
		for j = 1:nVals
			f[i,j] = i * j
		end
	end
	f = f .* x
	return f;
end





