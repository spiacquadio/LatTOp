using Zygote
using ChainRulesCore
using ReverseDiff

function getf2(x)
	f2 = zeros(2);
	for i = 1:2
		f2[i] = x;
	end
	return f2;
end;
Zygote.@adjoint getf2(x) = getf2(x), c̄ -> (sum(c̄),)

function funcMutate(x)
	f1 = [x^2, x];
	f2 = getf2(x)
	f3 = f1 .+ f2
	return f3;
end;

function centralDifference(funcPtr, x)
	esp = 1e-6;
	l = funcPtr(x - esp);
	r = funcPtr(x + esp);
	return (r - l) / (2 * esp);
end;
g = (x) -> Zygote.jacobian(funcMutate, x);