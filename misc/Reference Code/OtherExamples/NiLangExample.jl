using NiLang

@i function f(res, y, x)
	for i=1:size(x)[1]
		y += x[i] ^ 2
	end
	res += sqrt(y)
end

res_out, y_out, x_out = f(0.0, 0.0, [1, 2, 3.0])
(~f)(res_out, y_out, x_out)
∂res, ∂y, ∂x = NiLang.AD.gradient(Val(1), f, (0.0, 0.0, [1, 2, 3.0]))