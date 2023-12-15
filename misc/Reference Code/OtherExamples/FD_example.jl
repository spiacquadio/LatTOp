using FiniteDifferences

f(x) = (2.0 * x[1]^2) + (4.0 * x[2]^2)
x = [2.0, 2.0]

grad(central_fdm(2, 1), f, x)