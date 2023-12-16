using Base
using Plots

# create array of normalised thermal compliance
cT_norm=Array{Float64}(undef, 100)
for i in 1:100
    cT_norm[i]=i/100
end
# create array of normalised structural compliance
cS_norm=Array{Float64}(undef, 100)
for i in 1:100
    cS_norm[i]=1-i/100
end
plot(size=(650,400))
plot!(cT_norm, cS_norm, xlabel="cT/cTmax", ylabel="cS/cSmax", label="cT/cS", legend=:bottomright, title="cT/cS vs cT/cTmax", xlim=(0,1), ylim=(0,1))
