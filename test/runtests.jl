using RMPBT
using Test

@testset "RMPBT.jl" begin
    # --- This is the test function for the RMPBT

using RCall, DataStructures, LinearAlgebra, Distributions, Distributed, PDMats, StatsPlots, StatsBase, DataFrames, Plots, Kronecker, Statistics, Match, Dates, Random, ProximalOperators
#@rimport highmean as highmean ## Import the package highmean from r
#@rlibrary highmean
#---Checking the assumption of Matrix F Distributions

using LinearAlgebra
n1 = 50
n2 = 50
p = 200
α = 0.05
eye(p) = Matrix{Float64}(I,p, p)
#-- We consider a Block diagonal covariance matrix
    Σ = zeros(p,p)
    B = .85 .*eye(25) + 0.15 .*ones(25,25)
    num_b = floor(Integer, maximum(p)./25)
    for i in 1:num_b
     idx = Int.(collect((25 .*(i-1)+1):(i .*25)))
     Σ[idx,idx] = B
   end


#---
#-- Mean vector for group 1 (NUll group)
Random.seed!(2021)
μ1 = zeros(p)
X = zeros(p, n1)
rand!(MvNormal(μ1, Σ), X)
X = 1. .*X'


μ2 = zeros(p)
id = rand(range(1,p,step=1), 2) ## Selecy 10 entries to be non-zeros
μ2[id] = randn(length(id))
Y = zeros(p, n2)
rand!(MvNormal(μ2, Σ), Y)
Y = 1. .*Y'

#--- Compute the NULL distribution and the
N = 1000  # number of RPs
M = 1000  # number of samples

@time H00 = RMPBT.F_H0F(X, Y, N, M, α)

#--- compute all the other stuff!

@time Res1 = RMPBT.RMPBT_test(X , Y , 1. .*H00, 1, N, α) ## R2Mat projection matrix (sparse)
@time Res2 = RMPBT.RMPBT_test(X , Y , 1. .*H00, 2, N, α) ## Dense proejction matrix (iid normal(0,1))
@time Res3 = RMPBT.RMPBT_test(X , Y , 1. .*H00, 3, N, α) ## Proejction based on QR decomposition
end
