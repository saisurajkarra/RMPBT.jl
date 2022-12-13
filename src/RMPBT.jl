module RMPBT

#-- This module contains the main functions
#--
export RMPBT
export R2Mat, RupR1, DensR, F_H0F, BayeF_TestR2Mat, BayeF_TestRupR1 , PvalFunc , BayeF_TestDenseMat, Tauγ
using RCall, DataStructures, LinearAlgebra, Distributions, Distributed, PDMats, StatsPlots, StatsBase, DataFrames, Plots, Kronecker, Statistics, Match, Dates, Random, ProximalOperators
#@rimport highmean #as highmean ## Import the package highmean from r
#@rlibrary highmean
#---Checking the assumption of Matrix F Distributions
#using LinearAlgebra, Distributions, Distributed, PDMats, StatsPlots, StatsBase
#- specify constraints

eye(p) = Matrix{Float64}(I,p, p)
#--covariance function
#function ReducRank(Y,X,r,Gam)
function toeplitz_filter!(A,x,σ)
    n = length(x)
    for i in 1:n
        A[i,i:end] = σ.* x[1:n - i + 1]
        A[i:end,i] = σ .*x[1:n - i + 1]
    end
end

# R2 matrix from Sri et al.
function ArmaMaT!(Σ, σ, γ, ρ)
        p = size(Σ,1)
        x0 = collect([1.0; 1.0;ρ.^(1:(p-2))].*[1.0;repeat(collect(γ); inner = p-1)])
        toeplitz_filter!(Σ, x0, σ)
end

function Tauγ(n1::Int64, n2::Int64, α::Float64)
  n = n1+n2
  mval = collect(2:(n-3))
  val = [cquantile(FDist(mval[i],n-mval[i]-1),α) for i in 1:length(mval)]
  m = mval[val .== minimum(val)][1]
  n0 = 1/(1/n1 + 1/n2)
  Fal = minimum(val)
  C0 = m*Fal/(m*Fal +n -m -1)
  tp0 = (C0 - m/(n-1))/((1-C0)*m/(n-1))
  τ0 = n0/tp0
   γ = exp(-.5*m*log1p(tp0) - .5*(n-1)*log1p(-tp0*C0/(1 + tp0)))
 return [τ0 γ m]
end

function R2Mat(p::Int64, m::Int64)
    val = zeros(p)
    randn!(val)
    Rfin = zeros(p,m)
    m0 = div(p,m) ## number of elements in each columns
      m1 = p%m ## remainder after we devide by m
          id1 = repeat(collect(1:m);inner=m0) ## inner = each and outer=times
      normval = zeros(m)
      if(m1 != 0)
       id1 = collect([id1; sample(collect(1:m), m1)])
      end
      for i in 1:p
          Rfin[i,id1[i]] = val[i]
      end
       Rfin = Rfin[randperm(p),:]

      for i in 1:m
         normval[i] = norm(Rfin[:,i])
      end
      #scale!(Rfin, 1 ./normval)
      Rfin*diagm(1 ./normval)
   #return Rfin
end

function RupR1(p::Int64,m::Int64)
    return  1. .*(qr(randn(p,m)).Q)*diagm(sign.(diag(qr(randn(p,m)).Q)))[1:p,1:m]
end

function DenseR(p::Int64,m::Int64)
    R = zeros(p,m)
    randn!(R)
    return R
end

function F_func1(X,Y, Im)
        n1, p = size(X)
        n2 = size(Y,1)
        mu1 = zeros(p)
        mu2 = zeros(p)
        mdiff = zeros(p)
        msum = zeros(p)
        eta = 1.0
        nu0 = 2*p
        covx = zeros(p,p)
        covy = zeros(p,p)
        n = n1 + n2 - 2
        mean!(mu1,X') ## column means
        mean!(mu2,Y') ## column means
        ndel = (n1+n2)/(n1*n2)
        n0 = (n1*n2)/(n1+n2)

        for i=1:p
           mdiff[i] = mu2[i] - mu1[i]
        end

        mul!(covx, cov(X), Im .*(n1-1))
        mul!(covy, cov(Y), Im .*(n2-1))

        spool = zeros(p,p)

        for i=1:length(spool)
          spool[i] = (covx[i] + covy[i])./n
        end

        return ((n - p + 1)/(n*p*ndel)).*dot(mdiff, inv(spool), mdiff) # improper prior
end

function F_H0F(Xd, Yd, N, M, α)
    n1, p = size(Xd)
    n2 = size(Yd,1)
    τ0, γ0, m = Tauγ(n1 , n2 , α)
     m = Int(m)
    n = n1 + n2
    n_del0 = n1*n2/(n1+n2)
    n_del1 = 1/(1/n_del0 + 1/τ0)
    t0 = n_del0/τ0
    c1 = m/((1 + t0)*(n-m-1)) ## f[1]
    c3 = m/(n-m-1)
    c4 = .5*(n - 1)
    c5 = .5*m*log(1 + t0)

    Σ0 = eye(p)
    nu0 = 2.0*p
    R0 = [RupR1(p,m) for _ in 1:N]
    Im = Matrix{Float64}(I,m,m)

    X = zeros(n1, p) #Array{Float64}[]
    Y = zeros(n2, p) #Array{Float64}[]

    f = zeros(1)
    outF1 = zeros(N)
    outF2 = zeros(N)
    res = zeros(8, M)
    Fdi = FDist(m, n-m-1)
    println("<---- Start the loops --->")
    println("Start time: ",now())
    for i in 1:M
       randn!(Y)
       randn!(X)
       g1 = Array{Float64}[]
       g2 = Array{Float64}[]
       [push!(g1, X*R0[j]) for j in 1:N]
       [push!(g2, Y*R0[j]) for j in 1:N]
         for j in 1:N
            f = F_func1(g1[j], g2[j], Im)
            outF1[j] = -c5 - c4*(log(1 + c1*f[1]) - log(1 + c3*f[1])) ## based on f[1]
            outF2[j] = ccdf(Fdi,f[1])
         end
       res[1,i] = n1
       res[2,i] = n2
       res[3,i] = p
       res[4,i] = γ0
       res[5,i] = τ0
       res[6,i] = mean(outF1)
       res[7,i] = mean(outF1 .> log(γ0))
       res[8,i] = mean(outF2)

      i%1000 == 0 ? println("iter $(i)") : nothing
    end
    println("End time: ",now())
    println("<---- Done! --->")
   return res'
end

function BayeF_TestR2Mat(X, Y, N, α)
    Random.seed!(124567)
    n1 = size(X,1)
    n2 = size(Y,1)
    p =  size(Y,2)
    pam = Tauγ(n1, n2, α)
    m = floor(Integer, pam[3])
    τ0 = pam[1]
    γ0 = pam[2]
    Im = Matrix{Float64}(I, m, m)
    M = 1
    n = n1 + n2
    n_del0 = n1*n2/(n1+n2)
    n_del1 = 1/(1/n_del0 + 1/τ0)
    t0 = n_del0/τ0
    c1 = m/((1 + t0)*(n-m-1)) ## f[1]
    c3 = m/(n-m-1)
    c4 = .5*(n - 1)
    c5 = .5*m*log(1 + t0)
    nu0 = 2.0*p
    R0 = [R2Mat(p,m) for i in 1:N]
    f = zeros(1)
    outF1 = zeros(N)
    outF2 = zeros(N)
    res = zeros(11, M)
    Fdi = FDist(m, n-m-1)
    for i in 1:1
       g1 = Array{Float64}[]
       g2 = Array{Float64}[]
       [push!(g1, X*R0[j]) for j in 1:N]
       [push!(g2, Y*R0[j]) for j in 1:N]
         for j in 1:N
            f = F_func1(g1[j],g2[j], Im)
            outF1[j] = -c5 - c4*(log(1 + c1*f[1]) - log(1 + c3*f[1])) ## based on f[1]
            outF2[j] = ccdf(Fdi, f[1])
         end
       res[1,i] = n1
       res[2,i] = n2
       res[3,i] = p
       res[4,i] = γ0
       res[5,i] = τ0
       res[6,i] = mean(outF1)
       res[7,i] = mean(outF1 .> log(γ0))
       res[8,i] = mean(outF2)
        res[9,i] = rcopy(getfield(highmean, Symbol("apval_Bai1996"))(Y, X))[:pval]
      # res[10,i] = rcopy(getfield(highmean, Symbol("apval_Sri2008"))(Y, X))[:pval]
       res[11,i] = rcopy(getfield(highmean, Symbol("apval_Chen2010"))(Y, X))[:pval]
    end
   return res'
end

function BayeF_TestRupR1(X, Y, N, α)
  Random.seed!(124567)
    n1 = size(X,1)
    n2 = size(Y,1)
    p =  size(Y,2)
    pam = Tauγ(n1, n2, α)
    m = Int(pam[3])
    τ0 = pam[1]
    γ0 = pam[2]
    Im = Matrix{Float64}(I, m, m)
    M =1
    n = n1 + n2
    n_del0 = n1*n2/(n1+n2)
    n_del1 = 1/(1/n_del0 + 1/τ0)
    t0 = n_del0/τ0
    c1 = m/((1 + t0)*(n-m-1)) ## f[1]
    c3 = m/(n-m-1)
    c4 = .5*(n - 1)
    c5 = .5*m*log(1 + t0)
    nu0 = 2.0*p
    R0 = [RupR1(p,m) for i in 1:N]
    f = zeros(1)
    outF1 = zeros(N)
    outF2 = zeros(N)
    res = zeros(11, M)
    Fdi = FDist(m, n-m-1)
    for i in 1:1
       g1 = Array{Float64}[]
       g2 = Array{Float64}[]
       [push!(g1, X*R0[j]) for j in 1:N]
       [push!(g2, Y*R0[j]) for j in 1:N]
         for j in 1:N
            f = F_func1(g1[j],g2[j], Im)
            outF1[j] = -c5 - c4*(log(1 + c1*f[1]) - log(1 + c3*f[1])) ## based on f[1]
            outF2[j] = ccdf(Fdi,f[1])
         end
       res[1,i] = n1
       res[2,i] = n2
       res[3,i] = p
       res[4,i] = γ0
       res[5,i] = τ0
       res[6,i] = mean(outF1)
       res[7,i] = mean(outF1 .> log(γ0))
       res[8,i] = mean(outF2)
        res[9,i] = rcopy(getfield(highmean, Symbol("apval_Bai1996"))(Y, X))[:pval]
      # res[10,i] = rcopy(getfield(highmean, Symbol("apval_Sri2008"))(Y, X))[:pval]
       res[11,i] = rcopy(getfield(highmean, Symbol("apval_Chen2010"))(Y, X))[:pval]
    end
   return res'
end

function BayeF_TestDenseMat(X, Y, N, α)
  Random.seed!(124567)
    n1 = size(X,1)
    n2 = size(Y,1)
    p =  size(Y,2)
    pam = Tauγ(n1, n2, α)
    m = floor(Integer, pam[3])
    τ0 = pam[1]
    γ0 = pam[2]
    Im = Matrix{Float64}(I, m, m)
    M =1
    n = n1 + n2
    n_del0 = n1*n2/(n1+n2)
    n_del1 = 1/(1/n_del0 + 1/τ0)
    t0 = n_del0/τ0
    c1 = m/((1 + t0)*(n-m-1)) ## f[1]
    c3 = m/(n-m-1)
    c4 = .5*(n - 1)
    c5 = .5*m*log(1 + t0)
    nu0 = 2.0*p
    R0 = [DenseR(p,m) for i in 1:N]
    f = zeros(1)
    outF1 = zeros(N)
    outF2 = zeros(N)
    res = zeros(11, M)
    Fdi = FDist(m, n-m-1)
    for i in 1:1
       g1 = Array{Float64}[]
       g2 = Array{Float64}[]
       [push!(g1, X*R0[j]) for j in 1:N]
       [push!(g2, Y*R0[j]) for j in 1:N]

         for j in 1:N
            f = F_func1(g1[j],g2[j], Im)
            outF1[j] = -c5 - c4*(log(1 + c1*f[1]) - log(1 + c3*f[1])) ## based on f[1]
            outF2[j] = ccdf(Fdi,f[1])
         end
       res[1,i] = n1
       res[2,i] = n2
       res[3,i] = p
       res[4,i] = γ0
       res[5,i] = τ0
       res[6,i] = mean(outF1)
       res[7,i] = mean(outF1 .> log(γ0))
       res[8,i] = mean(outF2)
        res[9,i] = rcopy(getfield(highmean, Symbol("apval_Bai1996"))(Y, X))[:pval]
      # res[10,i] = rcopy(getfield(highmean, Symbol("apval_Sri2008"))(Y, X))[:pval]
       res[11,i] = rcopy(getfield(highmean, Symbol("apval_Chen2010"))(Y, X))[:pval]
    end
    return res'
end

function RMPBT_test(X::Array, Y::Array, H00::Array, met::Int64, N::Int64, α::Float64)
    n1 = size(X,1)
    n2 = size(Y,1)
    p = size(Y,2)
    Res2 = zeros(9)
    pam00 = Tauγ(n1,n2,α)
    fun =  @match met begin
       1 => BayeF_TestR2Mat
       2 => BayeF_TestDenseMat
       3 => BayeF_TestRupR1
    end
    ot = fun(X, Y, N, α)
    Res2[1] = mean(ot[7] .< H00[:,7]) # pvalue
    Res2[2] = ot[7]
    Res2[3] = quantile(H00[:,7], .95)
    Res2[4] = mean(H00[:,8] .< ot[8]) #pvalue
    Res2[5] = ot[8]
    Res2[6] = quantile(H00[:,8], .05)
    Res2[7:end] = ot[9:end]
    println("done")
    out = DataStructures.SortedDict("τ0"=>pam00[1],"γ" => pam00[2],"m" => pam00[3],"RMPBT_pval"=>Res2[1],"RMPBT_phi"=>Res2[2],"RMPBT_phi0"=>Res2[3],
      "RAPTT_pval" => Res2[4],"RAPTT_phi" => Res2[5],"RAPTT_u0" => Res2[6],"Bai96" => Res2[7],"Chen10" => Res2[9])
    return out
end

println("**** Module all loaded ****")

end


