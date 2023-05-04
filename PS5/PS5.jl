using Plots
# using LateXStrings # Pkg.add("LaTeXStrings") # https://github.com/stevengj/LaTeXStrings.jl
using Interpolations # Pkg.add("Interpolations") # https://github.com/JuliaMath/Interpolations.jl
using ForwardDiff # Pkg.add("ForwardDiff") # https://github.com/JuliaDiff/ForwardDiff.jl
using Optim # Pkg.add("Optim") # https://julianlsolvers.github.io/Optim.jl/stable/
using Optim: converged, maximum, maximizer, minimizer, iterations
using Roots # Pkg.add("Roots") # https://github.com/JuliaMath/Roots.jl
using Parameters # Pkg.add("Parameters") # https://github.com/mauro3/Parameters.jl
using Distributions #Pkg.add("Distributions")
using QuadGK # Pkg.add("QuadGK") # https://juliamath.github.io/QuadGK.jl/latest/
using LinearAlgebra
using Random
using Statistics
# Call Scaled Interpolation Functions
include("/Users/juanholguin/Documents/GitHub/Computational-Macro/Ps3/ScaledInter.jl")
println(" ")
println("------------------------")
println("Integration in Julia")
println("PWD: ",pwd())
println("This code uses Plots, Interpolations, Dierckx, ForwardDiff, Optim, ")
println("   Roots, Parameters, ScaledInterpolation, Distributions, QuadGK, LinearAlgebra, Random")
println("Optimization in the context of the neoclassical growth model")
println("------------------------")
println(" ")

#-----------------------------------------------------------
#-----------------------------------------------------------
# Set random seed
Random.seed!(3112);

#-----------------------------------------------------------
# Define a markov process struct
    # Generate structure for markov processes using Parameters module
    @with_kw struct MP
        # Model Parameters
        N::Int64 # Number of states
        grid     # Grid of discrete markov process
        Π        # Transition matrix
        PDF      # Stationary distribution
        CDF      # Stationary distribution
    end


#-----------------------------------------------------------
#-----------------------------------------------------------
# Tauchen (1986)
    # Objective is to discretize AR(1) process: z'=ρz+η, η~N(0,σ)
    # Code from Kopecky & Suen (2010)
    # Inputs:
        # ρ - Process persisntence
        # σ - Innovation standard deviation
        # N - Size of the grid
        # Ω - Grid expansion in number of standard devaitions (Optional)
    # Outputs:
        # z - Grid of N equally spaced points covering [-Ωσ,Ωσ]
        # Π - Transition matrix, a stochastic matrix (sums to 1 across columns)
        # PDF_z, CDF_z - Stationary distribution of z
function Tauchen86(ρ,σ,N,Ω::Any=3)
    # Create z grid
        z = range(-Ω*σ/sqrt(1-ρ^2),Ω*σ/sqrt(1-ρ^2),length=N)
    # Define intermediate step length
        h = (z[2]-z[1])/2
    # Define auxiliary matrices
        z_0 = repeat(z ,1,N) # Matrix of today's z each row is a value, columns are equal
        z_1 = repeat(z',N,1) # Matrix of tomorrow's z each column is a value, rows are equal
    # Define intervals
        z_lim = zeros(N,N,2) # First matrix is lower bounds. Second matrix is uppor bounds.
        z_lim[:,1      ,1] .= -Inf
        z_lim[:,2:end  ,1] .=  ( z_1[:,2:end  ] - ρ*z_0[:,2:end  ] .- h )./σ
        z_lim[:,1:end-1,2] .=  ( z_1[:,1:end-1] - ρ*z_0[:,1:end-1] .+ h )./σ
        z_lim[:,end    ,2] .=  Inf
    # Define reference distribution
        # This line uses "Distributions"
        F(x) = cdf.(Normal(),x)
    # Fill in transition matrix
        Π_z = F.(z_lim[:,:,2]) - F.(z_lim[:,:,1])
        Π_z = Π_z./repeat(sum(Π_z,dims=2),1,N)
    # Get stationary distribution of markov chain
        PDF_z = real(eigvecs(Π_z')[:,end]); PDF_z = PDF_z/sum(PDF_z) ;
        CDF_z = cumsum(PDF_z)
    # Return
        return MP(N=N,grid=z,Π=Π_z,PDF=PDF_z,CDF=CDF_z)
end


#-----------------------------------------------------------
#-----------------------------------------------------------
# Rouwenhorst (1995)
    # Objective is to discretize AR(1) process: z'=ρz+η, η~N(0,σ)
    # Code from Kopecky & Suen (2010)
    # Inputs:
        # ρ - Process persisntence
        # σ - Innovation standard deviation
        # N - Size of the grid
    # Outputs:
        # z - Grid of N equally spaced points covering [-ψ,ψ]
        # Π - Transition matrix, a stochastic matrix (sums to 1 across columns)
        # PDF_z, CDF_z - Stationary distribution of z
function Rouwenhorst95(ρ,σ,N)
    # Define paramters for Rouwenhorst's approximation
        p = (1+ρ)/2
        q = p                   # Note: I am leaving q here for comparability with source
        ψ = σ*sqrt((N-1)/(1-ρ^2))
        s = (1-q)/(2-(p+q))     # Note: s=0.5, I leave it for comparability with source
    # Fill in transition matrix
    if N==2
        Π_z = [p 1-p ; 1-q q]
    else
        MP_aux = Rouwenhorst95(ρ,σ,N-1)
        o = zeros(N-1)
        Π_z = p*[MP_aux.Π o ; o' 0] + (1-p)*[o MP_aux.Π ; 0 o'] + (1-q)*[o' 0 ; MP_aux.Π o] + q*[0 o' ; o MP_aux.Π]
        # Adjust scale for double counting
        Π_z = Π_z./repeat(sum(Π_z,dims=2),1,N)
    end
    # Distribution
        PDF_z = pdf.(Binomial(N-1,1-s),(0:N-1))
        CDF_z = cumsum(PDF_z)
    # Create z grid
        z    = range(-ψ,ψ,length=N)
    # Return
        return MP(N=N,grid=z,Π=Π_z,PDF=PDF_z,CDF=CDF_z)
end


#-----------------------------------------------------------
#-----------------------------------------------------------
# Simulation of Markov processes

# Simulation function for discrete Markov process
# The result of the simulation is a Markov chain
    # Inputs:
        # Ns - Number of simulated periods
        # z  - Grid with levels of Markov process
        # Π  - Transition matrix of Markov process
        # N_0- Number of periods to throw before reporting simulation (Optional)
    # Output:
        # z_MC - Vector of Ns elemenrts with Markov Chain
function Simulation_MC(Ns,MP::MP,N_0=1000)
    # Compute conditional CDF
    Γ = cumsum(MP.Π,dims=2)
    # Allocate simulated vector
    z_ind    = zeros(Int64,N_0+Ns)
    z_MC     = zeros(N_0+Ns)
    # Starting value for simulation
    z_ind[1] = Int(ceil(length(MP.grid)/2))
    z_MC[1]  = MP.grid[z_ind[1]]
    # Simulate
    for i=2:Ns+N_0
        #= Option 1
        # Draw a uniform random number (r). Compare it with conditional CDF.
        # Conditional CDF given by cumsum of current row of Π, call it Γ
        # Γ[i,j] gives the probability z'<z_j
        # We are looking for j s.t Γ[i,j-1]<r<Γ[i,j]
        # Equivalently, the lowest j s.t. Γ[i,j]-r>0
            z_ind[i] = findmax(sign.(Γ[z_ind[i-1],:] .- rand()))[2]
        =#
        #= Option 2
        # Alternatively we can draw directly from conditional distributional
        # Distributional is categorical with P[z'=z_j] = Π[i,j]
        =#
        z_ind[i] = rand(Categorical(MP.Π[z_ind[i-1],:]))
        z_MC[i]  = MP.grid[z_ind[i]]
    end
    # Throw out first N_0 elements
    z_MC = z_MC[N_0+1:end]
    # Return result
    return z_MC
end

# Moments function for sample
function Moments_MC(z_MC)
    mean_MC      = mean(z_MC)
    std_MC       = std(z_MC)
    auto_corr_MC = cor(z_MC[1:end-1],z_MC[2:end])
    return mean_MC, std_MC, auto_corr_MC
end

ρ  = 0.95
σ  = 0.2
Ns = 10000
MP_T5  = Tauchen86(ρ,σ,5)
MP_T15 = Tauchen86(ρ,σ,15)

MP_R5  = Rouwenhorst95(ρ,σ,5)
MP_R15 = Rouwenhorst95(ρ,σ,15)


T_s_5  = Simulation_MC(Ns,MP_T5 ,0)
T_s_15 = Simulation_MC(Ns,MP_T15,0)

R_s_5  = Simulation_MC(Ns,MP_R5 ,0)
R_s_15 = Simulation_MC(Ns,MP_R15,0)

momt5=Moments_MC(T_s_5)
momt15=Moments_MC(T_s_15)
momr5=Moments_MC(R_s_5)
momr15=Moments_MC(R_s_15)

