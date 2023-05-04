using Parameters
@with_kw mutable struct Par
z::Float64 = 1    ; # Productivity
α::Float64 = 1/3  ; # Production function
β::Float64 = 0.98 ; # Discount factor
η::Float64=1
σ::Float64=2
δ::Float64=1
χ::Float64=5.0
# VFI Paramters
max_iter::Int64   = 2000  ; # Maximum number of iterations
dist_tol::Float64 = 1E-9  ; # Tolerance for distance
# Howard's Policy Iterations
H_tol::Float64    = 1E-9  ; # Tolerance for policy function iteration
end
p=Par()
####### Utility Function ####

function chi(t,p)
    @unpack z, α, β, η, σ, δ = p
    l_ss1 = (((((z*(z*β*α)^(α/(1-α))-δ*(z*β*α)^(1/(1-α))))^(-σ))*z*(1-α)*(z*β*α)^(α/(1-α)))/t)^(1/(η+σ))-(0.4)
    return l_ss1
end


chi1=(t->chi(t,p))

using Plots
plot(chi1,0:100)

using Roots

p = reconstruct(p,χ=find_zero(chi1,0.1))




function utility(k,kp,l)
    @unpack z, α, β, η, σ,δ,χ = p
    c=(z*k.^(α)*l.^(1-α))+((k*(1-δ))-kp)
    if c>0
        return ((c.^(1-σ))/(1-σ)) - χ*((l.^(1+η)))/(1+η)
    else
        return -Inf
    end
end

###Steady State
function ss(p)
    @unpack z, α, β, η, σ,δ,χ = p
    l_ss = ((((z*(z*β*α)^(α/(1-α))-δ*(z*β*α)^(1/(1-α))))^(-σ)*z*(1-α)*(z*β*α)^(α/(1-α)))/χ)^(1/(η+σ));
    k_ss=l_ss*((z*β*α)^(1/(1-α)));
    y_ss=z*k_ss^α*l_ss^(1-α);
    c_ss=y_ss-k_ss;
    w_ss=z*(1-α)k_ss^α*l_ss^(-α);
    r_ss=α*k_ss^(α-1)*l_ss^(1-α);
    return k_ss,y_ss,c_ss,r_ss,w_ss,l_ss
end
ss(p)


#### Grid Function ###

function kgrid(lengthk,p)
    #I want the SS to define the max value of the Grid
    k_ss,y_ss,c_ss,r_ss,w_ss,l_ss=ss(p)
    grid=range(1E-3,2*k_ss,length=lengthk)
    return grid
end



##### Labor FOC ####
function labor(l,k,kp)
    @unpack z, α, β, η, σ,δ,χ = p
    t=l-(((((z*k^(α)*l^(1-α)-kp + k*(1-δ))^(-σ))*z*(1-α)*l^(-α)*k^(α))/χ)^(1/η))
    return t
end


nk=50

##### Optimal Labor for each k and k'
function labor_dec(length,p)

    @unpack z, α, β, η, σ,δ,χ  = p
    gridk=kgrid(length,p)
    labord=zeros(nk,nk)
    for i=1:nk
        for j=1:nk
            labor1=(l->labor(l,gridk[i],gridk[j]))
            labord[i,j]=try
            find_zero(labor1, 0.1)
            catch
                -Inf
            end   
              
        end
    end
    return labord
end

labord=labor_dec(50,p)
#### Checking Labor SS

kss=ss(p)

labor11=(l->labor(l,kss[1],kss[1]))
find_zero(labor11,0.5)  #### Correct Steady State of Labor if I use KSS

######

for i in 1:nk
    for j in 1:nk
        if labord[i, j] > 1
            labord[i, j] = 1
        end
    end
end

for i in 1:nk
    for j in 1:nk
        if labord[i, j] == -Inf
            labord[i, j] = 0
        end
    end
end



####### Bellman Operator ####

function Bellman(V_old,gridk)

    @unpack z, α, β, η, σ,δ,χ = p
    n_k = length(gridk)
    V = zeros(n_k)
    G_kp = fill(0,n_k)
    G_c = zeros(n_k)
    G_n = zeros(n_k)
    Gp = zeros(n_k)
    for i = 1:n_k
        V_aux = zeros(n_k)
        for j = 1:n_k
            V_aux[j]=try
             utility(gridk[i],gridk[j],labord[i,j]) + β*V_old[j]
            catch
                -Inf
            end
        end
        V[i], G_kp[i] = findmax(V_aux)
        G_n[i]=labord[i,G_kp[i]]
        G_c[i] = z * gridk[i]^(α)*G_n[i]^(1-α) - gridk[G_kp[i]]+gridk[i]*(1-δ)
        
        Gp[i]=gridk[G_kp[i]]
    end
    return V, G_kp, G_c, G_n,Gp
end



function VFI(T::Function,gridk)
    @unpack z, α, β, η, σ, δ,max_iter,dist_tol = p
    n_k    = length(gridk) ; # Number of grid nodes
    V_old  = zeros(n_k)     ; # Initial value, a vector of zeros
    V_dist = 1              ; # Initialize distance
    for iter=1:max_iter
        # Update value function
        V_new, G_kp, G_c,G_n,Gp = T(V_old)
        # Update distance and iterations
        V_dist = maximum(abs.(V_new./V_old.-1))
        # Update old function
        V_old  = V_new
        # Report progress
        if mod(iter,100)==0
            println("   VFI Loop: iter=$iter, dist=",100*V_dist,"%")
        end
        # Check convergence and return results
        if V_dist<=dist_tol
            println("VFI - Grid Search - n_k=$n_k")
            println("Iterations = $iter and Distance = ",100*V_dist,"%")
            println("------------------------")
            println(" ")
            return V_new, G_kp, G_c,G_n,Gp 
        end
    end
    # If loop ends there was no convergence -> Error!
    error("Error in VFI - Grid Search - Solution not found")
end

function Solve_VFI_loop(n_k)
    # Get Grid
    gridk=kgrid(nk,p)
    # Solve VFI
    V, G_kp, G_c,G_n,Gp = VFI(x->Bellman(x,gridk),gridk)
    # Return Solution
    return V, G_kp, G_c,G_n,Gp
end

Vd,Gd,Cd,gn,gkp=Solve_VFI_loop(200)
