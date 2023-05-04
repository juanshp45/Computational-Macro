using Plots

function utility(c)
    u=log(c)
    return u
end

nc=10


gridc=Vector(range(0.05,2,length=nc))

##### Interpolation #####

gridu=map(utility,gridc) 

function LinearInt(x,grid::Array,y::Array)
    ind = findmax(sign.(grid .- x))[2] - 1 
    A_x = (grid[ind+1] - x)/(grid[ind+1]-grid[ind])
    y_hat=A_x*y[ind]+(1-A_x)*y[ind+1]
    return y_hat
end

x1=range(0.05,2,100)
x2=range(0.05,2,1000)
x3=range(0.05,2,10000)
y_hat1 = [LinearInt(xi, gridc, gridu) for xi in x1]
y_hat2 = [LinearInt(xi, gridc, gridu) for xi in x2]
y_hat3 = [LinearInt(xi, gridc, gridu) for xi in x3]

gr()
plot(x1,y_hat1)
plot!(x1,map(utility,x1))

gr()
plot(x2,y_hat2)
plot!(x2,map(utility,x2))

gr()
plot(x3,y_hat3,label="Interpolation")
plot!(x3,map(utility,x3),title="Grid=10",label="Function Value")
savefig("grid10.png")



#### Assesing Accuracy ####
using Statistics

function MSE1(int::Array,truev::Array)
    MSE=mean((truev - int).^2)
    return MSE
end
###Creating vector for different grids
gridnumber=[5,10,50,100,1000]

asses=zeros(5)
for i=1:5
    
    gridc=Vector(range(0.05,2,length=gridnumber[i]))
    gridu=map(utility,gridc)
    x=range(0.05,2,1000) ### Grid to evaluate Interpolation
    y_hat = [LinearInt(xi, gridc, gridu) for xi in x]
    asses[i]=MSE1(y_hat,map(utility,x))
end


########## Cubic Splines #####
### All the functions here were taken from Sergio's Code #######

# Cubic spline interpolation
#-----------------------------------------------------------
#-----------------------------------------------------------

#-----------------------------------------------------------
# Solve tri-diagonal system
    # Solves for a vector u of size N the tridiagonal linear set given by equation (2.4.1) in NR
    # using a serial algorithm.
    # Input vectors b (diagonal elements) and r (right-hand sides) have size N,
    # while a and c (off-diagonal elements) are size N − 1.
function tridag(a,b,c,r)
    # Check that sizes agree
    n = length(a)+1 # Lenght of vectors
    if any( [length(b) length(c)+1 length(r)].!=n)
        error("Interpolation requires length(x)==length(y)")
    end
    # Check boundary conditions
    bet = b[1]
    if (bet == 0)
        error("tridag: Error at code stage 1")
        # If this happens then you should rewrite your equations as a set of order N − 1,
        # with u2 trivially eliminated.
    end
    # Solution of system
    u = zeros(n)
    g = zeros(n)
    u[1]=r[1]/bet
    for j=2:n
        g[j] = c[j-1]/bet
        bet  = b[j]-a[j-1]*g[j]
        if (bet == 0)
        error("tridag_ser: Error at code stage 2")
            # Decomposition and forward substitution.
            # Algorithm fails; see below routine in Vol. 1. of NR
        end
        u[j]=(r[j]-a[j-1]*u[j-1])/bet
    end
    for j=n-1:-1:1
        u[j]=u[j]-g[j+1]*u[j+1]
    end
    return u
end


#-----------------------------------------------------------
# Get second derivatives
    function spline_ypp(x,y,yp1=nothing,ypn=nothing)
        # Check that x grid is ordered
        if any(diff(x).<0)
            error("Grid for x must be increasing")
        end
        # Check that sizes agree
        if length(x)!=length(y)
            error("Interpolation requires length(x)==length(y)")
        end
        n = length(x) # Lenght of vectors
        # Set up tri-diagonal system
        a = Array{Float64}(undef,n)
        b = Array{Float64}(undef,n)
        c = Array{Float64}(undef,n)
        r = Array{Float64}(undef,n)
        # Fill in elements for tri-diagonal system
        c[1:n-1].= x[2:n].-x[1:n-1]
        r[1:n-1].= 6*((y[2:n].-y[1:n-1])./c[1:n-1])
        r[2:n-1].= r[2:n-1].-r[1:n-2]
        a[2:n-1].= c[1:n-2]
        b[2:n-1].= 2*(c[2:n-1].+a[2:n-1])
        b[1]     = 1
        b[n]     = 1
        # Lower Boundary
        if yp1==nothing # Use natural spline
            r[1] = 0
            a[1] = 0
        else # User supplied a derivative
            r[1] = (3/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1)
            c[1] = 0.5
        end
        # Upper Boundary
        if ypn==nothing # Use natural spline
            r[n] = 0
            c[n] = 0
        else # User supplied a derivative
            r[n] = (-3/(x[n]-x[n-1]))*((y[n]-y[n-1])/(x[n]-x[n-1])-ypn)
            a[n] = 0.5
        end
        # Compute second derivatives from tri-diagonal system
        ypp = tridag(a[2:n],b[1:n],c[1:n-1],r[1:n])
        return ypp
    end
#-----------------------------------------------------------

#-----------------------------------------------------------
# Get interpolation
    # Given the arrays x and y, which tabulate a function
    # (with the xi’s in increasing order),
    # and given the array ypp, which is the output from spline above,
    # and given a value of x, this routine returns a cubic-spline interpolated value.
    # The arrays x, y and ypp are all of the same size.
function spline_itp(x,y,ypp,z)
    # Check that z∈[x_min,x_max]
    if z<minimum(x) || z>maximum(x)
        error("Interpolation point z must be inside the grid bounds")
    end
    # Defien grid size and bracket z in grid
    n   = length(x) # Lenght of vectors
    khi = findmax(sign.(x .- z))[2] # Index of higher brakcet
    klo = khi-1                     # Index of lower braket
    h   = x[khi]-x[klo]             # Size of bracket
    # Define convexx weights
    a   = (x[khi]-z)/h
    b   = (z-x[klo])/h
    # Evaluate cubic spline
    itp = a*y[klo]+b*y[khi]+((a^3-a)*ypp[klo]+(b^3-b)*ypp[khi])*(h^2)/6
    return itp
end
#-----------------------------------------------------------

function spline_ditp(x,y,ypp,z)
    # Check that z∈[x_min,x_max]
    if z<minimum(x) || z>maximum(x)
        error("Interpolation point z must be inside the grid bounds")
    end
    # Defien grid size and bracket z in grid
    n   = length(x) # Lenght of vectors
    khi = findmax(sign.(x .- z))[2] # Index of higher brakcet
    klo = khi-1                     # Index of lower braket
    h   = x[khi]-x[klo]             # Size of bracket
    # Define convexx weights
    a   = (x[khi]-z)/h
    b   = (z-x[klo])/h
    # Evaluate cubic spline
    ditp = (y[khi]-y[klo])/(x[khi]-x[klo]) - ((3*a^2-1)*ypp[klo] + (3*b^2-1)*ypp[khi])*(x[khi]-x[klo])/6
    return ditp
end

function spline_NR(x,y,yp1=nothing,ypn=nothing)
    # Get second derivatives
    ypp = spline_ypp(x,y,yp1,ypn)
    # Define interpolation function
    F(z)  = spline_itp(x,y,ypp,z)
    # Define derivative interpolation function
    dF(z) = spline_ditp(x,y,ypp,z)
    return F,dF
end
#----

x=range(0.05,2,500000)

F,DF=spline_NR(gridc,gridu)

y_hat_spl=[F(xi) for xi in x]


asses=zeros(5,2)
for i=1:5
    y_hat=zeros(50000,2)
    for j=1:2
        gridc=Vector(range(0.05,2,length=gridnumber[i]))
        gridu=map(utility,gridc)
        x=range(0.05,2,50000) ### Grid to evaluate Interpolation
        F,DF=spline_NR(gridc,gridu)
        y_hat[:,1] = [LinearInt(xi, gridc, gridu) for xi in x]
        y_hat[:,2]=[F(xi) for xi in x]
        asses[i,j]=MSE1(y_hat[:,j],map(utility,x))
    end
end



############# Optim Curvature
function crra(c)
    σ=2
    u=(c^(1-σ))/(1-σ)
    return u
end



function optimcurv(N,lb,ub,θ)
    gridc=Vector(range(0,1,N))
    grid_values = collect(gridc)
    gridc = lb .+ (ub-lb) * grid_values .^ θ
    gridu=map(crra,gridc)
    F,DF=spline_NR(gridc,gridu)
    x=range(0.05,2,50000)
    y_hat_exp=[F(xi) for xi in x]
    mse=maximum(abs.(y_hat_exp-map(crra,x)))-10E-5
    return mse
end

using Roots
using Optim
optimcurv1=(θ->optimcurv(100,0.05,2,θ))
sol=optimize(optimcurv1,[0.5])  #### Correct Steady State of Labor if I use KSS
Optim.minimizer(sol)

##### Optimal θ is 1.197 #####

##### Point 3 (Extrapolation) with Gridsize=10 ####
#####(Incomplete)#####
