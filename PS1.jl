
cd("/Users/juanholguin/Documents/GitHub/Computational-Macro")
raw"Defining Parameters Problem A"

alpha=1/3
z=1
beta=0.90

"Steady state"

kss=(alpha*beta*z)^(1/(1-alpha))

k=zeros(Float64, 20, 5)

k[1,1]=kss
k[2,1]=0.8*kss

for i in 3:20
    k[i,1]=alpha*beta*z*(k[i-1,1])^(alpha)
end

for i in 1:20
    k[i,2]=z*k[i,1]^alpha
end

for i in 2:19
    k[i,3]=(z*k[i,1]^alpha)-k[i+1,1]
    k[20,3]=k[19,3]
    k[1,3]=z^(1/1-alpha)*((alpha*beta)^(alpha/(1-alpha)))-((alpha*beta)^(1/(1-alpha)))
end

for i in 1:20
    k[i,4]=z*alpha*k[i,1]^(alpha-1)  
end

for i in 1:20
    k[i,5]=z*(1-alpha)*k[i,1]^(alpha)  
end

k


using Plots
gr()
plot(k[:,1],linewidth=1,marker=(:diamond,1),markercolor=RGB(0.1,0.1,0.1),label=["Capital"])
hline!([k[1,1]],line=(:dash),labels="SS")


savefig("a_k.png") 

using Plots
gr()
plot(k[:,2],linewidth=1,marker=(:diamond,1),markercolor=RGB(0.1,0.1,0.1),label=["Production"])
hline!([k[1,2]],line=(:dash),labels="SS")
savefig("a_y.png") 

using Plots
gr()
plot(k[:,3],linewidth=1,marker=(:diamond,1),markercolor=RGB(0.1,0.1,0.1),label=["Consumption"])
hline!([k[1,3]],line=(:dash),labels="SS")
savefig("a_c.png") 

using Plots
gr()
plot(k[:,4],linewidth=1,marker=(:diamond,1),markercolor=RGB(0.1,0.1,0.1),label=["r"])
hline!([k[1,4]],line=(:dash),labels="SS")
savefig("a_r.png") 

using Plots
gr()
plot(k[:,5],linewidth=1,marker=(:diamond,1),markercolor=RGB(0.1,0.1,0.1),label=["w"])
hline!([k[1,5]],line=(:dash),labels="SS")
savefig("a_w.png") 

"Defining Parameters Problem B"

alpha=1/3
z=1.05
beta=0.90

"Steady state"

kss=(alpha*beta)^(1/(1-alpha))

k=zeros(Float64, 20, 5)

k[1,1]=kss


for i in 2:20
    k[i,1]=alpha*beta*z*(k[i-1,1])^(alpha)
end

for i in 1:20
    k[i,2]=z*k[i,1]^alpha
end

for i in 2:19
    k[i,3]=(z*k[i,1]^alpha)-k[i+1,1]
    k[20,3]=k[19,3]
    k[1,3]=z^(1/1-alpha)*((alpha*beta)^(alpha/(1-alpha)))-((alpha*beta)^(1/(1-alpha)))
end

for i in 1:20
    k[i,4]=z*alpha*k[i,1]^(alpha-1)  
end

for i in 1:20
    k[i,5]=z*(1-alpha)*k[i,1]^(alpha)  
end


using Plots
gr()
plot(k[:,1],linewidth=1,marker=(:diamond,1),markercolor=RGB(0.1,0.1,0.1),label=["Capital"])
savefig("b_k.png") 

hline!([k[1,1]],line=(:dash),labels="Initial SS")
hline!([k[20,1]],line=(:dash),labels="Final SS")

using Plots
gr()
plot(k[:,2],linewidth=1,marker=(:diamond,1),markercolor=RGB(0.1,0.1,0.1),label=["Production"])
hline!([k[1,2]],line=(:dash),labels="Initial SS")
hline!([k[20,2]],line=(:dash),labels="Final SS")
savefig("b_y.png") 
using Plots
gr()
plot(k[:,3],linewidth=1,marker=(:diamond,1),markercolor=RGB(0.1,0.1,0.1),label=["Consumption"])
hline!([k[1,3]],line=(:dash),labels="Initial SS")
hline!([k[20,3]],line=(:dash),labels="Final SS")
savefig("b_c.png") 
using Plots
gr()
plot(k[:,4],linewidth=1,marker=(:diamond,1),markercolor=RGB(0.1,0.1,0.1),label=["r"])
hline!([k[1,4]],line=(:dash),labels="Initial SS")
hline!([k[20,4]],line=(:dash),labels="Final SS")
savefig("b_r.png") 
using Plots
gr()
plot(k[:,5],linewidth=1,marker=(:diamond,1),markercolor=RGB(0.1,0.1,0.1),label=["w"])
hline!([k[1,5]],line=(:dash),labels="Initial SS")
hline!([k[20,5]],line=(:dash),labels="Final SS")
savefig("b_w.png") 