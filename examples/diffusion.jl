module DIFFUSION
using GenieFramework
using Clapeyron, Main.ThermoPlots
using SparseArrays, BandedMatrices, GCIdentifier, ChemicalIdentifiers
using ForwardDiff, OrdinaryDiffEq, UnPack, LinearAlgebra

import PlotlyBase, PlotlyJS, PlotlyKaleido
import Main.@timeout
@genietools

struct ScentDiffusionModel{T1, T2, T3, T4}
    ∂²::T1
    mixture::T2
    Diff::T3
    PVap::T4
end

cachey0 = [1.0; zeros(98)]
xmax = 2.0
h = xmax/100
domain = 0.0:h:100*h |> collect
C_T = 101325.0/(Clapeyron.R̄ * 298.15) #mol/m3
A_lg = 0.071 #m^2
    
function(f::ScentDiffusionModel)(du, u, p, t)

@unpack ∂², mixture, Diff, PVap = f

        y0 = activity_coefficient(mixture, 1.0, 298.15, [u[end, 1], u[end, 2], u[end, 3], u[end, 4]]).*PVap/101325.0.*
        u[end, :]./(u[end, 1] + u[end, 2] + u[end, 3] + u[end, 4])

        for i in 1:4 
        # Gas phase balance
        du[1:end - 1, i] = (Diff[i]) * (∂²*u[1:end - 1, i] + y0[i]*cachey0)/h^2 #Last parts account for boundary value
        
        # Liquid phase balance
        du[end , i] = (Diff[i]) * A_lg * C_T * (-u[2, i] + 4*u[1, i] - 3*y0[i])/(2*h)
        #du[end , i] = Diff[i] * A_lg * C_T * (u[1, i] - y0[i])/h
        
        end
end

@app begin
    @in x1 = 0.02
    @in x2 = 0.03
    @in x3 = 0.05
    @in simulate_button = false

    @out trace = []
    @out layout = PlotlyBase.Layout(autosize=false,width=700,height=470,
                            xaxis = PlotlyBase.attr(title = "Distance  / m", font_size=12, showgrid=false,            
                                          ticks="inside",mirror=true,showline=true,linecolor="black"),
                            yaxis = PlotlyBase.attr(title = "Concentration / (mol/dm³)", font_size=12, showgrid=false,       
                                          ticks="inside",mirror=true,showline=true,linecolor="black"),
                             showlegend=false, plot_bgcolor="white")

    @onbutton simulate_button begin
        species = ["limonene", "geraniol", "vanillin", "ethanol"]
        components = get_groups_from_name.(species,UNIFAC)
        mixture = UNIFAC(components;  puremodel = BasicIdeal)
        x₀ = [x1, x2, x3, 1-x1-x2-x3] #Limonene, geraniol, vanillin, ethanol
        N₀ = 1e-3 #mmol
        L = 2 # m

        #stencil
        n = 100
        ∂² = BandedMatrix{Float64}(undef, (n + 1, n + 1), (1,1))
        ∂²[band(0)] .= -2
        ∂²[band(1)] .= 1
        ∂²[band(-1)] .= 1
        ∂² = ∂²[2:end - 1, 2:end - 1]
        u0 = zeros(size(∂², 2) + 1, 4) .+ 0.0
        n₀x₀ = N₀*x₀ 
        u0[end, :] = n₀x₀


        D_quaternary = [6.010317767518301e-6, 5.822748524332796e-6, 6.803670097258762e-6, 1.2341364039719445e-5]
        PVap_quaternary = [277.0013821460195, 4.124352858268832, 0.053095574752943624, 7806.831390542604]
        C_T = 101325.0/(Clapeyron.R̄ * 298.15) #mol/m3
        A_lg = 0.071 #m^2
        L = 2.0 #m

        #Testing solution
        rhs_k = ScentDiffusionModel(∂², mixture, D_quaternary, PVap_quaternary)
        tspan_k = (0.0, 3600.0*10)
        #tspan_k = (0.0, 60.0*10)
        p = zeros(8)
        problem_k = ODEProblem(rhs_k, u0, tspan_k, p)
        solution_k = solve(problem_k, FBDF(autodiff = true),
        saveat = 100., abstol = 1e-8, reltol = 1e-8);

        solarray_k = Array(solution_k)

        for i in 1:size(solarray_k, 3)
            trace = [PlotlyBase.scatter(x = domain[1:end-2], y = solarray_k[1:end-1, j, i]*C_T, mode = "lines", line = PlotlyBase.attr(width = 2), name = species[j]) for j in 1:4]
            layout = PlotlyBase.Layout(autosize=false,width=700,height=470,
            xaxis = PlotlyBase.attr(title = "Distance  / m", font_size=12, showgrid=false,            
                          ticks="inside",mirror=true,showline=true,linecolor="black", range=[0, L]),
            yaxis = PlotlyBase.attr(title = "Concentration / (mol/m³)", font_size=12, showgrid=false,       
                          ticks="inside",mirror=true,showline=true,linecolor="black", range=[0, 0.02]),
             showlegend=true, plot_bgcolor="white")
            sleep(0.01)
        end
    end
end
end