function rhop_envelope(model;Npoints=200,color=:red,style=:solid)
    plt = plot(grid=:off,framestyle=:box,foreground_color_legend = nothing,legend_font=font(12))

    pmin = Inf
    pmax = 0.
    vmin = Inf
    vmax = 0.

    _rhop_envelope!(plt,model,vmin,vmax,pmin,pmax;Npoints=Npoints,color=color,style=style) 

    return plt
end

function rhop_envelope!(plt,model;Npoints=200,color=:red,style=:solid)
    p = plt.series_list[1].plotattributes[:y_extrema]
    pmax = p[2]
    pmin = p[1]

    v = plt.series_list[1].plotattributes[:x_extrema]
    vmax = 1e-3/v[1]
    vmin = 1e-3/v[2]

    _rhop_envelope!(plt,model,vmin,vmax,pmin,pmax;Npoints=Npoints,color=color,style=style) 

    return plt

end

function _rhop_envelope!(plt,model,vmin,vmax,pmin,pmax;Npoints=200,color=:red,style=:solid)
    vsat = zeros(2*Npoints)
    psat = zeros(2*Npoints)

    (Tc,Pc,Vc) = crit_pure(model)
    T = LinRange(0.3*Tc,Tc,Npoints)
    sat = saturation_pressure.(model,T)
    vsat[1:Npoints] .= [sat[i][2] for i in 1:Npoints]
    vsat[Npoints+1:2*Npoints] .= [sat[i][3] for i in Npoints:-1:1]

    psat[1:Npoints] .= [sat[i][1] for i in 1:Npoints]
    psat[Npoints+1:2*Npoints] .= [sat[i][1] for i in Npoints:-1:1]

    
    plot!(plt,1e-3 ./vsat,psat./1e5,color=color,line = (:path, 3),label = false)
    plot!(plt,[1e-3/Vc],[Pc]./1e5,seriestype=:scatter,color=color,markerstrokecolor=color, line = (:scatter, 0.5),label = false)

    Pmin = minimum([psat[1]./1e5,pmin])
    Pmax = maximum([Pc*2.0./1e5,pmax])

    vmin = minimum([vsat[1]*0.95,vmin])
    vmax = maximum([vsat[end]*1.05,vmax])

    ylabel!(plt,"Pressure / bar",xguidefontsize=12)
    xlabel!(plt,"Density / (mol/dm³)",yguidefontsize=12)
    xlims!(plt,(1e-3/vmax,1e-3/vmin))
    ylims!(plt,(Pmin,Pmax))

end
    
export rhop_envelope, rhop_envelope!