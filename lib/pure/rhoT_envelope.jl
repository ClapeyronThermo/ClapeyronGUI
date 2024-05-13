function rhoT_envelope(model;Npoints=200,color=:red,style=:solid)
    plt = plot(grid=:off,framestyle=:box,foreground_color_legend = nothing,legend_font=font(12))

    vmin = Inf
    vmax = 0.
    Tmin = Inf
    Tmax = 0.

    _rhoT_envelope!(plt,model,vmin,vmax,Tmin,Tmax;Npoints=Npoints,color=color,style=style)
    return plt
end

function rhoT_envelope!(plt,model;Npoints=200,color=:red,style=:solid)
    v = plt.series_list[1].plotattributes[:x_extrema]
    vmax = 1e-3/v[1]
    vmin = 1e-3/v[2]
    T = plt.series_list[1].plotattributes[:y_extrema]
    Tmax = T[2]
    Tmin = T[1]

    _rhoT_envelope!(plt,model,vmin,vmax,Tmin,Tmax;Npoints=Npoints,color=color,style=style)
    return plt
end

function _rhoT_envelope!(plt,model,vmin,vmax,Tmin,Tmax; Npoints=200, color=:red, style=:solid)
    vsat = zeros(2*Npoints)
    T = zeros(2*Npoints)

    (Tc,pc,Vc) = crit_pure(model)
    t = LinRange(0.3*Tc,Tc,Npoints)
    sat = saturation_pressure.(model,t)
    vsat[1:Npoints] .= [sat[i][2] for i in 1:Npoints]
    vsat[Npoints+1:2*Npoints] .= [sat[i][3] for i in Npoints:-1:1]

    T[1:Npoints] .= t
    T[Npoints+1:2*Npoints] .= reverse(t)

    plot!(plt,1e-3 ./vsat,T,color=color,line = (:path, 3),label = false)
    plot!(plt,[1e-3/Vc],[Tc],seriestype=:scatter,color=color,markerstrokecolor=color, line = (:scatter, 0.5),label = false)

    Tmin = maximum([T[1],Tmin])
    Tmax = maximum([Tc*1.05,Tmax])

    vmin = minimum([vsat[1]*0.95,vmin])
    vmax = maximum([vsat[end]*1.05,vmax])

    ylabel!(plt,"Temperature / K",xguidefontsize=12)
    xlabel!(plt,"Density / (mol/dm³)",yguidefontsize=12)
    xlims!(plt,(1e-3/vmax,1e-3/vmin))
    ylims!(plt,(Tmin,Tmax))
    return plt
end
    
export rhoT_envelope, rhoT_envelope!