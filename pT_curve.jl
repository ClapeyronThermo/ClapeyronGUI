function pT_curve(model;Npoints=200,color=:red, style=:solid)
    plt = Plots.plot(grid=:off,framestyle=:box,foreground_color_legend = nothing,legend_font=Plots.font(12), tickfontsize=12)

    pmax = 0.
    pmin = Inf
    Tmax = 0.
    Tmin = Inf

    _pT_curve!(plt,model,pmin,pmax,Tmin,Tmax;Npoints=Npoints,color=color,style=style)

    return plt
end

function pT_curve!(plt,model;Npoints=200,color=:red, style=:solid)
    p = plt.series_list[1].plotattributes[:y_extrema]
    pmax = p[2]
    pmin = p[1]

    T = plt.series_list[1].plotattributes[:x_extrema]
    Tmax = T[2]
    Tmin = T[1]

    _pT_curve!(plt,model,pmin,pmax,Tmin,Tmax;Npoints=Npoints,color=color,style=style)

    return plt
end

function _pT_curve!(plt,model,pmin,pmax,Tmin,Tmax;Npoints=200, color=:red, style=:solid)
    (Tc,Pc,vc) = crit_pure(model)
    T = LinRange(0.3*Tc,Tc,Npoints)
    sat = saturation_pressure.(model,T)
    psat = [sat[i][1] for i in 1:Npoints]
    
    Plots.plot!(plt,T,psat./1e5,color=color, style=style,line = (:path, 3),label = false)
    Plots.plot!(plt,[Tc],[Pc]./1e5,seriestype=:scatter,color=color,markerstrokecolor=color, line = (:scatter, 0.5),label = false)

    Tmin = minimum([minimum(T),Tmin])
    Tmax = maximum([Tc*1.05,Tmax])

    pmin = minimum([minimum(psat)./1e5,pmin])
    pmax = maximum([Pc*2.0./1e5,pmax])

    Plots.xlabel!(plt,"Temperature / K",xguidefontsize=12)
    Plots.ylabel!(plt,"Pressure / bar",yguidefontsize=12)
    Plots.ylims!(plt,(pmin,pmax))
    Plots.xlims!(plt,(Tmin,Tmax))
end

export pT_curve, pT_curve!