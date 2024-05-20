function isothermal_bulk(models,property,T,p_lims,z=[1.];phase=:unknown,Npoints=200,colors=[:red,:blue,:green,:purple,:black],styles=[:solid,:dash, :dot, :dashdot, :dashdotdot])
    nmodels = length(models)
    ntemps = length(T)
    plt = plot(grid=:off,framestyle=:box,foreground_color_legend = nothing,legend_font=font(12))

    p = LinRange(p_lims[1],p_lims[2],Npoints)
    Y_max = -Inf
    Y_min = Inf
    for i in 1:nmodels
        for j in 1:ntemps
            Y = property.(models[i],p,T[j],Ref(z);phase=phase)
            plot!(plt,p./1e5,Y,color=colors[j],linestyle=styles[i],line = (:path, 2),label = false)
            Y_max = maximum((Y_max,maximum(Y)))
            Y_min = minimum((Y_min,minimum(Y)))
        end
    end

    xlabel!(plt,"Pressure / bar",xguidefontsize=12)
    xlims!(plt,(p_lims[1],p_lims[2]))
    ylims!(plt,(Y_min,Y_max))

    return plt
end

function isobaric_bulk(models,property,p,T_lims,z=[1.];phase=:unknown,Npoints=200,colors=[:red,:blue,:green,:purple,:black],styles=[:solid,:dash, :dot, :dashdot, :dashdotdot])
    nmodels = length(models)
    npres = length(p)
    plt = plot(grid=:off,framestyle=:box,foreground_color_legend = nothing,legend_font=font(12))

    T = LinRange(T_lims[1],T_lims[2],Npoints)
    Y_max = -Inf
    Y_min = Inf
    for i in 1:nmodels
        for j in 1:npres
            Y = property.(models[i],p[j],T,Ref(z);phase=phase)
            plot!(plt,T,Y,color=colors[j],linestyle=styles[i],line = (:path, 2),label = false)
            Y_max = maximum((Y_max,maximum(Y)))
            Y_min = minimum((Y_min,minimum(Y)))
        end
    end

    xlabel!(plt,"Temperature / K",xguidefontsize=12)
    xlims!(plt,(T_lims[1],T_lims[2]))
    ylims!(plt,(Y_min,Y_max))

    return plt
end
    
export isobaric_bulk, isothermal_bulk