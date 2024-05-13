
function pT_projection(models;Npoints=200,colors=[:red,:blue,:green,:purple,:black],styles=[:solid,:dash, :dot, :dashdot, :dashdotdot])
    nmodels = length(models)
    npres = length(p)

    # Basic settings for the plot
    plt = plot(grid=:off,framestyle=:box,foreground_color_legend = nothing,legend_font=font(12))

    for i in 1:nmodels
        # Step 1: Saturation pressure of each component
        pures = split_model(models[i])
        crit = crit_pure(pures[j])
        tc = [crit[k][1] for k in 1:length(pures)]
        pc = [crit[k][2] for k in 1:length(pures)]
        for j in 1:2
            T .= LinRange(0.3*tc[j],tc[j],Npoints)

            sat = saturation_pressure.(models[i],T)
            psat = [sat[k][1] for k in 1:Npoints]
            plot!(plt,T,psat./1e5,color=colors[i],line = (:path, 3),label = false)        
            plot!(plt,[tc[j]],[pc[j]]./1e5,seriestype=:scatter,color=colors[i],markerstrokecolor=colors[i], line = (:scatter, 0.5),label = false)
        end

        # Step 2: Critical points of the mixture
        if tc[1]<tc[2]
            x = LinRange(0.,1.,Npoints)
        else
            x = LinRange(1.,0.,Npoints)
        end

        Tc = zeros(Npoints)
        Pc = zeros(Npoints)
        v0 = Clapeyron.x0_crit_mix(models[i],[x[1],1-x[1]])
        for k in 1:Npoints
            (Tc[k],Pc[k],vc) = crit_mix(models[i],[x[k],1-x[k]];v0=v0)
            v0 = [log10(vc),Tc[k]]
        end

        plot!(plt,Tc,Pc./1e5,color=colors[j],linestyle=styles[i],line = (:path, 2),label = false)
    end

    ylabel!(plt,"Pressure / bar",xguidefontsize=12)
    xlabel!(plt,"Temperature / K",yguidefontsize=12)
    ylims!(plt,(Tmin/1.1,Tmax*1.1))
    xlims!(plt,(0,1))

    return plt
end
    
export pT_projection