function pT_curve(model;Npoints=200,color="red", style="solid")
    # plt = plot(grid=:off,framestyle=:box,foreground_color_legend = nothing,legend_font=font(12))
    layout = Layout(xaxis = attr(title = "Temperature  / K", font_size=12, showgrid=false, ticks="inside",mirror=true,showline=true,linecolor="black"),
                    yaxis = attr(title = "Pressure / bar", font_size=12, showgrid=false, ticks="inside",mirror=true,showline=true,linecolor="black"),
                    showlegend=false, plot_bgcolor="white")
    plt = plot(scatter(),layout)

    pmax = 0.
    pmin = Inf
    Tmax = 0.
    Tmin = Inf

    _pT_curve!(plt,model,pmin,pmax,Tmin,Tmax;Npoints=Npoints,color=color,style=style)

    return plt
end

function pT_curve!(plt,model;Npoints=200,color="red", style="solid")
    p = plt.plot.layout[:yaxis][:range]
    pmax = p[2]
    pmin = p[1]

    T = plt.plot.layout[:xaxis][:range]
    Tmax = T[2]
    Tmin = T[1]

    _pT_curve!(plt,model,pmin,pmax,Tmin,Tmax;Npoints=Npoints,color=color,style=style)

    return plt
end

function _pT_curve!(plt,model,pmin,pmax,Tmin,Tmax;Npoints=200, color="red", style="solid")
    (Tc,Pc,vc) = crit_pure(model)
    T = LinRange(0.3*Tc,Tc,Npoints)
    sat = saturation_pressure.(model,T)
    psat = [sat[i][1] for i in 1:Npoints]

    T = T[.!isnan.(psat)]
    psat = psat[.!isnan.(psat)]
    
    line_sat = scatter(x=T,y=psat./1e5,mode="lines",line=attr(color=color, dash=style, width=3),name="Saturation curve")
    scatter_sat = scatter(x=[Tc],y=[Pc/1e5],mode="markers",marker=attr(color=color, size=6),name="Critical point")
    addtraces!(plt,line_sat,scatter_sat)

    Tmin = minimum([minimum(T),Tmin])
    Tmax = maximum([Tc*1.05,Tmax])

    update!(plt,layout=Layout(xaxis = attr(range = [Tmin,Tmax])))

    pmin = minimum([minimum(psat)./1e5,pmin])
    pmax = maximum([Pc*1.1./1e5,pmax])

    update!(plt,layout=Layout(yaxis = attr(range = [pmin,pmax])))
end

export pT_curve, pT_curve!