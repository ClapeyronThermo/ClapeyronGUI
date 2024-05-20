function rhoT_envelope(model;Npoints=200,color="red", style="solid")
    layout = Layout(autosize=false,width=700,height=470,
    xaxis = attr(title = "Density / (mol/dm³)", font_size=12, showgrid=false, ticks="inside",mirror=true,showline=true, linecolor="black"),
    yaxis = attr(title = "Temperature / K", font_size=12, showgrid=false, ticks="inside",mirror=true,showline=true, linecolor="black"),
    showlegend=false, plot_bgcolor="white")
    plt = plot(scatter(),layout)
    vmin = Inf
    vmax = 0.
    Tmin = Inf
    Tmax = 0.

    _rhoT_envelope!(plt,model,vmin,vmax,Tmin,Tmax;Npoints=Npoints,color=color,style=style)
    return plt
end

function rhoT_envelope!(plt,model;Npoints=200,color="red", style="solid")
    v = plt.plot.layout[:xaxis][:range]
    vmax = 1e-3/v[1]
    vmin = 1e-3/v[2]
    T = plt.plot.layout[:yaxis][:range]
    Tmax = T[2]
    Tmin = T[1]

    _rhoT_envelope!(plt,model,vmin,vmax,Tmin,Tmax;Npoints=Npoints,color=color,style=style)
    return plt
end

function _rhoT_envelope!(plt,model,vmin,vmax,Tmin,Tmax; Npoints=200,color="red", style="solid")
    vsat = zeros(2*Npoints)
    T = zeros(2*Npoints)

    (Tc,pc,Vc) = crit_pure(model)
    t = LinRange(0.3*Tc,Tc,Npoints)
    sat = saturation_pressure.(model,t)
    vsat[1:Npoints] .= [sat[i][2] for i in 1:Npoints]
    vsat[Npoints+1:2*Npoints] .= [sat[i][3] for i in Npoints:-1:1]

    T[1:Npoints] .= t
    T[Npoints+1:2*Npoints] .= reverse(t)

    T = T[.!isnan.(vsat)]
    vsat = vsat[.!isnan.(vsat)]

    line_sat = scatter(x=1e-3 ./vsat,y=T,mode="lines",line=attr(color=color, dash=style, width=3),name="VLE curve")
    scatter_sat = scatter(x=[1e-3/Vc],y=[Tc],mode="markers",marker=attr(color=color, size=6),name="Critical point")
    addtraces!(plt,line_sat,scatter_sat)
    
    Tmin = maximum([T[1],Tmin])
    Tmax = maximum([Tc*1.05,Tmax])

    update!(plt,layout=Layout(yaxis = attr(range = [Tmin,Tmax])))

    vmin = minimum([vsat[1]*0.95,vmin])
    vmax = maximum([vsat[end]*1.05,vmax])

    update!(plt,layout=Layout(xaxis = attr(range = [1e-3/vmax,1e-3/vmin])))
    return plt
end
    
export rhoT_envelope, rhoT_envelope!