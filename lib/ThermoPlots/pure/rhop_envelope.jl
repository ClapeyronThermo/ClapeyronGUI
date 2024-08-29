function rhop_envelope(model;Npoints=200,color="red", style="solid")
    layout = Layout(autosize=false,width=700,height=470,
    xaxis = attr(title = "Density / (mol/dm³)", font_size=12, showgrid=false, ticks="inside",mirror=true,showline=true,linecolor="black"),
    yaxis = attr(title = "Pressure / bar", font_size=12, showgrid=false, ticks="inside",mirror=true,showline=true,linecolor="black"),
    showlegend=false, plot_bgcolor="white")
    plt = plot(scatter(),layout)
    pmin = Inf
    pmax = 0.
    vmin = Inf
    vmax = 0.

    _rhop_envelope!(plt,model,vmin,vmax,pmin,pmax;Npoints=Npoints,color=color,style=style) 

    return plt
end

function rhop_envelope!(plt,model;Npoints=200,color="red", style="solid")
    p = plt.layout[:yaxis][:range]
    pmax = p[2]
    pmin = p[1]

    v = plt.layout[:xaxis][:range]
    vmax = 1e-3/v[1]
    vmin = 1e-3/v[2]

    _rhop_envelope!(plt,model,vmin,vmax,pmin,pmax;Npoints=Npoints,color=color,style=style) 

    return plt

end

function _rhop_envelope!(plt,model,vmin,vmax,pmin,pmax;Npoints=200,color="red", style="solid")
    vsat = zeros(2*Npoints)
    psat = zeros(2*Npoints)

    (Tc,Pc,Vc) = crit_pure(model)
    T = LinRange(0.3*Tc,Tc,Npoints)
    sat = saturation_pressure.(model,T)
    vsat[1:Npoints] .= [sat[i][2] for i in 1:Npoints]
    vsat[Npoints+1:2*Npoints] .= [sat[i][3] for i in Npoints:-1:1]

    psat[1:Npoints] .= [sat[i][1] for i in 1:Npoints]
    psat[Npoints+1:2*Npoints] .= [sat[i][1] for i in Npoints:-1:1]

    T = T[.!isnan.(psat)]
    vsat = vsat[.!isnan.(psat)]
    psat = psat[.!isnan.(psat)]

    line_sat = scatter(x=1e-3 ./vsat,y=psat./1e5,mode="lines",line=attr(color=color, dash=style, width=3),name="VLE curve")
    scatter_sat = scatter(x=[1e-3/Vc],y=[Pc/1e5],mode="markers",marker=attr(color=color, size=6),name="Critical point")
    addtraces!(plt,line_sat,scatter_sat)

    Pmin = minimum([psat[1]./1e5,pmin])
    Pmax = maximum([Pc*1.1./1e5,pmax])

    update!(plt,layout=Layout(yaxis = attr(range = [Pmin,Pmax])))

    vmin = minimum([vsat[1]*0.95,vmin])
    vmax = maximum([vsat[end]*1.05,vmax])

    update!(plt,layout=Layout(xaxis = attr(range = [1e-3/vmax,1e-3/vmin])))
end
    
export rhop_envelope, rhop_envelope!