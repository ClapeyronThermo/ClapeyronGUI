import PlotlyBase

function pT_projection(model;Npoints=200,color=:red,style=:solid,check_ucep=false)
    layout = PlotlyBase.Layout(autosize=false,width=700,height=470,
                    xaxis = PlotlyBase.attr(title = "Temperature / K", font_size=12, showgrid=false, ticks="inside",mirror=true,showline=true,linecolor="black"),
                    yaxis = PlotlyBase.attr(title = "Pressure / bar", font_size=12, showgrid=false, ticks="inside",mirror=true,showline=true,linecolor="black"),
                    showlegend=false, plot_bgcolor="white")
    plt = PlotlyBase.Plot(PlotlyBase.scatter(),layout)
    Tmin = Inf
    Tmax = 0.

    pmin = Inf
    pmax = 0.

    _pT_projection!(plt,model,Tmin,Tmax,pmin,pmax;Npoints=Npoints,color=color,style=style,check_ucep=check_ucep)
    return plt
end

function pT_projection!(plt,model;Npoints=200,color=:red,style=:solid,check_ucep=false)
    p = plt.layout[:yaxis][:range]
    pmax = p[2]
    pmin = p[1]

    T = plt.layout[:xaxis][:range]
    Tmax = T[2]
    Tmin = T[1]

    _pT_projection!(plt,model,Tmin,Tmax,pmin,pmax;Npoints=Npoints,color=color,style=style,check_ucep=check_ucep)
    return plt
end

function _pT_projection!(plt,model,Tmin,Tmax,pmin,pmax;Npoints=200,check_ucep=false,color=:red,style=:solid)
    present_ucep = false
        # Step 1: Saturation pressure of each component
        pures = split_model(model)
        crit = crit_pure.(pures)
        tc = [crit[k][1] for k in 1:length(pures)]
        pc = [crit[k][2] for k in 1:length(pures)]
        for j in 1:2
            T = LinRange(0.6*tc[j],tc[j],Npoints)

            sat = saturation_pressure.(pures[j],T)
            psat = [sat[k][1] for k in 1:Npoints]

            # remove nans
            T = T[.!isnan.(psat)]
            psat = psat[.!isnan.(psat)]
            line_sat = PlotlyBase.scatter(x=T,y=psat./1e5,mode="lines",line=PlotlyBase.attr(color=color, dash=style, width=3),name="Saturation curve")
            scatter_sat = PlotlyBase.scatter(x=[tc[j]],y=[pc[j]/1e5],mode="markers",marker=PlotlyBase.attr(color=color, size=6),name="Critical point")
            PlotlyBase.addtraces!(plt,line_sat,scatter_sat)

            Tmin = min(Tmin,minimum(T))
            Tmax = max(Tmax,maximum(T))

            pmin = min(pmin,minimum(psat))
            pmax = max(pmax,maximum(psat))
        end

        # Step 2: Critical points of the mixture
        x = zeros(Npoints,2)
        if tc[1]<tc[2]
            x[:,1] = LinRange(0.,1.,Npoints)
            x[:,2] = LinRange(1.,0.,Npoints)
        else
            x[:,1] = LinRange(1.,0.,Npoints)
            x[:,2] = LinRange(0.,1.,Npoints)
        end

        npaths = 2
        Tc = zeros(Npoints,2)
        Pc = zeros(Npoints,2)
        v0 = Clapeyron.x0_crit_mix(model,[x[1],1-x[1]])
        i = 1
        while i <= npaths
            idxend = Npoints
            for k in 1:Npoints
                (Tc[k,i],Pc[k,i],vc) = crit_mix(model,[x[k,i],1-x[k,i]];v0=v0)
                v0 = [log10(vc),Tc[k,i]]
                if Pc[k,i]<= minimum(pc) && check_ucep && !present_ucep
                    present_ucep=true

                    println("Mixture is not type I. Will need to search for UCEP")
                elseif Pc[k,i] <= 0. 
                    idxend = k-1
                    
                    break
                elseif k > 2
                    if (Pc[k,i]-Pc[k-1,i])/(Pc[k-1,i]-Pc[k-2,i]) < 0  && (Tc[k,i]-Tc[k-1,i])/(Tc[k-1,i]-Tc[k-2,i]) < 0
                        idxend = k-1
                        break
                    end
                elseif k == Npoints
                    println("Mixture is type I or II.")
                    npaths = 1
                    break
                end
            end
            pmin = min(pmin,minimum(Pc[1:idxend,i]))
            pmax = max(pmax,maximum(Pc[1:idxend,i]))

            Tmin = min(Tmin,minimum(Tc[1:idxend,i]))
            Tmax = max(Tmax,maximum(Tc[1:idxend,i]))
            i += 1
        end

        

        for i in 1:npaths
            line_crit = PlotlyBase.scatter(x=Tc[Tc[:,i].>0,i],y=Pc[Tc[:,i].>0,i]./1e5,mode="lines",line=PlotlyBase.attr(color=color, dash="dash", width=2),name="Critical curve")
            PlotlyBase.addtraces!(plt,line_crit)
        end

        # if present_ucep && check_ucep
        #     try 
        #         ucep = UCEP_mix(model)
        #         if abs(ucep[3]-ucep[4])/ucep[3] > 1e-2 && ucep[1]<maximum(tc)
        #             plot!(plt,[ucep[1]],[ucep[2]./1e5],seriestype=:scatter,color=color,markerstrokecolor=color, line = (:scatter, 0.5),label = false)
        #             println("UCEP found. Now checking for VLLE.")
        #             T_VLLE = LinRange(min(ucep[1]-10,Tmin),ucep[1],Npoints)
        #             p_VLLE = zeros(Npoints)
        #             v0 = Clapeyron.x0_VLLE_pressure(model,Tmin)
        #             for k in 1:Npoints
        #                 VLLE = VLLE_pressure(model,T_VLLE[k];v0=v0)
        #                 p_VLLE[k] = VLLE[1]
        #                 v0 = [log10(VLLE[2]),log10(VLLE[3]),log10(VLLE[4]),VLLE[5][1],VLLE[6][1],VLLE[7][1]]
        #             end
        #             plot!(plt,T_VLLE,p_VLLE./1e5,color=color,linestyle=style,line = (:path, 2),label = false)
        #         else 
        #             println("UCEP not found")
        #         end
        #     catch
        #         println("UCEP not found")
        #     end
        # end

        PlotlyBase.update!(plt,layout=PlotlyBase.Layout(xaxis = PlotlyBase.attr(range = [Tmin/1.1,Tmax*1.1])))
        PlotlyBase.update!(plt,layout=PlotlyBase.Layout(yaxis = PlotlyBase.attr(range = [pmin/1.1/1e5,pmax*1.1/1e5])))
    return plt
end
    
export pT_projection