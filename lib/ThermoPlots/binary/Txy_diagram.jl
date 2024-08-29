# TO DOs for Txy_diagram
# 1. Add check for VLLE when it does not arise from an azeotrope
# 2. Add critical point evaluators
# 3. Add the case where both components are supercritical
import PlotlyBase
function Txy_diagram(model,p;iscrit=nothing,check_lle=false,lle_present=false,Npoints=200,color=:red,style=:solid)
    layout = PlotlyBase.Layout(autosize=false,width=700,height=470,
    xaxis = PlotlyBase.attr(title = "Molar composition of "*model.components[1], font_size=12, showgrid=false, ticks="inside",mirror=true,showline=true,linecolor="black"),
    yaxis = PlotlyBase.attr(title = "Temperature / K", font_size=12, showgrid=false, ticks="inside",mirror=true,showline=true,linecolor="black"),
    showlegend=false, plot_bgcolor="white")
plt = PlotlyBase.Plot(PlotlyBase.scatter(),layout)

    Tmax = 0.
    Tmin = Inf

    _Txy_diagram(plt,model,p,Tmax,Tmin;iscrit=iscrit,check_lle=check_lle,lle_present=lle_present,Npoints=Npoints,color=color,style=style)
    return plt
end

function Txy_diagram!(plt,model,p;iscrit=nothing,check_lle=false,lle_present=false,Npoints=200,color=:red,style=:solid)
    T = plt.layout[:yaxis][:range]
    Tmax = T[2]
    Tmin = T[1]
    _Txy_diagram(plt,model,p,Tmax,Tmin;iscrit=iscrit,check_lle=check_lle,lle_present=lle_present,Npoints=Npoints,color=color,style=style)
    return plt
end

function _Txy_diagram(plt,model,p,Tmax,Tmin;iscrit=nothing,check_lle=false, check_vlle = false,lle_present=false,Npoints=200,color=:red,style=:solid)
    T_lle = []
    Klle = []
    zlle = [0.5,0.5]
    stop_vlle = false

        if isnothing(iscrit)
            pures = split_model(model)
            crit = crit_pure.(pures)
            pc = [crit[k][2] for k in 1:length(pures)]
            crit = pc.<p
        else
            crit = iscrit[j]
        end

        # Check whether or not none, one or both components are supercritical
        if all(crit)
            @error "Critical point is below the temperature of interest for all species"
        elseif any(crit)
            x = zeros(Npoints,1)
            # Depending on which component is supercritical, we scan in different directions
            if crit[1]
                x[:,1] = LinRange(0.,1.,Npoints)
            else
                x[:,1] = LinRange(1.,0.,Npoints)
            end
            v0 = [Clapeyron.x0_bubble_temperature(model,p,[x[1],1-x[1]])]
            npaths = 1
        else
            sat = saturation_temperature.(pures,p)
            x = ones(Npoints,2)
            # When both components are subcritical, we need to scan in both directions
            # This is because there is a possibility that there are two critical points at a given temperature
            if sat[1][1]>sat[2][1]
                x[:,1] = LinRange(1e-6,1.0-1e-6,Npoints)
                x[:,2] = LinRange(1.0-1e-6,1e-6,Npoints)
            else
                x[:,1] = LinRange(1.0-1e-6,1e-6,Npoints)
                x[:,2] = LinRange(1e-6,1.0-1e-6,Npoints)
            end

            v0 = [Clapeyron.x0_bubble_temperature(model,p,[x[1,k],1-x[1,k]]) for k in 1:2]
            npaths = 2
        end


        # Scan through the VLE region
        y = zeros(Npoints,npaths)
        T = zeros(Npoints,npaths)
        for l in 1:npaths
            idxend = Npoints
            dx = x[2,l]-x[1,l]
            k=1
            while k <= Npoints

                # Obtain the bubble points
                bub = bubble_temperature(model,p,[x[k,l],1-x[k,l]];v0=v0[l])
                v0[l] = [bub[1],log10(bub[2]),log10(bub[3]),bub[4][1]+dx,bub[4][2]-dx]
                T[k,l] = bub[1]
                y[k,l] = bub[4][1]

                if k>1
                    # If there's a sudden change in direction along the x-y and p space, then there's VLLE
                    # This step can only help identify VLLE when it is also part of an azeotrope i.e. x < y < xx
                    # We need something to handle x < xx < y
                    if (y[k,l]-y[k-1,l])/dx < 0 && (T[k,l] - T[k-1,l])/(y[k,l]-y[k-1,l]) > 0 && !stop_vlle
                        println("Vapour-Liquid-Liquid equilibrium is present")
                        x0 = [x[k,l],1-x[k,l]]
                        y0 = bub[4]
                        xx0 = 1 .- x0

                        vll0 = log10(Clapeyron.volume(model,p,bub[1],xx0;phase=:l))
                        v0_dew = [bub[1], vll0, log10(bub[3]), xx0[1], xx0[2]]
                        dew = dew_temperature(model,p,y0;v0=v0_dew)
                        
                        T0 = (dew[1]+bub[1])/2

                        Klle = dew[4]./x0
                        lle = tp_flash(model,p,T0,(x0.+dew[4])./2,MichelsenTPFlash(K0=Klle,equilibrium=:lle))
                        Kvle = bub[4]./x0
                        vle = tp_flash(model,p,T0,(x0.+y0)./2,MichelsenTPFlash(K0=Kvle,equilibrium=:vle))

                        x0 = lle[1][1,1]
                        y0 = vle[1][2,1]
                        xx0 = lle[1][2,1]

                        vl0 = Clapeyron.volume(model,p,T0,[x0,1-x0];phase=:l)
                        vll0 = Clapeyron.volume(model,p,T0,[xx0,1-xx0];phase=:l)
                        vv0 = Clapeyron.volume(model,p,T0,[y0,1-y0];phase=:v)

                        v0_vlle = [T0,log10(vl0),log10(vll0),log10(vv0),x0,xx0,y0]

                        vlle = VLLE_temperature(model,p;v0=v0_vlle)
                        T_vlle = vlle[1]
                        x_vlle = vlle[5][1]
                        xx_vlle = vlle[6][1]
                        y_vlle = vlle[7][1]

                        line_vlle = PlotlyBase.scatter(x=sort([xx_vlle,x_vlle,y_vlle]),y=[T_vlle,T_vlle,T_vlle],mode="lines",line=PlotlyBase.attr(color=color, dash=style, width=3),name="VLLE curve")
                        PlotlyBase.addtraces!(plt,line_vlle)

                        idx_vlle = sum(T[1:k,l].>T_vlle)
                        T[idx_vlle+1,l] = T_vlle
                        x[idx_vlle+1,l] = x_vlle
                        y[idx_vlle+1,l] = y_vlle

                        T[idx_vlle+2,l] = T_vlle
                        x[idx_vlle+2,l] = xx_vlle
                        y[idx_vlle+2,l] = y_vlle

                        T[idx_vlle+3:end,l] .= 0
                        x[idx_vlle+3:end,l] .= LinRange(xx_vlle+dx,x[end,l],Npoints-idx_vlle-2)
                        y[idx_vlle+3:end,l] .= 0
                        k = idx_vlle+2

                        v0[l] = [vlle[1],log10(vlle[3]),log10(vlle[4]),vlle[7][1]+dx,vlle[7][2]-dx]

                        # Given there is VLLE, there must also be LLE
                        lle_present = true
                        check_lle = true
                        zlle = (xx_vlle+x_vlle)/2
                        zlle = [zlle, 1-zlle]
                        Klle = [xx_vlle/x_vlle,(1-xx_vlle)/(1-x_vlle)]
                        T_lle = LinRange(T_vlle,maximum((250,T_vlle-100)),Npoints)
                        stop_vlle = true
                    end
                end

                # If the Clapeyron.volumes of the two phases are similar or get reversed, we've reached a critical point
                # We can find that critical point and break the while loop and go to the next branch
                if bub[3]-bub[2]<1e-6
                    # x0 = (x[k-1,l]+y[k-1,l])/2
                    # v0_crit = [log10((bub[2]+bub[3])/2),[x0,1-x0]]
                    # crit_point = UCST_mix(model,T[j];v0=v0_crit)
                    # x[k,l] = crit_point[3][1]
                    # y[k,l] = crit_point[3][1]
                    # p[k,l] = crit_point[1]
                    idxend = k
                    break
                end
                k+=1
            end

            # Find the bounds of the plot
            Tmax = maximum(vcat(Tmax,T[:,l]))
            Tmin = minimum(vcat(Tmin,T[:,l]))

            # If we need to check for LLE, we need to look at higher pressures.
            # For visual clarity, we only really need to consider LLE when it takes up half the phase diagram
            if check_lle && !lle_present
                T_test = maximum([Tmin-50,250])
                x_test = [0.5,0.5]
                tpd = Clapeyron.tpd(model,p,T_test,x_test)
                if length(tpd[1])>1
                    Klle = tpd[1][1]./tpd[1][2]
                    lle_present=true
                    T_lle = LinRange(T_test,Tmin,Npoints)
                end
            end

            # If LLE is present, use a flash algorithm to find the LLE curve
            if lle_present && check_lle
                println("Liquid-Liquid Equilibria is present")
                idxend = Npoints
                x_lle= zeros(Npoints)
                xx_lle = zeros(Npoints)
                for k in 1:Npoints
                    flash = tp_flash(model,p,T_lle[k],zlle,MichelsenTPFlash(equilibrium=:lle))
                    x_lle[k] = flash[1][1,1]
                    xx_lle[k] = flash[1][2,1]
                    Klle = flash[1][1,:]./flash[1][2,:]
                    zlle = (flash[1][1,:]+flash[1][2,:])./2
                    # If the compositions of the two phase become similar, we have reached a critical point
                    if abs(x_lle[k]-xx_lle[k])<1e-4
                        idxend = k
                        # x0 = (x_lle[k,l]+xx_lle[k,l])/2
                        # vl0 = Clapeyron.volume(model,p_lle[k],T[j],[x0,1-x0];phase=:l)
                        # v0_crit = [log10(vl0),[x0,1-x0]]
                        # crit_point = UCST_mix(model,T[j];v0=v0_crit)
                        # x_lle[k,l] = crit_point[3][1]
                        # xx_lle[k,l] = crit_point[3][1]
                        # p_lle[k,l] = crit_point[1]
                        break
                    end
                end
                line_lle1 = PlotlyBase.scatter(x=x_lle[1:idxend],y=T_lle[1:idxend],mode="lines",line=PlotlyBase.attr(color=color, dash=style, width=3),name="LLE curve")
            line_lle2 = PlotlyBase.scatter(x=xx_lle[1:idxend],y=T_lle[1:idxend],mode="lines",line=PlotlyBase.attr(color=color, dash=style, width=3),name="LLE curve")
            PlotlyBase.addtraces!(plt,line_lle1,line_lle2)
                check_lle = false
                Tmin = minimum(vcat(Tmin,T_lle))
            end

            if idxend == Npoints
                X = vcat(x[:,l],reverse(y[:,l]))
                Y = vcat(T[:,l],reverse(T[:,l]))
                line_vle = PlotlyBase.scatter(x=X,y=Y,mode="lines",line=PlotlyBase.attr(color=color, dash=style, width=3),name="VLE curve")
                PlotlyBase.addtraces!(plt,line_vle)
                break
            else
                X = vcat(x[1:idxend,l],reverse(y[1:idxend,l]))
                Y = vcat(T[1:idxend,l],reverse(T[1:idxend,l]))
                line_vle = PlotlyBase.scatter(x=X,y=Y,mode="lines",line=PlotlyBase.attr(color=color, dash=style, width=3),name="VLE curve")
                PlotlyBase.addtraces!(plt,line_vle)
            end 
        end

    PlotlyBase.update!(plt,layout=PlotlyBase.Layout(xaxis = PlotlyBase.attr(range = [0.,1.])))
    if lle_present
        PlotlyBase.update!(plt,layout=PlotlyBase.Layout(yaxis = PlotlyBase.attr(range = [Tmin,Tmax*1.05])))
    else
        PlotlyBase.update!(plt,layout=PlotlyBase.Layout(yaxis = PlotlyBase.attr(range = [Tmin/1.05,Tmax*1.05])))
    end
    return plt
end
export Txy_diagram