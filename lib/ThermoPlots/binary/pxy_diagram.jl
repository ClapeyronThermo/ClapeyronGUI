# TO DOs for pxy_diagram
# 1. Add check for VLLE when it does not arise from an azeotrope
import PlotlyJS


function pxy_diagram(model,T;iscrit=nothing,check_lle=false,lle_present=false,Npoints=1000,color=:red,style=:solid)

    # Basic settings for the plot
    layout = PlotlyJS.Layout(autosize=false,width=700,height=470,
                    xaxis = PlotlyJS.attr(title = "Molar composition of "*model.components[1], font_size=12, showgrid=false, ticks="inside",mirror=true,showline=true,linecolor="black"),
                    yaxis = PlotlyJS.attr(title = "Pressure / bar", font_size=12, showgrid=false, ticks="inside",mirror=true,showline=true,linecolor="black"),
                    showlegend=false, plot_bgcolor="white")
    plt = PlotlyJS.plot(PlotlyJS.scatter(),layout)
    pmax = 0.
    pmin = Inf

    _pxy_diagram!(plt,model,T,pmax,pmin;iscrit=iscrit,check_lle=check_lle,lle_present=lle_present,Npoints=Npoints,color=color,style=style)
    return plt
end

function pxy_diagram!(plt,model,T;iscrit=nothing,check_lle=false,lle_present=false,Npoints=1000,color=:red,style=:solid)
    p = plt.plot.layout[:yaxis][:range]
    pmax = p[2]*1e5
    pmin = p[1]*1e5
    _pxy_diagram!(plt,model,T,pmax,pmin;iscrit=iscrit,check_lle=check_lle,lle_present=lle_present,Npoints=Npoints,color=color,style=style)
    return plt
end

function _pxy_diagram!(plt,model,T,pmax,pmin;iscrit=nothing,check_lle=false, check_vlle = false,lle_present=false,Npoints=200,color=:red,style=:solid)
    stop_vlle = false
    p_lle = []
    Klle = []
    zlle = [0.5,0.5]
    if isnothing(iscrit)
        pures = split_model(model)
        crit = crit_pure.(pures)
        Tc = [crit[k][1] for k in 1:length(pures)]
        crit = Tc.<T
    else
        crit = iscriT
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
        v0 = [Clapeyron.x0_bubble_pressure(model,T,[x[1],1-x[1]])]
        npaths = 1
    else
        sat = saturation_pressure.(pures,T)
        x = ones(Npoints,2)
        # When both components are subcritical, we need to scan in both directions
        # This is because there is a possibility that there are two critical points at a given temperature
        if sat[1][1]>sat[2][1]
            x[:,1] = LinRange(0.,1.,Npoints)
            x[:,2] = LinRange(1.,0.,Npoints)
        else
            x[:,1] = LinRange(1.,0.,Npoints)
            x[:,2] = LinRange(0.,1.,Npoints)
        end

        v0 = [Clapeyron.x0_bubble_pressure(model,T,[x[1,k],1-x[1,k]]) for k in 1:2]
        npaths = 2
    end


    # Scan through the VLE region
    y = zeros(Npoints,npaths)
    p = zeros(Npoints,npaths)
    for l in 1:npaths
        idxend = Npoints
        idx_vlle = [1,1]
        dx = x[2,l]-x[1,l]
        k=1
        n=1
        while k <= Npoints
            # println(k)

            # Obtain the bubble points
            try
                bub = bubble_pressure(model,T,[x[k,l],1-x[k,l]];v0=v0[l])
                v0[l] = [log10(bub[2]),log10(bub[3]),bub[4][1]+dx,bub[4][2]-dx]
                p[k,l] = bub[1]
                y[k,l] = bub[4][1]
            

                if k>2
                    # If there's a sudden change in direction along the x-y and p space, then there's VLLE
                    # This step can only help identify VLLE when it is also part of an azeotrope i.e. x < y < xx
                    # We need something to handle x < xx < y
                    if !check_vlle && (y[k,l]-y[k-1,l])*(y[k-1,l]-y[k-2,l]) < 0 && !stop_vlle                        
                        idx_vlle[n] = k-1
                        if n==2
                            check_vlle = true
                            n=1
                        else
                            println("Vapour-Liquid-Liquid equilibrium may be present")
                            n+=1
                        end
                        # break
                    elseif check_vlle && y[k]>y[idx_vlle[1]] && !stop_vlle
                        println("Vapour-Liquid-Liquid equilibrium is present")

                        idx = argmin(abs.(y[1:idx_vlle[1],l].-y[idx_vlle[2]:k-1,l]'))
                        y0 = y[idx[1],l]
                        x0 = x[idx[1],l]
                        xx0 = x[idx_vlle[2]+idx[2],l]
                        p0 = p[idx[1],l]

                        vll0 = log10(Clapeyron.volume(model,p0,T,[xx0,1-xx0];phase=:l))
                        vl0 = log10(Clapeyron.volume(model,p0,T,[x0,1-x0];phase=:l))
                        vv0 = log10(Clapeyron.volume(model,p0,T,[y0,1-y0];phase=:v))

                        v0_vlle = [vl0,vll0,vv0,x0,xx0,y0]
                        vlle = VLLE_pressure(model,T;v0=v0_vlle)
                        p_vlle = vlle[1]
                        x_vlle = vlle[5][1]
                        xx_vlle = vlle[6][1]
                        y_vlle = vlle[7][1]

                        line_vlle = PlotlyJS.scatter(x=sort([xx_vlle,x_vlle,y_vlle]),y=[p_vlle,p_vlle,p_vlle]./1e5,mode="lines",line=PlotlyJS.attr(color=color, dash=style, width=3),name="VLLE curve")
                        PlotlyJS.addtraces!(plt,line_vlle)

                        idx_vlle = sum(p[1:idx_vlle[1],l].<p_vlle)
                        p[idx_vlle+1,l] = p_vlle
                        x[idx_vlle+1,l] = x_vlle
                        y[idx_vlle+1,l] = y_vlle

                        p[idx_vlle+2,l] = p_vlle
                        x[idx_vlle+2,l] = xx_vlle
                        y[idx_vlle+2,l] = y_vlle

                        p[idx_vlle+3:end,l] .= 0
                        x[idx_vlle+3:end,l] .= LinRange(xx_vlle+dx,x[end,l],Npoints-idx_vlle-2)
                        y[idx_vlle+3:end,l] .= 0
                        k = idx_vlle+2

                        v0[l] = [log10(vlle[3]),log10(vlle[4]),vlle[7][1]+dx,vlle[7][2]-dx]

                        # Given there is VLLE, there must also be LLE
                        lle_present = true
                        zlle = (xx_vlle+x_vlle)/2
                        zlle = [zlle, 1-zlle]
                        Klle = [xx_vlle/x_vlle,(1-xx_vlle)/(1-x_vlle)]
                        p_lle = LinRange(p_vlle,1.5*p_vlle,Npoints)
                        stop_vlle = true

                        # idxend = idx[1]
                        # break
                    end
                end

                # If the Clapeyron.volumes of the two phases are similar or get reversed, we've reached a critical point
                # We can find that critical point and break the while loop and go to the next branch
                if bub[3]-bub[2]<1e-6
                    x0 = (x[k-1,l]+y[k-1,l])/2
                    v0_crit = [log10((bub[2]+bub[3])/2),[x0,1-x0]]
                    crit_point = UCST_mix(model,T;v0=v0_crit)
                    x[k,l] = crit_point[3][1]
                    y[k,l] = crit_point[3][1]
                    p[k,l] = crit_point[1]
                    idxend = k
                    break
                end
                k+=1
            catch
                if check_vlle 
                    stop_vlle = true
                end
                dx = dx/2
                x[k,l] = x[k-1,l]+dx
                v0[l][3] = v0[l][3]-dx
                v0[l][4] = v0[l][4]+dx
                continue
            end
        end

        # Find the bounds of the plot
        pmax = maximum(vcat(pmax,p[:,l]))
        pmin = minimum(vcat(pmin,p[:,l]))

        # If we need to check for LLE, we need to look at higher pressures.
        # For visual clarity, we only really need to consider LLE when it takes up half the phase diagram
        if check_lle
            p_test = 10*pmax
            x_test = [0.5,0.5]
            tpd = Clapeyron.tpd(model,p_test,T,x_test)
            if length(tpd[1])>1
                Klle = tpd[1][1]./tpd[1][2]
                lle_present=true
                p_lle = LinRange(p_test,pmax,Npoints)
                p_lle = [p_lle[k] for k in 1:Npoints]
            end
        end

        # If LLE is present, use a flash algorithm to find the LLE curve
        if lle_present
            println("Liquid-Liquid Equilibria is present")
            idxend = Npoints
            x_lle= zeros(Npoints)
            xx_lle = zeros(Npoints)
            for k in 1:Npoints
                flash = tp_flash(model,p_lle[k],T,zlle,MichelsenTPFlash(K0=Klle,equilibrium=:lle))
                x_lle[k] = flash[1][1,1]
                xx_lle[k] = flash[1][2,1]
                Klle = flash[1][1,:]./flash[1][2,:]
                zlle = (flash[1][1,:]+flash[1][2,:])./2
                # If the compositions of the two phase become similar, we have reached a critical point
                if abs(x_lle[k]-xx_lle[k])<1e-4
                    idxend = k
                    x0 = (x_lle[k,l]+xx_lle[k,l])/2
                    vl0 = Clapeyron.volume(model,p_lle[k],T,[x0,1-x0];phase=:l)
                    v0_crit = [log10(vl0),[x0,1-x0]]
                    crit_point = UCST_mix(model,T;v0=v0_crit)
                    x_lle[k] = crit_point[3][1]
                    xx_lle[k] = crit_point[3][1]
                    p_lle[k] = crit_point[1]
                    break
                end
            end
            line_lle1 = PlotlyJS.scatter(x=x_lle[1:idxend],y=p_lle[1:idxend]./1e5,mode="lines",line=PlotlyJS.attr(color=color, dash=style, width=3),name="LLE curve")
            line_lle2 = PlotlyJS.scatter(x=xx_lle[1:idxend],y=p_lle[1:idxend]./1e5,mode="lines",line=PlotlyJS.attr(color=color, dash=style, width=3),name="LLE curve")
            PlotlyJS.addtraces!(plt,line_lle1,line_lle2)

            lle_present = false
        end

        if idxend == Npoints
            X = vcat(x[:,l],reverse(y[:,l]))
            Y = vcat(p[:,l],reverse(p[:,l]))
            line_vle = PlotlyJS.scatter(x=X,y=Y./1e5,mode="lines",line=PlotlyJS.attr(color=color, dash=style, width=3),name="VLE curve")
            PlotlyJS.addtraces!(plt,line_vle)
            break
        else
            X = vcat(x[1:idxend,l],reverse(y[1:idxend,l]))
            Y = vcat(p[1:idxend,l],reverse(p[1:idxend,l]))
            line_vle = PlotlyJS.scatter(x=X,y=Y./1e5,mode="lines",line=PlotlyJS.attr(color=color, dash=style, width=3),name="VLE curve")
            PlotlyJS.addtraces!(plt,line_vle)
        end 
    end
    
    PlotlyJS.update!(plt,layout=PlotlyJS.Layout(yaxis = PlotlyJS.attr(range = [pmin/1.1/1e5,pmax*1.1/1e5])))
    PlotlyJS.update!(plt,layout=PlotlyJS.Layout(xaxis = PlotlyJS.attr(range = [0.,1.])))
    return plt
end
    
export pxy_diagram, pxy_diagram!