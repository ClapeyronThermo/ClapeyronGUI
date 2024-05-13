# TO DOs for Txy_diagram
# 1. Add check for VLLE when it does not arise from an azeotrope
# 2. Add critical point evaluators
# 3. Add the case where both components are supercritical
function Txy_diagram(models,p;iscrit=nothing,check_lle=false,lle_present=false,Npoints=200,colors=[:red,:blue,:green,:purple,:black],styles=[:solid,:dash, :dot, :dashdot, :dashdotdot])
    nmodels = length(models)
    npres = length(p)

    # Basic settings for the plot
    plt = plot(grid=:off,framestyle=:box,foreground_color_legend = nothing,legend_font=font(12))

    Tmax = 0.
    Tmin = Inf
    T_lle = []
    Klle = []
    zlle = [0.5,0.5]

    for i in 1:nmodels
        for j in 1:npres
            if isnothing(iscrit)
                pures = split_model(models[i])
                crit = crit_pure.(pures)
                pc = [crit[k][2] for k in 1:length(pures)]
                crit = pc.<p[j]
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
                v0 = [Clapeyron.x0_bubble_temperature(models[i],p[j],[x[1],1-x[1]])]
                npaths = 1
            else
                sat = saturation_temperature.(pures,p[j])
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

                v0 = [Clapeyron.x0_bubble_temperature(models[i],p[j],[x[1,k],1-x[1,k]]) for k in 1:2]
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
                    bub = bubble_temperature(models[i],p[j],[x[k,l],1-x[k,l]];v0=v0[l])
                    v0[l] = [bub[1],log10(bub[2]),log10(bub[3]),bub[4][1]+dx,bub[4][2]-dx]
                    T[k,l] = bub[1]
                    y[k,l] = bub[4][1]

                    if k>1
                        # If there's a sudden change in direction along the x-y and p space, then there's VLLE
                        # This step can only help identify VLLE when it is also part of an azeotrope i.e. x < y < xx
                        # We need something to handle x < xx < y
                        if (y[k,l]-y[k-1,l])/dx < 0 && (T[k,l] - T[k-1,l])/(y[k,l]-y[k-1,l]) > 0
                            println("Vapour-Liquid-Liquid equilibrium is present")
                            x0 = [x[k,l],1-x[k,l]]
                            y0 = bub[4]
                            xx0 = 1 .- x0

                            vll0 = log10(volume(models[i],p[j],bub[1],xx0;phase=:l))
                            v0_dew = [bub[1], vll0, log10(bub[3]), xx0[1], xx0[2]]
                            dew = dew_temperature(models[i],p[j],y0;v0=v0_dew)
                            
                            T0 = (dew[1]+bub[1])/2

                            Klle = dew[4]./x0
                            lle = tp_flash(models[i],p[j],T0,(x0.+dew[4])./2,MichelsenTPFlash(K0=Klle,equilibrium=:lle))
                            Kvle = bub[4]./x0
                            vle = tp_flash(models[i],p[j],T0,(x0.+y0)./2,MichelsenTPFlash(K0=Kvle,equilibrium=:vle))

                            x0 = lle[1][1,1]
                            y0 = vle[1][2,1]
                            xx0 = lle[1][2,1]

                            vl0 = volume(models[i],p[j],T0,[x0,1-x0];phase=:l)
                            vll0 = volume(models[i],p[j],T0,[xx0,1-xx0];phase=:l)
                            vv0 = volume(models[i],p[j],T0,[y0,1-y0];phase=:v)

                            v0_vlle = [T0,log10(vl0),log10(vll0),log10(vv0),x0,xx0,y0]

                            vlle = VLLE_temperature(models[i],p[j];v0=v0_vlle)
                            T_vlle = vlle[1]
                            x_vlle = vlle[5][1]
                            xx_vlle = vlle[6][1]
                            y_vlle = vlle[7][1]

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
                            zlle = (xx_vlle+x_vlle)/2
                            zlle = [zlle, 1-zlle]
                            Klle = [xx_vlle/x_vlle,(1-xx_vlle)/(1-x_vlle)]
                            T_lle = LinRange(T_vlle,maximum((250,T_vlle-100)),Npoints)
                        end
                    end

                    # If the volumes of the two phases are similar or get reversed, we've reached a critical point
                    # We can find that critical point and break the while loop and go to the next branch
                    if bub[3]-bub[2]<1e-6
                        # x0 = (x[k-1,l]+y[k-1,l])/2
                        # v0_crit = [log10((bub[2]+bub[3])/2),[x0,1-x0]]
                        # crit_point = UCST_mix(models[i],T[j];v0=v0_crit)
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
                    T_test = maximum([Tmin-100,250])
                    x_test = [0.5,0.5]
                    tpd = Clapeyron.tpd(models[i],p[j],T_test,x_test)
                    if length(tpd[1])>1
                        Klle = tpd[1][1]./tpd[1][2]
                        lle_present=true
                        T_lle = LinRange(T_test,Tmin,Npoints)
                    end
                end

                # If LLE is present, use a flash algorithm to find the LLE curve
                if lle_present
                    println("Liquid-Liquid Equilibria is present")
                    idxend = Npoints
                    x_lle= zeros(Npoints)
                    xx_lle = zeros(Npoints)
                    for k in 1:Npoints
                        flash = tp_flash(models[i],p[j],T_lle[k],zlle,MichelsenTPFlash(K0=Klle,equilibrium=:lle))
                        x_lle[k] = flash[1][1,1]
                        xx_lle[k] = flash[1][2,1]
                        Klle = flash[1][1,:]./flash[1][2,:]
                        zlle = (flash[1][1,:]+flash[1][2,:])./2

                        # If the compositions of the two phase become similar, we have reached a critical point
                        if abs(x_lle[k]-xx_lle[k])<1e-4
                            idxend = k
                            # x0 = (x_lle[k,l]+xx_lle[k,l])/2
                            # vl0 = volume(models[i],p_lle[k],T[j],[x0,1-x0];phase=:l)
                            # v0_crit = [log10(vl0),[x0,1-x0]]
                            # crit_point = UCST_mix(models[i],T[j];v0=v0_crit)
                            # x_lle[k,l] = crit_point[3][1]
                            # xx_lle[k,l] = crit_point[3][1]
                            # p_lle[k,l] = crit_point[1]
                            break
                        end
                    end
                    plot!(plt,x_lle[1:idxend],T_lle[1:idxend],color=colors[j],linestyle=styles[i],line = (:path, 2),label = false)
                    plot!(plt,xx_lle[1:idxend],T_lle[1:idxend],color=colors[j],linestyle=styles[i],line = (:path, 2),label = false)
                    lle_present = false
                end

                if idxend == Npoints
                    X = vcat(x[:,l],reverse(y[:,l]))
                    Y = vcat(T[:,l],reverse(T[:,l]))
                    plot!(plt,X,Y,color=colors[j],linestyle=styles[i],line = (:path, 2),label = false)
                    break
                else
                    X = vcat(x[1:idxend,l],reverse(y[1:idxend,l]))
                    Y = vcat(T[1:idxend,l],reverse(T[1:idxend,l]))
                    plot!(plt,X,Y,color=colors[j],linestyle=styles[i],line = (:path, 2),label = false)
                end 
            end
        end
    end

    xlabel!(plt,"molar composition of "*models[1].components[1],xguidefontsize=12)
    ylabel!(plt,"Temperature / K",yguidefontsize=12)
    ylims!(plt,(Tmin/1.1,Tmax*1.1))
    xlims!(plt,(0,1))

    return plt
end
    
export Txy_diagram