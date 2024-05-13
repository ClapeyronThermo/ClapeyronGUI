# TO DOs for pxy_diagram
# 1. Add check for VLLE when it does not arise from an azeotrope

function ternary_diagrams(model,p,T,method=MichelsenTPFlash();Npoints=200,colors=[:red,:blue,:green,:purple,:black],styles=[:solid,:dash, :dot, :dashdot, :dashdotdot])

    # Basic settings for the plot
    plt = plot(grid=:off,framestyle=:box,foreground_color_legend = nothing,legend_font=font(12))

    # Check each side for a phase split
    z0 = ones(3,3)*0.5-eye(3)*(0.5-1e-10)
    for i in 1:3
        xtest = zeros(2,3)
        itest = 1
        idxtest = 1:2

        # perform stability analysis to verify if there is a phase split
        (xcand,G,phase0,phase1) = Clapeyron.tpd(model,p,T,z0[i,:])
        gL0 = gibbs_free_energy(model,p,T,z0[i,:],phase=:liquid)
        gV0 = gibbs_free_energy(model,p,T,z0[i,:],phase=:vapour)
        ntp = length(x)
        for j in 1:ntp
            if phase0[j]==:vapour && gV0<gL0
                xtest[itest] = xcand[j]
                idxtest[itest] = j
                itest+=1
            elseif phase0[j]==:liquid && gL0<gV0
                xtest[itest] = xcand[j]
                idxtest[itest] = j
                itest+=1
            end
        end
        if all(xtest[:,1].!=0)
            # We have a phase split on the i-th side
            K0 = xtest[1,:]./xtest[2,:]
            if phase1[idxtest]==[:vapour,:liquid] || phase1[idxtest]==[:liquid,:vapour]
                equilibrium=:vle
            else
                equilibrium=:lle
            end
            x_eq = zeros(Npoints,6)
            (x,n,G) = tp_flash(model,p,T,z0[i,:],method(K0=K0,equilibrium=equilibrium))
            x_eq[1,1:3] = x[1,:]
            x_eq[1,4:6] = x[2,:]
            z0[i,i] += 1/Npoints
            idxend = Npoints
            for k in 2:Npoints
                (x,n,G) = tp_flash(model,p,T,z0[i,:],method(K0=K0,equilibrium=equilibrium))
                x_eq[k,1:3] = x[1,:]
                x_eq[k,4:6] = x[2,:]
                dx1 = x_eq[k,1:3]-x_eq[k-1,1:3]
                dx2 = x_eq[k,4:6]-x_eq[k-1,4:6]
                dz0 = (dx1+dx2)/2
                z0[i,:] += dz0
                K0 = x[1,:]./x[2,:]
                if abs(x[1,1]-x[2,1])<1e-5
                    # if the two-phase region closes, break the loop
                    idxend = k-1
                    break
                end
            end
            x_eq = hcat(x_eq[1:idxend,1:3],x_eq[1:idxend,4:6])
        end
    end
    return plt
end
    
export pxy_diagram