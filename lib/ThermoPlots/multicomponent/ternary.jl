norm(z) = sqrt(sum(z.^2))


function make_ax(title, tickangle)
    PlotlyBase.attr(title=title, titlefont_size=20, tickangle=tickangle,
        tickfont_size=15, tickcolor="rgb(0, 0, 0)", ticklen=5,gridcolor="rgba(0, 0, 0)",linecolor="rgba(0, 0, 0)",
        showline=true, showgrid=true)
end

function ternary_diagram(model,p,T; Npoints=200,color=:red,style=:solid,check_three_phase=false)
    components = model.components
    layout = PlotlyBase.Layout(
        ternary=PlotlyBase.attr(
            sum=1,
            aaxis=make_ax(components[1], 0),
            baxis=make_ax(components[2], 45),
            caxis=make_ax(components[3], -45),
            bgcolor="#ffffff",
        ))
    plt = plot(PlotlyBase.scatterternary(),layout)
    if typeof(model)<:Clapeyron.ActivityModel
        method = RRTPFlash
    else
        method = MichelsenTPFlash
    end

    return _ternary_diagram(plt,model,p,T,method;Npoints=Npoints,color=color,style=style,check_three_phase=check_three_phase)
end

function _ternary_diagram(plt,model,p,T,method=MichelsenTPFlash;Npoints=200,color=:red,style=:solid,check_three_phase=false)

    # Basic settings for the plot
    pures = split_model(model)
    vol = Clapeyron.volume.(pures,p,T)
    Π = Clapeyron.pip.(pures,vol,T)
    phases = [if Π[i] > 1. "liquid" else "vapour" end for i in 1:3]
    if all(phases .== "vapour")
        println("All components are vapour. There will be no phase splits.")
    end

    # Check each side for a phase split
    is_closed = true
    for i in 1:3 # Go through each side
        xtest = ones(3)
        xtest[i] = 0

        if phases[1:3 .!= i] == ["liquid","liquid"] || typeof(model) <: Clapeyron.ActivityModel
            equilibrium = :lle
        else
            equilibrium = :vle
        end

        model_bin = index_reduction(model,xtest)[1]
        (x,n,G) = tp_flash(model_bin,p,T,[0.5,0.5],method(equilibrium=equilibrium))
        # println(x)

        if abs(x[1,1]-x[2,1])./x[1,1] .> 1e-5
            # println("Phase split detected on side $model.components[$i] with equilibrium $equilibrium")
            (x,n,G) = tp_flash(model_bin,p,T,[0.5,0.5],method(equilibrium=equilibrium))
            K0_bin = x[1,:]./x[2,:]
            K0 = zeros(3)
            K0[1:3 .!==i] = K0_bin
            K0[i] = 1.

            x_eq = zeros(Npoints,6)
            z0 = ones(3)*1e-5
            z0[1:3 .!==i] .= 0.5
            z0 = z0/sum(z0)
            (x,n,G) = tp_flash(model,p,T,z0,method(equilibrium=equilibrium))
            x_eq[1,1:3] = x[1,:]
            x_eq[1,4:6] = x[2,:]
            # println(x_eq[1,:])
            idxend = Npoints
            z0[i] += 1.5/Npoints
            z0 = z0/sum(z0)
            if equilibrium == :lle
                K0 = x[1,:]./x[2,:]
            else
                K0 = nothing
            end

            for i in 2:Npoints
                (x,n,G) = tp_flash(model,p,T,z0,method(K0=K0,equilibrium=equilibrium))
                # println(x)
                x_eq[i,1:3] = x[1,:]
                x_eq[i,4:6] = x[2,:]
                # println(x_eq[i,:])
                dx1 = x_eq[i,1:3]-x_eq[i-1,1:3]
                dx2 = x_eq[i,4:6]-x_eq[i-1,4:6]
                dz0 = (dx1+dx2)/2
                dz0 = dz0/norm(dz0)*1.5/Npoints
                if norm(dz0) > 1.5/Npoints
                    # println(norm(dz0))
                    dz0 = dz0/norm(dz0)*1.5/Npoints
                    # println("big")
                elseif norm(dz0) < 1e-3
                    dz0 = dz0/norm(dz0)*1.5/Npoints
                    # println("small")
                end
                z0 += dz0
                # if sqrt(sum(dz0.^2)) < 1/Npoints
                #     dz0 = dz0/sqrt(sum(dz0.^2))*1/Npoints
                # end
                z0 = z0/sum(z0)
                if equilibrium == :lle
                    K0 = x[1,:]./x[2,:]
                else
                    K0 = nothing
                end
                if maximum(abs.(x[1,:].-x[2,:]))<1e-5
                    # if the two-phase region closes, break the loop
                    # println("met tolerance")
                    is_closed = true
                    idxend = i-1
                    break
                elseif any(x.==0)
                    # if the two-phase region closes, break the loop
                    # println("met zero")
                    is_closed = false
                    idxend = i-1
                    break
                end
            end
            
            # println(idxend)
            if !is_closed
                line1 = PlotlyBase.scatterternary( mode="lines",
                                        a=x_eq[1:idxend,1],
                                        b=x_eq[1:idxend,2],
                                        c=x_eq[1:idxend,3],
                                        line=PlotlyBase.attr(color=color, dash=style, width=3),
                                        name="")
                line2 = PlotlyBase.scatterternary( mode="lines",
                                        a=x_eq[1:idxend,4],
                                        b=x_eq[1:idxend,5],
                                        c=x_eq[1:idxend,6],
                                        line=PlotlyBase.attr(color=color, dash=style, width=3),
                                        name="")  
                PlotlyBase.addtraces!(plt,line1,line2)  
                break
            else
                X = vcat(x_eq[1:idxend,1],reverse(x_eq[1:idxend,4]))
                Y = vcat(x_eq[1:idxend,2],reverse(x_eq[1:idxend,5]))
                Z = vcat(x_eq[1:idxend,3],reverse(x_eq[1:idxend,6]))
                line = PlotlyBase.scatterternary( mode="lines",
                                        a=X,
                                        b=Y,
                                        c=Z,
                                        line=PlotlyBase.attr(color=color, dash=style, width=3),
                                        name="")
                PlotlyBase.addtraces!(plt,line)
            end
        end
    end
    return plt
end

export ternary_diagram