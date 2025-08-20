module POLYMER
using GenieFramework
using Clapeyron, Main.ThermoPlots
using CoolProp
import PlotlyBase, PlotlyKaleido
@genietools
import Main.@timeout

epsilon = Dict("PBD" => 288.84, 
               "PS" => 348.2, 
               "PDMS" => 204.9,
               "LDPE" => 249.5, 
               "PBMA" => 264.7, 
               "PIB" => 267.6, 
               "PαMS" => 354.05, 
               "PMMA" => 264.60)
sigma = Dict("PBD" => 4.097, 
             "PS" => 4.152, 
             "PDMS" => 3.531,
             "LDPE" => 4.022, 
             "PBMA" => 3.884, 
             "PIB" => 4.117, 
             "PαMS" => 4.204, 
             "PMMA" => 3.553)

segment = Dict("PBD" => 0.0245, 
               "PS" => 0.0205, 
               "PDMS" => 0.0324,
               "LDPE" => 0.0263, 
               "PBMA" => 0.0241, 
               "PIB" => 0.0233, 
               "PαMS" => 0.0204, 
               "PMMA" => 0.027)

SAFTgammaMieGroups = Dict(
    "PBD" => ["CH="=>2,"CH2"=>1],
    "PS" => ["CH2"=>1,"aCH"=>5,"aCCH"=>1],
    "LDPE" => ["CH2"=>2],
    "PBMA" => ["CH3"=>2,"CH2"=>4,"C"=>1,"COO"=>1],
    "PIB" => ["CH2"=>1,"CH3"=>2,"C"=>1],
    "PMMA" => ["CH2"=>1,"CH3"=>2,"COO"=>1,"C"=>1]
)

UNIFACGroups = Dict(
    "PBD" => ["CH=CH"=>1,"CH2"=>1],
    "PS" => ["CH2"=>1,"ACH"=>5,"ACCH"=>1],
    "LDPE" => ["CH2"=>2],
    "PBMA" => ["CH3"=>2,"CH2"=>4,"C"=>1,"COO"=>1],
    "PIB" => ["CH2"=>1,"CH3"=>2,"C"=>1],
    "PMMA" => ["CH2"=>1,"CH3"=>2,"COO"=>1,"C"=>1]
)

Monomer_MW = Dict(
    "PBD" => 54.09044,
    "PS" => 104.14912,
    "PDMS" => 74.0,
    "LDPE" => 28.05316,
    "PBMA" => 142.19556,
    "PIB" => 56.10632,
    "PαMS" => 118.1757,
    "PMMA" => 100.11582
)

@app begin
    @in tab_selected = "gas_solubility"
    @in Select_eos_1 = "PCSAFT"
    @out Select_eos_list_1 = ["PCSAFT","SAFTγMie"]
    @in Select_eos_2 = "PCSAFT"
    @out Select_eos_list_2 = ["PCSAFT","SAFTγMie","UNIFAC","FH"]
    @in polymer = "PBD"
    @in Mw_poly = 1000
    @in polymer1 = "PBD"
    @in Mw_poly1 = 2000
    @in polymer2 = "PS"
    @in Mw_poly2 = 4000
    @in gas = "carbon dioxide"

    @in temp = 298.15

    @in new_T_button = false
    @in new_p_button = false

    @out polymer_list = ["PBD","PS","PDMS","LDPE","PBMA","PIB","PαMS","PMMA"]
    @out gas_list = ["carbon dioxide","nitrogen","oxygen","methane","ethane"]

    @out trace_T = []
    @out layout_T = PlotlyBase.Layout(autosize=false,width=700,height=470,
    yaxis = PlotlyBase.attr(title = "Pressure  / bar", font_size=12, showgrid=false,            
                                      ticks="inside",mirror=true,showline=true,linecolor="black"),
                         xaxis = PlotlyBase.attr(title = "mass fraction / (kg/kg)", font_size=12, showgrid=false,       
                                      ticks="inside",mirror=true,showline=true,linecolor="black"),
                         showlegend=false, plot_bgcolor="white")

    @out trace_p = []
    @out layout_p = PlotlyBase.Layout(autosize=false,width=700,height=470,
    xaxis = PlotlyBase.attr(title = "mass fraction  / (kg/kg)", font_size=12, showgrid=false,            
                                    ticks="inside",mirror=true,showline=true,linecolor="black"),
                        yaxis = PlotlyBase.attr(title = "Temperature / K", font_size=12, showgrid=false,       
                                    ticks="inside",mirror=true,showline=true,linecolor="black"),
                        showlegend=false, plot_bgcolor="white")
    
    @onbutton new_T_button begin
        Npoints = 200
        # try
        if Select_eos_1 == "PCSAFT"
            polymer_param = ParamTable(:like,(
                species = [polymer],
                Mw = [Mw_poly],
                epsilon = [epsilon[polymer]],
                sigma = [sigma[polymer]],
                segment = [segment[polymer]*Mw_poly],
                n_H = [0],
                n_e = [0]))
            model = PCSAFT([polymer,gas]; userlocations=[polymer_param])
        elseif Select_eos_1 == "SAFTγMie"
            N = round(Mw_poly / Monomer_MW[polymer])
            monomer_groups = SAFTgammaMieGroups[polymer]
            poly_groups = [(monomer_groups[i][1]=>Int(monomer_groups[i][2]*N)) for i in 1:length(monomer_groups)]
            model = SAFTγMie([(polymer,poly_groups),gas])
        else
            notify(__model__, "Equation of State $Select_eos_1 is not implemented.", :warning)
            throw(TypeError("Equation of State $Select_eos_1 is not implemented."))
        end
        # catch
        #     notify(__model__, "Species $gas is not available in $Select_eos_1.", :warning)
        #     throw(TypeError("Species $gas is not available in $Select_eos_1."))
        # end
        Mw = Clapeyron.mw(model)
        w = LinRange(1,0.8,Npoints)
        x = @. w*Mw[1]/(w*Mw[1]+(1-w)*Mw[2])
        
        method = FugBubblePressure(nonvolatiles=[polymer])
        p = zeros(Npoints)

        for i in 2:Npoints
            bub = bubble_pressure(model,temp,[x[i],1-x[i]],method)
            p[i] = bub[1]
        end

        trace_T = [PlotlyBase.scatter(x=1 .-w,y=p./1e5,mode="lines",line=PlotlyBase.attr(color="red", dash="solid", width=3),name="")]
        layout_T = PlotlyBase.Layout(autosize=false,width=700,height=470,
        yaxis = PlotlyBase.attr(title = "Pressure / bar", font_size=12, showgrid=false,            
                                          ticks="inside",mirror=true,showline=true,linecolor="black"),
                             xaxis = PlotlyBase.attr(title = gas*" mass fraction / (kg/kg)", font_size=12, showgrid=false,       
                                          ticks="inside",mirror=true,showline=true,linecolor="black"),
                             showlegend=false, plot_bgcolor="white")
    end

    @onbutton new_p_button begin
        Npoints = 1000
        if Select_eos_2 == "PCSAFT"
            polymer_param = ParamTable(:like,(
                species = [polymer1,polymer2],
                Mw = [Mw_poly1,Mw_poly2],
                epsilon = [epsilon[polymer1],epsilon[polymer2]],
                sigma = [sigma[polymer1],sigma[polymer2]],
                segment = [segment[polymer1]*Mw_poly1,segment[polymer2]*Mw_poly2],
                n_H = [0,0],
                n_e = [0,0]))
            model = PCSAFT([polymer1,polymer2]; userlocations=[polymer_param])
        elseif Select_eos_2 == "SAFTγMie"
            N1 = round(Mw_poly1 / Monomer_MW[polymer1])
            N2 = round(Mw_poly2 / Monomer_MW[polymer2])
            monomer_groups1 = SAFTgammaMieGroups[polymer1]
            monomer_groups2 = SAFTgammaMieGroups[polymer2]
            poly_groups1 = [(monomer_groups1[i][1]=>Int(monomer_groups1[i][2]*N1)) for i in 1:length(monomer_groups1)]
            poly_groups2 = [(monomer_groups2[i][1]=>Int(monomer_groups2[i][2]*N2)) for i in 1:length(monomer_groups2)]
            model = SAFTγMie([(polymer1,poly_groups1),(polymer2,poly_groups2)])
        elseif Select_eos_2 == "FH"
            N1 = round(Mw_poly1 / Monomer_MW[polymer1])
            N2 = round(Mw_poly2 / Monomer_MW[polymer2])
            model = FloryHuggins([polymer1,polymer2],[N1,N2])
        elseif Select_eos_2 == "UNIFAC"
            N1 = round(Mw_poly1 / Monomer_MW[polymer1])
            N2 = round(Mw_poly2 / Monomer_MW[polymer2])
            monomer_groups1 = UNIFACGroups[polymer1]
            monomer_groups2 = UNIFACGroups[polymer2]
            poly_groups1 = [(monomer_groups1[i][1]=>Int(monomer_groups1[i][2]*N1)) for i in 1:length(monomer_groups1)]
            poly_groups2 = [(monomer_groups2[i][1]=>Int(monomer_groups2[i][2]*N2)) for i in 1:length(monomer_groups2)]
            model = UNIFAC([(polymer1,poly_groups1),(polymer2,poly_groups2)]; puremodel=BasicIdeal)
        elseif Select_eos_2 == "UNIFACFV"
            N1 = round(Mw_poly1 / Monomer_MW[polymer1])
            N2 = round(Mw_poly2 / Monomer_MW[polymer2])
            monomer_groups1 = UNIFACGroups[polymer1]
            monomer_groups2 = UNIFACGroups[polymer2]
            poly_groups1 = [(monomer_groups1[i][1]=>Int(monomer_groups1[i][2]*N1)) for i in 1:length(monomer_groups1)]
            poly_groups2 = [(monomer_groups2[i][1]=>Int(monomer_groups2[i][2]*N2)) for i in 1:length(monomer_groups2)]
            model = UNIFACFV([(polymer1,poly_groups1),(polymer2,poly_groups2)]; puremodel=BasicIdeal)
        end
        Mw = [Mw_poly1,Mw_poly2]

        T = LinRange(300,500,Npoints)
        p = 1e5 
        w = zeros(Npoints,2)
        w0 = [0.5,0.5]
        z0 = w0./Mw ./ sum(w0./Mw)
        K0 = [1e3,1e-3]
        idxend = Npoints

        (x,n,G) = tp_flash(model,p,T[1],z0,RRTPFlash(equilibrium=:lle))
        # println(x)
        K0 = nothing
        if abs(x[1,1]-x[2,1])/x[1,1] > 1e-4
            if x[1,1] < x[2,1]
                w[1,1] = x[2,1].*Mw[1] ./ sum(x[2,:].*Mw)
                w[1,2] = x[1,1].*Mw[1] ./ sum(x[1,:].*Mw)
                K0 = x[1,:]./x[2,:]
            else
                w[1,1] = x[1,1].*Mw[2] ./ sum(x[1,:].*Mw)
                w[1,2] = x[2,1].*Mw[2] ./ sum(x[2,:].*Mw)
                K0 = x[2,:]./x[1,:]
            end
            z0 = (x[1,:]+x[2,:])/2
            for i in 2:Npoints
                (x,n,G) = tp_flash(model,p,T[i],z0,RRTPFlash(K0=K0,equilibrium=:lle))
                if x[1,1] > x[2,1]
                    w[i,1] = x[1,1].*Mw[1] ./ sum(x[1,:].*Mw)
                    w[i,2] = x[2,1].*Mw[1] ./ sum(x[2,:].*Mw)
                    K0 = x[2,:]./x[1,:]
                else
                    w[i,1] = x[2,1].*Mw[1] ./ sum(x[2,:].*Mw)
                    w[i,2] = x[1,1].*Mw[1] ./ sum(x[1,:].*Mw)
                    K0 = x[1,:]./x[2,:]
                end
                z0 = (x[1,:]+x[2,:])/2
                if abs(x[1,1]-x[2,1])/x[1,1] < 1e-4
                    # println(x)
                    idxend = i
                    break
                end
            end
            if idxend == Npoints
                trace_p = [PlotlyBase.scatter(x=w[:,1],y=T,mode="lines",line=PlotlyBase.attr(color="green", dash="solid", width=3),name=""),
                           PlotlyBase.scatter(x=w[:,2],y=T,mode="lines",line=PlotlyBase.attr(color="green", dash="solid", width=3),name="")]
            else
                w = vcat(w[1:idxend,1],reverse(w[1:idxend,2]))
                T = vcat(T[1:idxend],reverse(T[1:idxend]))
                trace_p = [PlotlyBase.scatter(x=w,y=T,mode="lines",line=PlotlyBase.attr(color="green", dash="solid", width=3),name="")]
            end
            layout_p = PlotlyBase.Layout(autosize=false,width=700,height=470,
            yaxis = PlotlyBase.attr(title = "Temperature / K", font_size=12, showgrid=false,            
                                            ticks="inside",mirror=true,showline=true,linecolor="black",range=[300,500]),
                                xaxis = PlotlyBase.attr(title = polymer1*" mass fraction / (kg/kg)", font_size=12, showgrid=false,       
                                            ticks="inside",mirror=true,showline=true,linecolor="black",range=[0,1]),
                                showlegend=false, plot_bgcolor="white")

        end



        
        
    end
end
end