module ELECTROLYTES
using GenieFramework, StippleUI
using Clapeyron, Main.ThermoPlots
using CoolProp
import PlotlyBase, PlotlyKaleido
import Main.@timeout

@app begin
    @in tab_selected = "miac"
    @in cation = "sodium"
    @in anion = "chloride"
    @in temp = 298.15
    @in pre = 1e5
    @in Select_eos = "ePCSAFT"
    @in new_miac_button = false
    @in new_oc_button = false
    @in new_sat_button = false
    @out Select_eos_list = ["ePCSAFT","SAFTVREMie","eSAFTVRMie"]
    @out Select_cation_list = ["sodium", "lithium", "potassium","magnesium","calcium"] 
    @out Select_anion_list = ["chloride", "fluoride", "bromide","iodide"]
    @out trace_T = []
    @out layout_T = PlotlyBase.Layout(autosize=false,width=700,height=470,
    xaxis = PlotlyBase.attr(title = "sodium chloride molality / (mol/kg)", font_size=12, showgrid=false,            
                                      ticks="inside",mirror=true,showline=true,linecolor="black"),
                         yaxis = PlotlyBase.attr(title = "Mean Ionic Activity Coefficient", font_size=12, showgrid=false,       
                                      ticks="inside",mirror=true,showline=true,linecolor="black"),
                         showlegend=false, plot_bgcolor="white")
    @out trace_p = []
    @out layout_p = PlotlyBase.Layout(autosize=false,width=700,height=470,
    xaxis = PlotlyBase.attr(title = "sodium chloride molality / (mol/kg)", font_size=12, showgrid=false,            
                                      ticks="inside",mirror=true,showline=true,linecolor="black"),
                         yaxis = PlotlyBase.attr(title = "Osmotic Coefficient", font_size=12, showgrid=false,       
                                      ticks="inside",mirror=true,showline=true,linecolor="black"),
                         showlegend=false, plot_bgcolor="white")
    @out trace_pT = []
    @out layout_pT = PlotlyBase.Layout(autosize=false,width=700,height=470,
    xaxis = PlotlyBase.attr(title = "sodium chloride molality / (mol/kg)", font_size=12, showgrid=false,            
                                      ticks="inside",mirror=true,showline=true,linecolor="black"),
                         yaxis = PlotlyBase.attr(title = "Saturation Pressure / bar", font_size=12, showgrid=false,       
                                      ticks="inside",mirror=true,showline=true,linecolor="black"),
                         showlegend=false, plot_bgcolor="white")
    @onbutton new_miac_button begin
        N = 200
        eos = Symbol(Select_eos)
        if Select_eos == "ePCSAFT"
            try
                model = @eval $eos(["water08"],[$cation,$anion])
            catch
                notify(__model__, "Species $cation or $anion are not available in $Select_eos.", :warning)
                throw(TypeError("Species $cation or $anion are not available in $Select_eos."))
            end
            model = @eval $eos(["water08"],[$cation,$anion])
        else
            try
                model = @eval $eos(["water"],[$cation,$anion])
            catch
                notify(__model__, "Species $cation or $anion are not available in $Select_eos.", :warning)
                throw(TypeError("Species $cation or $anion are not available in $Select_eos."))
            end

            model = @eval $eos(["water"],[$cation,$anion])
        end

        if temp < 0
            notify(__model__, "Temperature must be greater than 273.15 K.", :warning)
            throw(TypeError("Temperature must be greater than 273.15 K."))
        end


        charges = model.charge
        salts = [(cation*" "*anion,(cation=>abs(charges[3]),anion=>abs(charges[2])))]

        m = LinRange(1e-5,6,N)
        # z = molality_to_composition.(model,salts,m)
        γ = zeros(N)
        for i in 1:N
            γ[i] = mean_ionic_activity_coefficient(model,salts,pre,temp,m[i])[1]
        end

        trace_T = [PlotlyBase.scatter(x=m,y=γ,mode="lines",line=PlotlyBase.attr(color="green", dash="solid", width=3),name="")]
        layout_T = PlotlyBase.Layout(autosize=false,width=700,height=470,
    xaxis = PlotlyBase.attr(title = "$cation $anion molality / (mol/kg)", font_size=12, showgrid=false,            
                                      ticks="inside",mirror=true,showline=true,linecolor="black",range = [0.,6.]),
                         yaxis = PlotlyBase.attr(title = "Mean Ionic Activity Coefficient", font_size=12, showgrid=false,       
                                      ticks="inside",mirror=true,showline=true,linecolor="black"),
                         showlegend=false, plot_bgcolor="white")
    end

    @onbutton new_oc_button begin
        N = 200
        eos = Symbol(Select_eos)
        if Select_eos == "ePCSAFT"
            try
                model = @eval $eos(["water08"],[$cation,$anion])
            catch
                notify(__model__, "Species $cation or $anion are not available in $Select_eos.", :warning)
                throw(TypeError("Species $cation or $anion are not available in $Select_eos."))
            end
            model = @eval $eos(["water08"],[$cation,$anion])
        else
            try
                model = @eval $eos(["water"],[$cation,$anion])
            catch
                notify(__model__, "Species $cation or $anion are not available in $Select_eos.", :warning)
                throw(TypeError("Species $cation or $anion are not available in $Select_eos."))
            end
            
            model = @eval $eos(["water"],[$cation,$anion])
        end

        if temp < 0
            notify(__model__, "Temperature must be greater than 273.15 K.", :warning)
            throw(TypeError("Temperature must be greater than 273.15 K."))
        end

        charges = model.charge
        salts = [(cation*" "*anion,(cation=>abs(charges[3]),anion=>abs(charges[2])))]

        m = LinRange(1e-5,6,N)
        # z = molality_to_composition.(model,salts,m)
        ϕ = zeros(N)
        for i in 1:N
            ϕ[i] = osmotic_coefficient(model,salts,pre,temp,m[i])[1]
        end

        trace_p = [PlotlyBase.scatter(x=m,y=ϕ,mode="lines",line=PlotlyBase.attr(color="green", dash="solid", width=3),name="")]
        layout_p = PlotlyBase.Layout(autosize=false,width=700,height=470,
    xaxis = PlotlyBase.attr(title = "$cation $anion molality / (mol/kg)", font_size=12, showgrid=false,            
                                      ticks="inside",mirror=true,showline=true,linecolor="black",range = [0.,6.]),
                         yaxis = PlotlyBase.attr(title = "Osmotic Coefficient", font_size=12, showgrid=false,       
                                      ticks="inside",mirror=true,showline=true,linecolor="black"),
                         showlegend=false, plot_bgcolor="white")
    end

    @onbutton new_sat_button begin
        N = 200
        eos = Symbol(Select_eos)
        if Select_eos == "ePCSAFT"
            try
                model = @eval $eos(["water08"],[$cation,$anion])
            catch
                notify(__model__, "Species $cation or $anion are not available in $Select_eos.", :warning)
                throw(TypeError("Species $cation or $anion are not available in $Select_eos."))
            end
            model = @eval $eos(["water08"],[$cation,$anion])
        else
            try
                model = @eval $eos(["water"],[$cation,$anion])
            catch
                notify(__model__, "Species $cation or $anion are not available in $Select_eos.", :warning)
                throw(TypeError("Species $cation or $anion are not available in $Select_eos."))
            end
            
            model = @eval $eos(["water"],[$cation,$anion])
        end

        if temp < 0
            notify(__model__, "Temperature must be greater than 273.15 K.", :warning)
            throw(TypeError("Temperature must be greater than 273.15 K."))
        end

        charges = model.charge
        salts = [(cation*" "*anion,(cation=>abs(charges[3]),anion=>abs(charges[2])))]
        nonvolatiles = [cation,anion]
        method = FugBubblePressure(nonvolatiles=nonvolatiles)

        m = LinRange(1e-5,6,N)
        z = molality_to_composition.(model,Ref(salts),m)
        p = zeros(N)
        for i in 1:N
            bub = bubble_pressure(model,temp,z[i],method)
            p[i] = bub[1]
        end

        trace_pT = [PlotlyBase.scatter(x=m,y=p./1e5,mode="lines",line=PlotlyBase.attr(color="green", dash="solid", width=3),name="")]
        layout_pT = PlotlyBase.Layout(autosize=false,width=700,height=470,
    xaxis = PlotlyBase.attr(title = "$cation $anion molality / (mol/kg)", font_size=12, showgrid=false,            
                                      ticks="inside",mirror=true,showline=true,linecolor="black",range = [0.,6.]),
                         yaxis = PlotlyBase.attr(title = "Saturation Pressure / bar", font_size=12, showgrid=false,       
                                      ticks="inside",mirror=true,showline=true,linecolor="black"),
                         showlegend=false, plot_bgcolor="white")
    end

end
end