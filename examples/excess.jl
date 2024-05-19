module EXCESS
using GenieFramework
using Clapeyron, Main.ThermoPlots
using CoolProp
import PlotlyBase, PlotlyJS, PlotlyKaleido
@genietools
@app begin
    @in species1 = "water"
    @in species2 = "ethanol"
    @in Select_eos = "PCSAFT"
    @in pre = 1.
    @in temp = 298.15
    @in new_button = false
    @in y_axis = "Volume"

    @out Select_eos_list = ["PCSAFT","SAFTVRMie","SAFTγMie","PR","RK","vdW","UNIFAC","MultiFluid"]
    @out Select_property = ["Volume",
                            "Internal Energy","Enthalpy","Entropy","Gibbs Free Energy",
                            "Isobaric Heat Capacity","Isochoric Heat Capacity",
                            "Speed of Sound"]
    @out trace = []
    @out layout = PlotlyBase.Layout(autosize=false,width=700,height=470,
    xaxis = PlotlyBase.attr(title = "water composition / (mol/mol)", font_size=12, showgrid=false,            
                                      ticks="inside",mirror=true,showline=true,linecolor="black",type="linear"),
                         yaxis = PlotlyBase.attr(title = "Excess Volume / (dm³/mol)", font_size=12, showgrid=false,       
                                      ticks="inside",mirror=true,showline=true,linecolor="black"),
                         showlegend=false, plot_bgcolor="white")
    @out v_E = []
    @out u_E = []
    @out h_E = []
    @out s_E = []
    @out g_E = []
    @out cp_E = []
    @out cv_E = []
    @out w_E = []

    @out y_label = "Excess Volume / (dm³/mol)"
    @out x = []
    @out y = []

    @onbutton new_button begin
        Npoints = 200
        eos = Symbol(Select_eos)
        model = @eval $eos([$species1,$species2])

        x = LinRange(0,1,Npoints)
        v_E = zeros(Npoints)
        u_E = zeros(Npoints)
        h_E = zeros(Npoints)
        s_E = zeros(Npoints)
        g_E = zeros(Npoints)
        cp_E = zeros(Npoints)
        cv_E = zeros(Npoints)
        w_E = zeros(Npoints)

        for i in 1:Npoints
            v_E[i] = excess(model,pre*1e5,temp,[x[i],1-x[i]],volume)
            u_E[i] = excess(model,pre*1e5,temp,[x[i],1-x[i]],internal_energy)
            h_E[i] = excess(model,pre*1e5,temp,[x[i],1-x[i]],enthalpy)
            s_E[i] = excess(model,pre*1e5,temp,[x[i],1-x[i]],entropy)
            g_E[i] = excess(model,pre*1e5,temp,[x[i],1-x[i]],gibbs_free_energy)
            cp_E[i] = excess(model,pre*1e5,temp,[x[i],1-x[i]],isobaric_heat_capacity)
            cv_E[i] = excess(model,pre*1e5,temp,[x[i],1-x[i]],isochoric_heat_capacity)
            w_E[i] = excess(model,pre*1e5,temp,[x[i],1-x[i]],speed_of_sound)
        end

        if y_axis == "Volume"
            y = v_E.*1e3
            y_label = "Volume / (dm³/mol)"
        elseif y_axis == "Internal Energy"
            y = u_E
            y_label = "Internal Energy / (J/mol)"
        elseif y_axis == "Enthalpy"
            y = h_E
            y_label = "Enthalpy / (J/mol)"
        elseif y_axis == "Entropy"
            y = s_E
            y_label = "Entropy / (J/mol K)"
        elseif y_axis == "Gibbs Free Energy"
            y = g_E
            y_label = "Gibbs Free Energy / (J/mol)"
        elseif y_axis == "Isobaric Heat Capacity"
            y = cp_E
            y_label = "Isobaric Heat Capacity / (J/mol K)"
        elseif y_axis == "Isochoric Heat Capacity"
            y = cv_E
            y_label = "Isochoric Heat Capacity / (J/mol K)"
        elseif y_axis == "Speed of Sound"
            y = w_E
            y_label = "Speed of Sound / (m/s)"
        end

        y_label = "Excess " * y_label

        trace = [PlotlyBase.scatter(x=x,y=y,mode="lines",line=PlotlyBase.attr(color="red", dash="solid", width=3),name="")]
        layout = PlotlyBase.Layout(autosize=false,width=700,height=470,
        xaxis = PlotlyBase.attr(title = species1*" composition / (mol/mol)", font_size=12, showgrid=false,            
                                          ticks="inside",mirror=true,showline=true,linecolor="black"),
                             yaxis = PlotlyBase.attr(title = y_label, font_size=12, showgrid=false,       
                                          ticks="inside",mirror=true,showline=true,linecolor="black"),
                             showlegend=false, plot_bgcolor="white")
    end
    @onchange y_axis begin
        if y_axis == "Volume"
            y = v_E.*1e3
            y_label = "Volume / (dm³/mol)"
        elseif y_axis == "Internal Energy"
            y = u_E
            y_label = "Internal Energy / (J/mol)"
        elseif y_axis == "Enthalpy"
            y = h_E
            y_label = "Enthalpy / (J/mol)"
        elseif y_axis == "Entropy"
            y = s_E
            y_label = "Entropy / (J/mol K)"
        elseif y_axis == "Gibbs Free Energy"
            y = g_E
            y_label = "Gibbs Free Energy / (J/mol)"
        elseif y_axis == "Isobaric Heat Capacity"
            y = cp_E
            y_label = "Isobaric Heat Capacity / (J/mol K)"
        elseif y_axis == "Isochoric Heat Capacity"
            y = cv_E
            y_label = "Isochoric Heat Capacity / (J/mol K)"
        elseif y_axis == "Speed of Sound"
            y = w_E
            y_label = "Speed of Sound / (m/s)"
        end

        y_label = "Excess " * y_label

        trace = [PlotlyBase.scatter(x=x,y=y,mode="lines",line=PlotlyBase.attr(color="red", dash="solid", width=3),name="")]
        layout = PlotlyBase.Layout(autosize=false,width=700,height=470,
        xaxis = PlotlyBase.attr(title = species1*" composition / (mol/mol)", font_size=12, showgrid=false,            
                                          ticks="inside",mirror=true,showline=true,linecolor="black"),
                             yaxis = PlotlyBase.attr(title = y_label, font_size=12, showgrid=false,       
                                          ticks="inside",mirror=true,showline=true,linecolor="black"),
                             showlegend=false, plot_bgcolor="white")
    end
end
end