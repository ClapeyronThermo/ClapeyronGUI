module PURERHOT
using GenieFramework
using Clapeyron, Main.ThermoPlots
using CoolProp
import PlotlyBase, PlotlyKaleido

@app begin
    @in species = "water"
    @in Select_eos = "PCSAFT"
    @in new_button = false
    @in add_button = false
    @out Select_eos_list = ["PCSAFT","SAFTVRMie","SAFTγMie","PR","RK","vdW","MultiFluid"]
    @out color = ["red","blue","green","purple","black"]
    @out i = 1
    @out trace = []
    @out layout = PlotlyBase.Layout(xaxis = PlotlyBase.attr(title = "Density / (mol/dm³)", font_size=12, showgrid=false,            
                                      ticks="inside",mirror=true,showline=true,linecolor="black"),
                         yaxis = PlotlyBase.attr(title = "Temperature  / K", font_size=12, showgrid=false,       
                                      ticks="inside",mirror=true,showline=true,linecolor="black"),
                         showlegend=false, plot_bgcolor="white")
    @onbutton new_button begin
        i = 1
        eos = Symbol(Select_eos)
        model = @eval $eos([$species])

        plt = rhoT_envelope(model; color=color[mod(7,i)+1])
        trace = plt.plot.data
        layout = plt.plot.layout
    end

    @onbutton add_button begin
        i += 1
        eos = Symbol(Select_eos)
        model = @eval $eos([$species])
        rhoT_envelope!(plt,model; color=color[mod(7,i)+1])
        trace = plt.plot.data
        layout = plt.plot.layout
    end
end
end