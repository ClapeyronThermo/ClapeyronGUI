module MIXPXY
using GenieFramework, StippleUI
using Clapeyron, Main.ThermoPlots
import PlotlyBase, PlotlyKaleido

@app begin
    @in species1 = "methanol"
    @in species2 = "hexane"
    @in temp = 298.15
    @in Select_eos = "PCSAFT"
    @in new_button = false
    @in add_button = false
    model = PCSAFT(["methanol","hexane"])
    @out Select_eos_list = ["PCSAFT","SAFTVRMie","SAFTγMie","PR","RK","vdW","GERG2008","Wilson","NRTL","UNIFAC","COSMOSAC"]
    @out color = ["red","blue","green","purple","black"]
    @out i = 1
    @out trace = []
    @out layout = PlotlyBase.Layout(xaxis = PlotlyBase.attr(title = "Temperature  / K", font_size=12, showgrid=false,            
                                      ticks="inside",mirror=true,showline=true,linecolor="black"),
                         yaxis = PlotlyBase.attr(title = "Pressure / bar", font_size=12, showgrid=false,       
                                      ticks="inside",mirror=true,showline=true,linecolor="black"),
                         showlegend=false, plot_bgcolor="white")
    @onbutton new_button begin
        i = 1
        eos = Symbol(Select_eos)
        model = @eval $eos([$species1,$species2])

        plt = pxy_diagram(model, temp; color=color[mod(7,i)+1])
        trace = plt.plot.data
        layout = plt.plot.layout
    end

    @onbutton add_button begin
        i += 1
        eos = Symbol(Select_eos)
        model = @eval $eos([$species1,$species2])
        pxy_diagram!(plt,model, temp; color=color[mod(7,i)+1])
        trace = plt.plot.data
        layout = plt.plot.layout
    end
end
end