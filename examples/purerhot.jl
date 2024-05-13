module PURERHOT
using GenieFramework, StippleUI
using Clapeyron, Main.ThermoPlots
import PlotlyBase, PlotlyKaleido, Plots
Plots.plotly()

@app begin
    @in species = "water"
    @in Select_eos = "PCSAFT"
    @in new_button = false
    @in add_button = false
    model = PCSAFT(["water"])
    plt = rhoT_envelope(model)
    # y log Scale
    @out Select_eos_list = ["PCSAFT","SAFTVRMie","SAFTγMie","PR","RK","vdW","GERG2008","IAPWS95"]
    @out color = [:red,:blue,:green,:purple,:black]
    @out i = 1
    @out trace = Plots.plotly_traces(plt)
    @out layout = Plots.plotly_layout(plt)
    @onbutton new_button begin
        plt = []
        i = 1
        eos = Symbol(Select_eos)
        model = @eval $eos([$species])

        plt = rhoT_envelope(model; color=color[mod(7,i)+1])
        trace = Plots.plotly_traces(plt)
        layout = Plots.plotly_layout(plt)
    end

    @onbutton add_button begin
        i += 1
        eos = Symbol(Select_eos)
        model = @eval $eos([$species])
        rhoT_envelope!(plt,model; color=color[mod(7,i)+1])
        trace = Plots.plotly_traces(plt)
        layout = Plots.plotly_layout(plt)
    end
end
end