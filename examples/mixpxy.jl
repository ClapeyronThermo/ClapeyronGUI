module MIXPXY
using GenieFramework, StippleUI
using Clapeyron, Main.ThermoPlots
import PlotlyBase, PlotlyKaleido, Plots
Plots.plotly()

@app begin
    @in species1 = "methanol"
    @in species2 = "hexane"
    @in temp = 298.15
    @in Select_eos = "PCSAFT"
    @in new_button = false
    @in add_button = false
    model = PCSAFT(["methanol","hexane"])
    T = 298.15
    plt = pxy_diagram(model, T)
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
        model = @eval $eos([$species1,$species2])

        plt = pxy_diagram(model, temp; color=color[mod(7,i)+1])
        trace = Plots.plotly_traces(plt)
        layout = Plots.plotly_layout(plt)
    end

    @onbutton add_button begin
        i += 1
        eos = Symbol(Select_eos)
        model = @eval $eos([$species1,$species2])
        pxy_diagram!(plt,model, temp; color=color[mod(7,i)+1])
        trace = Plots.plotly_traces(plt)
        layout = Plots.plotly_layout(plt)
    end
end
end