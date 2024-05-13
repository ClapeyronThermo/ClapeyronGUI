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
    T = 298.15
    plt = pxy_diagram(model, T)
    # y log Scale
    @out Select_eos_list = ["PCSAFT","SAFTVRMie","SAFTγMie","PR","RK","vdW","GERG2008","IAPWS95"]
    @out color = ["red","blue","green","purple","black"]
    @out i = 1
    @out trace = plt.plot.data
    @out layout = plt.plot.layout
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