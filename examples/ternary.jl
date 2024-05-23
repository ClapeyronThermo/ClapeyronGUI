module TERNARY_DIAGRAM
using GenieFramework
using Clapeyron, Main.ThermoPlots
using CoolProp
import PlotlyBase, PlotlyJS, PlotlyKaleido
@genietools

function make_ax(title, tickangle)
    PlotlyBase.attr(title=title, titlefont_size=20, tickangle=tickangle,
        tickfont_size=15, tickcolor="rgba(0, 0, 0, 0)", ticklen=5,
        showline=true, showgrid=true)
end

@app begin
    @in species1 = "water"
    @in species2 = "hexane"
    @in species3 = "ethanol"
    @in Select_eos = "SAFTgammaMie"
    @in pre = 1.
    @in temp = 298.15
    @in new_button = false

    @out Select_eos_list = ["PCSAFT","SAFTVRMie","SAFTγMie","PR","RK","UNIFAC"]

    @out trace = []
    @out layout = PlotlyBase.Layout(
                            ternary=PlotlyBase.attr(
                                sum=1,
                                aaxis=PlotlyBase.attr(title="water", titlefont_size=20, tickangle=0,
                                                    tickfont_size=15, tickcolor="rgb(0, 0, 0)", ticklen=5,
                                                    gridcolor="rgba(0, 0, 0)",linecolor="rgba(0, 0, 0)",
                                                    showline=true, showgrid=true),
                                baxis=PlotlyBase.attr(title="hexane", titlefont_size=20, tickangle=45,
                                                    tickfont_size=15, tickcolor="rgb(0, 0, 0)", ticklen=5,
                                                    gridcolor="rgba(0, 0, 0)",linecolor="rgba(0, 0, 0)",
                                                    showline=true, showgrid=true),
                                caxis=PlotlyBase.attr(title="ethanol", titlefont_size=20, tickangle=-45,
                                                    tickfont_size=15, tickcolor="rgb(0, 0, 0)", ticklen=5,
                                                    gridcolor="rgba(0, 0, 0)",linecolor="rgba(0, 0, 0)",
                                                    showline=true, showgrid=true),
                                bgcolor="#ffffff",
                            ))

    @onbutton new_button begin
        Npoints = 200
        eos = Symbol(Select_eos)
        model = @eval $eos([$species1,$species2,$species3])

        plt = ternary_diagram(model, pre*1e5, temp; Npoints=Npoints, color=:blue, style=:solid, check_three_phase=false)
        trace = plt.plot.data
        layout = plt.plot.layout
    end
end
end