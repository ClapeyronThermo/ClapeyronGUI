module TERNARY_DIAGRAM
using GenieFramework
using Clapeyron, Main.ThermoPlots
using CoolProp
import PlotlyBase, PlotlyKaleido
import Main.@timeout
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

        try
            model = @eval $eos([$species1,$species2,$species3])
        catch
            notify(__model__, "Species $species1 or $species2 or $species3 are not available in $Select_eos.", :warning)
            throw(TypeError("Species $species1 or $species2 or $species3 are not available in $Select_eos."))
        end

        model = @eval $eos([$species1,$species2,$species3])
        
        if temp < 0
            notify(__model__, "Temperature must be positive.", :warning)
            throw(TypeError("Temperature must be positive."))
        end

        if pre < 0
            notify(__model__, "Pressure must be positive.", :warning)
            throw(TypeError("Pressure must be positive."))
        end

        plt = @timeout 200 ternary_diagram(model, pre*1e5, temp; Npoints=Npoints, color=:blue, style=:solid, check_three_phase=false)
        trace = plt.data
        layout = plt.layout
    end
end
end