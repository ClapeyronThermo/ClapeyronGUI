module App
using GenieFramework
using Clapeyron
import PlotlyBase, PlotlyKaleido, Plots
Plots.plotly()
include("pt_curve.jl") 

@app begin
    @in species = "water"
    @in Select_eos = "PCSAFT"
    @in new_button = false
    @in add_button = false
    model = PCSAFT(["water"])
    plt = pT_curve(model)
    # y log Scale
    Plots.yaxis!(plt, :log10)
    @out Select_eos_list = ["PCSAFT","SAFTVRMie","SAFTγMie","PR","RK","vdW","GERG2008","IAPWS95"]
    @out trace = Plots.plotly_traces(plt)
    @out layout = Plots.plotly_layout(plt)
    @onbutton new_button begin
        plt = []
        eos = Symbol(Select_eos)
        model = @eval $eos([$species])

        plt = pT_curve(model)
        Plots.yaxis!(plt, :log10)
        trace = Plots.plotly_traces(plt)
        layout = Plots.plotly_layout(plt)
    end

    @onbutton add_button begin
        eos = Symbol(Select_eos)
        model = @eval $eos([$species])
        pT_curve!(plt,model)
        Plots.yaxis!(plt, :log10)
        trace = Plots.plotly_traces(plt)
        layout = Plots.plotly_layout(plt)
    end
end

function ui()
    # row([
    # cell(class="st-col col-3", [
    #     h1("Saturation curve")])])
    
    row([cell(class="st-col col-3", [
        select(:Select_eos, options = :Select_eos_list, label = "Equation of State"),
        plot(:trace, layout=:layout)
        ]), 
        cell(class="st-col col-3", [
        textfield("Species:", :species)
       ]),
       cell(class="st-col col-3", [
        btn("New", @click(:new_button), loading=:new_button)
       ]),
       cell(class="st-col col-3", [
        btn("Add", @click(:add_button), loading=:add_button)
       ])
       ])
    
    # row([
    # cell(class="st-col col-3", [
    #     plot(:trace, layout=:layout)
    # ])])
end

@page("/", ui)
end