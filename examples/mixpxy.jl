module MIXPXY
using GenieFramework, StippleUI
using Clapeyron, Main.ThermoPlots
using CoolProp
import PlotlyBase, PlotlyKaleido

@app begin
    @in tab_selected = "pxy"
    @in species1 = "methanol"
    @in species2 = "hexane"
    @in temp = 298.15
    @in pre = 1.
    @in check_lle = false
    @in Select_eos = "PCSAFT"
    @in new_p_button = false
    @in new_T_button = false
    @in new_pT_button = false
    @in log_y_pT = false
    model = PCSAFT(["methanol","hexane"])
    @out Select_eos_list = ["PCSAFT","SAFTVRMie","SAFTγMie","PR","RK","vdW","MultiFluid","Wilson","NRTL","UNIFAC","COSMOSAC"]
    @out color = ["red","blue","green","purple","black"]
    @out i = 1
    @out trace_T = []
    @out layout_T = PlotlyBase.Layout(autosize=false,width=700,height=470,
    xaxis = PlotlyBase.attr(title = "methanol composition / (mol/mol)", font_size=12, showgrid=false,            
                                      ticks="inside",mirror=true,showline=true,linecolor="black"),
                         yaxis = PlotlyBase.attr(title = "Pressure / bar", font_size=12, showgrid=false,       
                                      ticks="inside",mirror=true,showline=true,linecolor="black"),
                         showlegend=false, plot_bgcolor="white")
    @out trace_p = []
    @out layout_p = PlotlyBase.Layout(autosize=false,width=700,height=470,
    xaxis = PlotlyBase.attr(title = "methanol composition / (mol/mol)", font_size=12, showgrid=false,            
                                      ticks="inside",mirror=true,showline=true,linecolor="black"),
                         yaxis = PlotlyBase.attr(title = "Temperature / K", font_size=12, showgrid=false,       
                                      ticks="inside",mirror=true,showline=true,linecolor="black"),
                         showlegend=false, plot_bgcolor="white")
    @out trace_pT = []
    @out layout_pT = PlotlyBase.Layout(autosize=false,width=700,height=470,
    xaxis = PlotlyBase.attr(title = "Temperature / K", font_size=12, showgrid=false,            
                                      ticks="inside",mirror=true,showline=true,linecolor="black"),
                         yaxis = PlotlyBase.attr(title = "Pressure / bar", font_size=12, showgrid=false,       
                                      ticks="inside",mirror=true,showline=true,linecolor="black"),
                         showlegend=false, plot_bgcolor="white")
    @onbutton new_T_button begin
        i = 1
        eos = Symbol(Select_eos)
        model = @eval $eos([$species1,$species2])

        plt = pxy_diagram(model, temp; color="blue")
        trace_T = plt.plot.data
        layout_T = plt.plot.layout
    end

    @onbutton new_p_button begin
        i = 1
        eos = Symbol(Select_eos)
        model = @eval $eos([$species1,$species2])

        plt = Txy_diagram(model, pre*1e5; color="red", check_lle=check_lle)
        trace_p = plt.plot.data
        layout_p = plt.plot.layout
    end

    @onbutton new_pT_button begin
        i = 1
        eos = Symbol(Select_eos)
        model = @eval $eos([$species1,$species2])

        plt = pT_projection(model; color="purple")
        trace_pT = plt.plot.data
        layout_pT = plt.plot.layout
    end

end
end