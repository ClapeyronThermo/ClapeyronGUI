module CDFT
using GenieFramework
using Clapeyron, Main.cDFT
import PlotlyBase, PlotlyJS, PlotlyKaleido
import Main.@timeout
@genietools

@app begin
    @in tab_selected = "sft"
    @in species = "water"
    @in temp = 298.15
    @in surfactant = "ethanol"
    @in Select_eos = "PCSAFT"

    @in converge_sft_button = false

    @in converge_ift_button = false

    @out eos_list = ["PCSAFT","SAFTVRMie"]
    @out sft = 0.
    @out ift = 0.

    @out trace_sft = []
    @out layout_sft = PlotlyBase.Layout(autosize=false,width=700,height=470,
                            xaxis = PlotlyBase.attr(title = "Length  / Å", font_size=12, showgrid=false,            
                                          ticks="inside",mirror=true,showline=true,linecolor="black"),
                            yaxis = PlotlyBase.attr(title = "Density / (kg/m³)", font_size=12, showgrid=false,       
                                          ticks="inside",mirror=true,showline=true,linecolor="black"),
                             showlegend=false, plot_bgcolor="white")
    @out trace_ift = []
    @out layout_ift = PlotlyBase.Layout(autosize=false,width=700,height=470,
                            xaxis = PlotlyBase.attr(title = "Length  / Å", font_size=12, showgrid=false,            
                                        ticks="inside",mirror=true,showline=true,linecolor="black"),
                            yaxis = PlotlyBase.attr(title = "Density / (kg/m³)", font_size=12, showgrid=false,       
                                        ticks="inside",mirror=true,showline=true,linecolor="black"),
                            showlegend=false, plot_bgcolor="white")

    @onbutton converge_sft_button begin
        eos = Symbol(Select_eos)

        try 
            model = @eval $eos([$species])
        catch
            notify(__model__, "Species $species is not available in $Select_eos.", :warning)
            throw(TypeError("Species $species is not available in $Select_eos."))
        end

        model = @eval $eos([$species])

        if temp < 0
            notify(__model__, "Temperature must be positive.", :warning)
            throw(TypeError("Temperature must be positive."))
        end

        Mw = model.params.Mw.values

        L = length_scale(model)

        (p, vl, vv) = saturation_pressure(model, temp)
    
        structure = SurfaceTension1DCart((p, temp, [1.]),[-10L,10L], 101)

        system = DFTSystem(model, structure)
        converge!(system)
        sft = round(surface_tension(system)*1e3, digits=2)

        x = LinRange(-10L,10L,1001)
        y = system.profiles[1].(x).*Mw[1].*1e-3
        ymax = maximum(y)*1.05
        trace_sft = [PlotlyBase.scatter(x=x./1e-10, y=y, mode="lines", line=PlotlyBase.attr(color="blue", width=3), name="")]
        layout_sft = PlotlyBase.Layout(autosize=false,width=700,height=470,
                            xaxis = PlotlyBase.attr(title = "Length  / Å", font_size=12, showgrid=false,            
                                          ticks="inside",mirror=true,showline=true,linecolor="black"),
                            yaxis = PlotlyBase.attr(title = "Density / (kg/m³)", font_size=12, showgrid=false,       
                                          ticks="inside",mirror=true,showline=true,linecolor="black"),
                             showlegend=false, plot_bgcolor="white", yaxis_range=[0,ymax])
    end

    @onbutton converge_ift_button begin
        eos = Symbol(Select_eos)

        try
            model = @eval $eos(["water","hexane",$surfactant]; assoc_options=AssocOptions(combining=:elliott))
        catch
            notify(__model__, "Species $surfactant are not available in $Select_eos.", :warning)
            throw(TypeError("Species $surfactant are not available in $Select_eos."))
        end

        model = @eval $eos(["water","hexane",$surfactant]; assoc_options=AssocOptions(combining=:elliott))

        Mw = model.params.Mw.values

        L = length_scale(model)

        (x,n,G) = tp_flash(model, 1e5, 298.15, [0.5,0.5,1e-3], RRTPFlash(equilibrium=:lle))

        structure = InterfacialTension1DCart((1e5, 298.15, x[1,:]),[-10L,10L], 201, x[2,:])        
        system = DFTSystem(model, structure)
        converge!(system)
        ift = round(interfacial_tension(system)*1e3, digits=2)

        x = LinRange(-10L,10L,1001)
        y1 = system.profiles[1].(x).*Mw[1].*1e-3
        y2 = system.profiles[2].(x).*Mw[2].*1e-3
        y3 = system.profiles[3].(x).*Mw[3].*1e-3

        ymax = maximum([maximum(y1),maximum(y2),maximum(y3)])*1.05

        trace_ift = [PlotlyBase.scatter(x=x./1e-10, y=y1, mode="lines", line=PlotlyBase.attr(color="blue", width=3), name=""),
                     PlotlyBase.scatter(x=x./1e-10, y=y2, mode="lines", line=PlotlyBase.attr(color="red", width=3), name=""),
                     PlotlyBase.scatter(x=x./1e-10, y=y3, mode="lines", line=PlotlyBase.attr(color="green", width=3), name="")]

        layout_ift = PlotlyBase.Layout(autosize=false,width=700,height=470,
        xaxis = PlotlyBase.attr(title = "Length  / Å", font_size=12, showgrid=false,            
                    ticks="inside",mirror=true,showline=true,linecolor="black"),
        yaxis = PlotlyBase.attr(title = "Density / (kg/m³)", font_size=12, showgrid=false,       
                    ticks="inside",mirror=true,showline=true,linecolor="black"),
        showlegend=false, plot_bgcolor="white", yaxis_range=[0,ymax])
    end
end
end