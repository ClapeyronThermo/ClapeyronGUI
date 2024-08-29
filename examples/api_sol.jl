module API_SOL
using GenieFramework
using Clapeyron, Main.ThermoPlots
using CoolProp
import PlotlyBase, PlotlyKaleido
@genietools
@app begin
    @in tab_selected = "temperature"

    @in api = "paracetamol"
    @in solvent = "water"
    @in solvent1 = "water"
    @in solvent2 = "ethanol"
    @in temp = 298.15


    @in Select_eos = "PCSAFT"

    @in new_T_button = false
    @in new_p_button = false

    @in x_T_axis = "Pressure"
    @in y_T_axis = "Density"
    @in x_p_axis = "Temperature"
    @in y_p_axis = "Density"

    @in log_y_T = false
    @in log_y_p = false

    @out Select_eos_list = ["PCSAFT","SAFTVRMie","SAFTγMie","PR","RK","vdW","UNIFAC"]
    @out Select_api = ["ibuprofen","paracetamol"]
    @out Select_solvent = ["water","ethanol","methanol","acetone","benzene","toluene","hexane","octane","decane"]

    @out trace_T = []
    @out layout_T = PlotlyBase.Layout(autosize=false,width=700,height=470,
    xaxis = PlotlyBase.attr(title = "Temperature  / K", font_size=12, showgrid=false,            
                                      ticks="inside",mirror=true,showline=true,linecolor="black"),
                         yaxis = PlotlyBase.attr(title = "API Solubility / (mol/mol)", font_size=12, showgrid=false,       
                                      ticks="inside",mirror=true,showline=true,linecolor="black"),
                         showlegend=false, plot_bgcolor="white")

    @out trace_p = []
    @out layout_p = PlotlyBase.Layout(autosize=false,width=700,height=470,
    xaxis = PlotlyBase.attr(title = "water composition / (mol/mol)", font_size=12, showgrid=false,            
                                    ticks="inside",mirror=true,showline=true,linecolor="black"),
                        yaxis = PlotlyBase.attr(title = "API Solubility / (mol/mol)", font_size=12, showgrid=false,       
                                    ticks="inside",mirror=true,showline=true,linecolor="black"),
                        showlegend=false, plot_bgcolor="white")
    @in temp_start = 300.
    @in temp_end = 350.

    @onbutton new_T_button begin
        Npoints = 200
        eos = Symbol(Select_eos)
        try
            model = @eval CompositeModel([$solvent,$api];fluid=$eos,solid=SolidHfus)
        catch
            notify(__model__, "Species $solvent or $api are not available in $Select_eos.", :warning)
            throw(TypeError("Species $solvent or $api are not available in $Select_eos."))
        end

        model = @eval CompositeModel([$solvent,$api];fluid=$eos,solid=SolidHfus)

        p = 1e5

        if !(typeof(temp_start)<:Number) || !(typeof(temp_end)<:Number)
            notify(__model__, "Temperature range must be a number.", :warning)
            throw(TypeError("Temperature range must be a number."))
        elseif temp_start < 0 || temp_end < 0
            notify(__model__, "Temperature range must be positive.", :warning)
            throw(TypeError("Temperature range must be positive."))
        end

        T = LinRange(temp_start,temp_end,Npoints)

        s = zeros(length(T))

        for i in 1:length(T)
            s[i] = sle_solubility(model,p,T[i],[1.,0.];solute=[api])[2]
        end

        if log_y_T
            type_y = "log"
        else
            type_y = "linear"
        end

        trace_T = [PlotlyBase.scatter(x=T,y=s,mode="lines",line=PlotlyBase.attr(color="blue", dash="solid", width=3),name="")]
        layout_T = PlotlyBase.Layout(autosize=false,width=700,height=470,
        xaxis = PlotlyBase.attr(title = "Temperature / K", font_size=12, showgrid=false,            
                                          ticks="inside",mirror=true,showline=true,linecolor="black"),
                             yaxis = PlotlyBase.attr(title = "API Solubility / (mol/mol)", font_size=12, showgrid=false,       
                                          ticks="inside",mirror=true,showline=true,linecolor="black",type=type_y),
                             showlegend=false, plot_bgcolor="white")
    end

    @onchange log_y_T begin
        if log_y_T
            type_y = "log"
        else
            type_y = "linear"
        end

        layout_T = PlotlyBase.Layout(autosize=false,width=700,height=470,
        xaxis = PlotlyBase.attr(title = "Temperature / K", font_size=12, showgrid=false,            
                                          ticks="inside",mirror=true,showline=true,linecolor="black"),
                             yaxis = PlotlyBase.attr(title = "API Solubility / (mol/mol)", font_size=12, showgrid=false,       
                                          ticks="inside",mirror=true,showline=true,linecolor="black",type=type_y),
                             showlegend=false, plot_bgcolor="white")

        trace_T = deepcopy(trace_T)
    end

    @onbutton new_p_button begin
        Npoints = 200
        eos = Symbol(Select_eos)

        try
            model = @eval CompositeModel([$solvent1,$solvent2,$api];fluid=$eos,solid=SolidHfus)
        catch
            notify(__model__, "Species $solvent1, $solvent2 or $api are not available in $Select_eos.", :warning)
            throw(TypeError("Species $solvent1, $solvent2 or $api are not available in $Select_eos."))
        end

        model = @eval CompositeModel([$solvent1,$solvent2,$api];fluid=$eos,solid=SolidHfus)
        
        x = LinRange(0.,1.,Npoints)
        
        p = 1e5

        if !(typeof(temp)<:Number)
            notify(__model__, "Temperature must be a number.", :warning)
            throw(TypeError("Temperature range must be a number."))
        elseif temp < 0
            notify(__model__, "Temperature must be positive.", :warning)
            throw(TypeError("Temperature range must be positive."))
        end
        
        s = zeros(length(x))
        for i in 1:length(x)
            s[i] = sle_solubility(model,p,temp,[x[i],1-x[i],0.];solute=[api])[end]
        end

        if log_y_p
            type_y = "log"
        else
            type_y = "linear"
        end

        trace_p = [PlotlyBase.scatter(x=x,y=s,mode="lines",line=PlotlyBase.attr(color="blue", dash="solid", width=3),name="")]
        layout_p = PlotlyBase.Layout(autosize=false,width=700,height=470,
        xaxis = PlotlyBase.attr(title = solvent1*" composition / (mol/mol)", font_size=12, showgrid=false,            
                                          ticks="inside",mirror=true,showline=true,linecolor="black"),
                             yaxis = PlotlyBase.attr(title = "API Solubility / (mol/mol)", font_size=12, showgrid=false,       
                                          ticks="inside",mirror=true,showline=true,linecolor="black",type=type_y),
                             showlegend=false, plot_bgcolor="white")
    end

    @onchange log_y_p begin

        if log_y_p
            type_y = "log"
        else
            type_y = "linear"
        end

        layout_p = PlotlyBase.Layout(autosize=false,width=700,height=470,
        xaxis = PlotlyBase.attr(title = solvent1*" composition / (mol/mol)", font_size=12, showgrid=false,            
                                          ticks="inside",mirror=true,showline=true,linecolor="black"),
                             yaxis = PlotlyBase.attr(title = "API Solubility / (mol/mol)", font_size=12, showgrid=false,       
                                          ticks="inside",mirror=true,showline=true,linecolor="black",type=type_y),
                             showlegend=false, plot_bgcolor="white")

        trace_p = deepcopy(trace_p)
    end
end
end