module BULK
using GenieFramework
using Clapeyron, Main.ThermoPlots
using CoolProp
import PlotlyBase, PlotlyKaleido
@genietools
@app begin
    @in tab_selected = "isobaric"

    @in species = "water"
    @in Select_eos = "PCSAFT"

    @in new_T_button = false
    @in new_p_button = false

    @in x_T_axis = "Pressure"
    @in y_T_axis = "Density"
    @in x_p_axis = "Temperature"
    @in y_p_axis = "Density"

    @in log_x_T = false
    @in log_y_T = false
    @in log_x_p = false
    @in log_y_p = false

    @out Select_eos_list = ["PCSAFT","SAFTVRMie","SAFTγMie","PR","RK","vdW","MultiFluid"]
    @out Select_property = ["Density","Volume",
                            "Internal Energy","Enthalpy","Entropy",
                            "Isobaric Heat Capacity","Isochoric Heat Capacity",
                            "Speed of Sound","Joule-Thomson Coefficient",
                            "Isentropic Compressibility","Isothermal Compressibility",
                            "Isobaric Expansivity"]
    @out trace_T = []
    @out layout_T = PlotlyBase.Layout(autosize=false,width=700,height=470,
    xaxis = PlotlyBase.attr(title = "Pressure  / bar", font_size=12, showgrid=false,            
                                      ticks="inside",mirror=true,showline=true,linecolor="black"),
                         yaxis = PlotlyBase.attr(title = "Density / (mol/dm³)", font_size=12, showgrid=false,       
                                      ticks="inside",mirror=true,showline=true,linecolor="black"),
                         showlegend=false, plot_bgcolor="white")

    @out trace_p = []
    @out layout_p = PlotlyBase.Layout(autosize=false,width=700,height=470,
    xaxis = PlotlyBase.attr(title = "Temperature  / K", font_size=12, showgrid=false,            
                                    ticks="inside",mirror=true,showline=true,linecolor="black"),
                        yaxis = PlotlyBase.attr(title = "Density / (mol/dm³)", font_size=12, showgrid=false,       
                                    ticks="inside",mirror=true,showline=true,linecolor="black"),
                        showlegend=false, plot_bgcolor="white")
    @in temp = 298.15
    @in pre_start = 0.1
    @in pre_end = 10.
    @out v_T = []
    @out ρ_T = []
    @out u_T = []
    @out h_T = []
    @out s_T = []
    @out cp_T = []
    @out cv_T = []  
    @out w_T = []
    @out jt_T = []
    @out ks_T = []
    @out kt_T = []
    @out β_T = []
    @out y_T_label = "Density / (mol/dm³)"
    @out x_T = []
    @out y_T = []

    @in pre = 1.
    @in temp_start = 273.15
    @in temp_end = 423.15
    @out v_p = []
    @out ρ_p = []
    @out u_p = []
    @out h_p = []
    @out s_p = []
    @out cp_p = []
    @out cv_p = []
    @out w_p = []
    @out jt_p = []
    @out ks_p = []
    @out kt_p = []
    @out β_p = []
    @out y_p_label = "Density / (mol/dm³)"
    @out x_p = []
    @out y_p = []

    @onbutton new_T_button begin
        Npoints = 200
        eos = Symbol(Select_eos)
        try
            model = @eval $eos([$species])
        catch
            notify(__model__, "Species $species is not available in $Select_eos.", :warning)
            throw(TypeError("Species $species is not available in $Select_eos."))
        end

        model = @eval $eos([$species])

        if !(typeof(pre_start)<:Number) || !(typeof(pre_end)<:Number)
            notify(__model__, "Pressure range must be a number.", :warning)
            throw(TypeError("Pressure range must be a number."))
        elseif pre_start < 0 || pre_end < 0
            notify(__model__, "Pressure range must be positive.", :warning)
            throw(TypeError("Pressure range must be positive."))
        end

        if !(typeof(temp) <: Number)
            notify(__model__, "Temperature must be a number.", :warning)
            throw(TypeError("Temperature must be a number."))
        elseif temp < 0
            notify(__model__, "Temperature must be positive.", :warning)
            throw(TypeError("Temperature must be positive."))
        end

        if log_x_T
            p = LinRange(log(pre_start),log(pre_end),Npoints)
            p = exp.(p)
        else
            p = LinRange(pre_start,pre_end,Npoints)
        end
        p = p.*1e5
        v_T = volume.(model,p,temp)
        ρ_T = 1 ./v_T
        u_T = internal_energy.(model,p,temp)
        h_T = enthalpy.(model,p,temp)
        s_T = entropy.(model,p,temp)
        cp_T = isobaric_heat_capacity.(model,p,temp)
        cv_T = isochoric_heat_capacity.(model,p,temp)
        w_T = speed_of_sound.(model,p,temp)
        jt_T = joule_thomson_coefficient.(model,p,temp)
        ks_T = isentropic_compressibility.(model,p,temp)
        kt_T = isothermal_compressibility.(model,p,temp)
        β_T = isobaric_expansivity.(model,p,temp)

        x_T = p./1e5

        if y_T_axis == "Density"
            y_T = ρ_T.*1e-3
            y_T_label = "Density / (mol/dm³)"
        elseif y_T_axis == "Volume"
            y_T = v_T.*1e3
            y_T_label = "Volume / (dm³/mol)"
        elseif y_T_axis == "Internal Energy"
            y_T = u_T
            y_T_label = "Internal Energy / (J/mol)"
        elseif y_T_axis == "Enthalpy"
            y_T = h_T
            y_T_label = "Enthalpy / (J/mol)"
        elseif y_T_axis == "Entropy"
            y_T = s_T
            y_T_label = "Entropy / (J/mol K)"
        elseif y_T_axis == "Isobaric Heat Capacity"
            y_T = cp_T
            y_T_label = "Isobaric Heat Capacity / (J/mol K)"
        elseif y_T_axis == "Isochoric Heat Capacity"
            y_T = cv_T
            y_T_label = "Isochoric Heat Capacity / (J/mol K)"
        elseif y_T_axis == "Speed of Sound"
            y_T = w_T
            y_T_label = "Speed of Sound / (m/s)"
        elseif y_T_axis == "Joule-Thomson Coefficient"
            y_T = jt_T.*1e5
            y_T_label = "Joule-Thomson Coefficient / (K/bar)"
        elseif y_T_axis == "Isentropic Compressibility"
            y_T = ks_T./1e5
            y_T_label = "Isentropic Compressibility / bar"
        elseif y_T_axis == "Isothermal Compressibility"
            y_T = kt_T./1e5
            y_T_label = "Isothermal Compressibility / bar"
        elseif y_T_axis == "Isobaric Expansivity"
            y_T = β_T
            y_T_label = "Isobaric Expansivity / K⁻¹"
        end

        if log_x_T
            type_x = "log"
        else
            type_x = "linear"
        end

        if log_y_T
            type_y = "log"
        else
            type_y = "linear"
        end

        trace_T = [PlotlyBase.scatter(x=x_T,y=y_T,mode="lines",line=PlotlyBase.attr(color="red", dash="solid", width=3),name="")]
        layout_T = PlotlyBase.Layout(autosize=false,width=700,height=470,
        xaxis = PlotlyBase.attr(title = "Pressure / bar", font_size=12, showgrid=false,            
                                          ticks="inside",mirror=true,showline=true,linecolor="black",type=type_x),
                             yaxis = PlotlyBase.attr(title = y_T_label, font_size=12, showgrid=false,       
                                          ticks="inside",mirror=true,showline=true,linecolor="black",type=type_y),
                             showlegend=false, plot_bgcolor="white")
    end

    @onchange y_T_axis begin
        if y_T_axis == "Density"
            y_T = ρ_T.*1e-3
            y_T_label = "Density / (mol/dm³)"
        elseif y_T_axis == "Volume"
            y_T = v_T.*1e3
            y_T_label = "Volume / (dm³/mol)"
        elseif y_T_axis == "Internal Energy"
            y_T = u_T
            y_T_label = "Internal Energy / (J/mol)"
        elseif y_T_axis == "Enthalpy"
            y_T = h_T
            y_T_label = "Enthalpy / (J/mol)"
        elseif y_T_axis == "Entropy"
            y_T = s_T
            y_T_label = "Entropy / (J/mol K)"
        elseif y_T_axis == "Isobaric Heat Capacity"
            y_T = cp_T
            y_T_label = "Isobaric Heat Capacity / (J/mol K)"
        elseif y_T_axis == "Isochoric Heat Capacity"
            y_T = cv_T
            y_T_label = "Isochoric Heat Capacity / (J/mol K)"
        elseif y_T_axis == "Speed of Sound"
            y_T = w_T
            y_T_label = "Speed of Sound / (m/s)"
        elseif y_T_axis == "Joule-Thomson Coefficient"
            y_T = jt_T.*1e5
            y_T_label = "Joule-Thomson Coefficient / (K/bar)"
        elseif y_T_axis == "Isentropic Compressibility"
            y_T = ks_T./1e5
            y_T_label = "Isentropic Compressibility / bar"
        elseif y_T_axis == "Isothermal Compressibility"
            y_T = kt_T./1e5
            y_T_label = "Isothermal Compressibility / bar"
        elseif y_T_axis == "Isobaric Expansivity"
            y_T = β_T
            y_T_label = "Isobaric Expansivity / K⁻¹"
        end

        if log_x_T
            type_x = "log"
        else
            type_x = "linear"
        end

        if log_y_T
            type_y = "log"
        else
            type_y = "linear"
        end

        trace_T = [PlotlyBase.scatter(x=x_T,y=y_T,mode="lines",line=PlotlyBase.attr(color="red", dash="solid", width=3),name="")]
        layout_T = PlotlyBase.Layout(autosize=false,width=700,height=470,
        xaxis = PlotlyBase.attr(title = "Pressure / bar", font_size=12, showgrid=false,            
                                          ticks="inside",mirror=true,showline=true,linecolor="black",type=type_x),
                             yaxis = PlotlyBase.attr(title = y_T_label, font_size=12, showgrid=false,       
                                          ticks="inside",mirror=true,showline=true,linecolor="black",type=type_y),
                             showlegend=false, plot_bgcolor="white")
    end

    @onchange log_x_T, log_y_T begin
        if log_x_T
            type_x = "log"
        else
            type_x = "linear"
        end

        if log_y_T
            type_y = "log"
        else
            type_y = "linear"
        end

        layout_T = PlotlyBase.Layout(autosize=false,width=700,height=470,
        xaxis = PlotlyBase.attr(title = "Pressure / bar", font_size=12, showgrid=false,            
                                          ticks="inside",mirror=true,showline=true,linecolor="black",type=type_x),
                             yaxis = PlotlyBase.attr(title = y_T_label, font_size=12, showgrid=false,       
                                          ticks="inside",mirror=true,showline=true,linecolor="black",type=type_y),
                             showlegend=false, plot_bgcolor="white")

        trace_T = deepcopy(trace_T)
    end

    @onbutton new_p_button begin
        Npoints = 200
        eos = Symbol(Select_eos)
        try
            model = @eval $eos([$species])
        catch 
            notify(__model__, "Species $species is not available in $Select_eos.", :warning)
            throw(TypeError("Species $species is not available in $Select_eos."))
        end

        model = @eval $eos([$species])

        if !(typeof(temp_start)<:Number) || !(typeof(temp_end)<:Number)
            notify(__model__, "Temperature range must be a number.", :warning)
            throw(TypeError("Temperature range must be a number."))
        elseif temp_start < 0 || temp_end < 0
            notify(__model__, "Temperature range must be positive.", :warning)
            throw(TypeError("Temperature range must be positive."))
        end

        if !(typeof(pre) <: Number)
            notify(__model__, "Pressure must be a number.", :warning)
            throw(TypeError("Pressure must be a number."))
        elseif pre < 0
            notify(__model__, "Pressure must be positive.", :warning)
            throw(TypeError("Pressure must be positive."))
        end

        T = LinRange(temp_start,temp_end,Npoints)
        p = pre*1e5
        v_p = volume.(model,p,T)
        ρ_p = 1 ./v_p
        u_p = internal_energy.(model,p,T)
        h_p = enthalpy.(model,p,T)
        s_p = entropy.(model,p,T)
        cp_p = isobaric_heat_capacity.(model,p,T)
        cv_p = isochoric_heat_capacity.(model,p,T)
        w_p = speed_of_sound.(model,p,T)
        jt_p = joule_thomson_coefficient.(model,p,T)
        ks_p = isentropic_compressibility.(model,p,T)
        kt_p = isothermal_compressibility.(model,p,T)
        β_p = isobaric_expansivity.(model,p,T)

        x_p = T

        if y_p_axis == "Density"
            y_p = ρ_p.*1e-3
            y_p_label = "Density / (mol/dm³)"
        elseif y_p_axis == "Volume"
            y_p = v_p.*1e3
            y_p_label = "Volume / (dm³/mol)"
        elseif y_p_axis == "Internal Energy"
            y_p = u_p
            y_p_label = "Internal Energy / (J/mol)"
        elseif y_p_axis == "Enthalpy"
            y_p = h_p
            y_p_label = "Enthalpy / (J/mol)"
        elseif y_p_axis == "Entropy"
            y_p = s_p
            y_p_label = "Entropy / (J/mol K)"
        elseif y_p_axis == "Isobaric Heat Capacity"
            y_p = cp_p
            y_p_label = "Isobaric Heat Capacity / (J/mol K)"
        elseif y_p_axis == "Isochoric Heat Capacity"
            y_p = cv_p
            y_p_label = "Isochoric Heat Capacity / (J/mol K)"
        elseif y_p_axis == "Speed of Sound"
            y_p = w_p
            y_p_label = "Speed of Sound / (m/s)"
        elseif y_p_axis == "Joule-Thomson Coefficient"
            y_p = jt_p.*1e5
            y_p_label = "Joule-Thomson Coefficient / (K/bar)"
        elseif y_p_axis == "Isentropic Compressibility"
            y_p = ks_p./1e5
            y_p_label = "Isentropic Compressibility / bar"
        elseif y_p_axis == "Isothermal Compressibility"
            y_p = kt_p./1e5
            y_p_label = "Isothermal Compressibility / bar"
        elseif y_p_axis == "Isobaric Expansivity"
            y_p = β_p
            y_p_label = "Isobaric Expansivity / K⁻¹"
        end

        if log_x_p
            type_x = "log"
        else
            type_x = "linear"
        end

        if log_y_p
            type_y = "log"
        else
            type_y = "linear"
        end

        trace_p = [PlotlyBase.scatter(x=x_p,y=y_p,mode="lines",line=PlotlyBase.attr(color="blue", dash="solid", width=3),name="")]
        layout_p = PlotlyBase.Layout(autosize=false,width=700,height=470,
        xaxis = PlotlyBase.attr(title = "Temperature / K", font_size=12, showgrid=false,            
                                          ticks="inside",mirror=true,showline=true,linecolor="black",type=type_x),
                             yaxis = PlotlyBase.attr(title = y_p_label, font_size=12, showgrid=false,       
                                          ticks="inside",mirror=true,showline=true,linecolor="black",type=type_y),
                             showlegend=false, plot_bgcolor="white")
    end

    @onchange y_p_axis begin
        if y_p_axis == "Density"
            y_p = ρ_p.*1e-3
            y_p_label = "Density / (mol/dm³)"
        elseif y_p_axis == "Volume"
            y_p = v_p.*1e3
            y_p_label = "Volume / (dm³/mol)"
        elseif y_p_axis == "Internal Energy"
            y_p = u_p
            y_p_label = "Internal Energy / (J/mol)"
        elseif y_p_axis == "Enthalpy"
            y_p = h_p
            y_p_label = "Enthalpy / (J/mol)"
        elseif y_p_axis == "Entropy"
            y_p = s_p
            y_p_label = "Entropy / (J/mol K)"
        elseif y_p_axis == "Isobaric Heat Capacity"
            y_p = cp_p
            y_p_label = "Isobaric Heat Capacity / (J/mol K)"
        elseif y_p_axis == "Isochoric Heat Capacity"
            y_p = cv_p
            y_p_label = "Isochoric Heat Capacity / (J/mol K)"
        elseif y_p_axis == "Speed of Sound"
            y_p = w_p
            y_p_label = "Speed of Sound / (m/s)"
        elseif y_p_axis == "Joule-Thomson Coefficient"
            y_p = jt_p.*1e5
            y_p_label = "Joule-Thomson Coefficient / (K/bar)"
        elseif y_p_axis == "Isentropic Compressibility"
            y_p = ks_p./1e5
            y_p_label = "Isentropic Compressibility / bar"
        elseif y_p_axis == "Isothermal Compressibility"
            y_p = kt_p./1e5
            y_p_label = "Isothermal Compressibility / bar"
        elseif y_p_axis == "Isobaric Expansivity"
            y_p = β_p
            y_p_label = "Isobaric Expansivity / K⁻¹"
        end

        if log_x_p
            type_x = "log"
        else
            type_x = "linear"
        end

        if log_y_p
            type_y = "log"
        else
            type_y = "linear"
        end

        trace_p = [PlotlyBase.scatter(x=x_p,y=y_p,mode="lines",line=PlotlyBase.attr(color="blue", dash="solid", width=3),name="")]
        layout_p = PlotlyBase.Layout(autosize=false,width=700,height=470,
        xaxis = PlotlyBase.attr(title = "Temperature / K", font_size=12, showgrid=false,            
                                          ticks="inside",mirror=true,showline=true,linecolor="black",type=type_x),
                             yaxis = PlotlyBase.attr(title = y_p_label, font_size=12, showgrid=false,       
                                          ticks="inside",mirror=true,showline=true,linecolor="black",type=type_y),
                             showlegend=false, plot_bgcolor="white")
    end

    @onchange log_x_p, log_y_p begin
        if log_x_p
            type_x = "log"
        else
            type_x = "linear"
        end

        if log_y_p
            type_y = "log"
        else
            type_y = "linear"
        end

        layout_p = PlotlyBase.Layout(autosize=false,width=700,height=470,
        xaxis = PlotlyBase.attr(title = "Temperature / K", font_size=12, showgrid=false,            
                                          ticks="inside",mirror=true,showline=true,linecolor="black",type=type_x),
                             yaxis = PlotlyBase.attr(title = y_p_label, font_size=12, showgrid=false,       
                                          ticks="inside",mirror=true,showline=true,linecolor="black",type=type_y),
                             showlegend=false, plot_bgcolor="white")

        trace_p = deepcopy(trace_p)
    end
end
end