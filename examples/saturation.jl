module SATURATION
using GenieFramework
using Clapeyron, Main.ThermoPlots
using CoolProp
import PlotlyBase, PlotlyKaleido
@genietools
@app begin
    @in species = "water"
    @in Select_eos = "PCSAFT"
    @in new_button = false
    @in add_button = false
    @in x_axis = "Temperature"
    @in y_axis = "Pressure"
    @in log_x = false
    @in log_y = true
    @in exp_data = false

    @out Select_eos_list = ["PCSAFT","SAFTVRMie","SAFTγMie","PR","RK","vdW","MultiFluid"]
    @out Select_property = ["Temperature","Pressure","Density","Volume",
                            "Internal Energy","Enthalpy","Entropy",
                            "Isobaric Heat Capacity","Isochoric Heat Capacity",
                            "Speed of Sound","Joule-Thomson Coefficient",
                            "Isentropic Compressibility","Isothermal Compressibility",
                            "Isobaric Expansivity"]
    @out trace = []
    @out layout = PlotlyBase.Layout(autosize=false,width=700,height=470,
    xaxis = PlotlyBase.attr(title = "Temperature  / K", font_size=12, showgrid=false,            
                                      ticks="inside",mirror=true,showline=true,linecolor="black",type="linear"),
                         yaxis = PlotlyBase.attr(title = "Pressure / bar", font_size=12, showgrid=false,       
                                      ticks="inside",mirror=true,showline=true,linecolor="black",type="log"),
                         showlegend=false, plot_bgcolor="white")
    @out T = []
    @out psat = []
    @out vl_sat = []
    @out vv_sat = []
    @out ρl_sat = []
    @out ρv_sat = []
    @out ul_sat = []
    @out uv_sat = []
    @out hl_sat = []
    @out hv_sat = []
    @out sl_sat = []
    @out sv_sat = []
    @out cp_l_sat = []
    @out cp_v_sat = []
    @out cv_l_sat = []
    @out cv_v_sat = []
    @out w_l_sat = []
    @out w_v_sat = []
    @out jt_l_sat = []
    @out jt_v_sat = []
    @out ks_l_sat = []
    @out ks_v_sat = []
    @out kt_l_sat = []
    @out kt_v_sat = []
    @out β_l_sat = []
    @out β_v_sat = []

    @out Tc = 0.
    @out Pc = 0.
    @out vc = 0.
    @out ρc = 0.
    @out uc = 0.
    @out hc = 0.
    @out sc = 0.
    @out cp_c = 0.
    @out cv_c = 0.
    @out w_c = 0.
    @out jt_c = 0.
    @out ks_c = 0.
    @out kt_c = 0.
    @out β_c = 0.

    @out x_label = "Temperature  / K"
    @out y_label = "Pressure / bar"
    @out x = []
    @out y = []
    @out x_crit = 0.
    @out y_crit = 0.

    @onbutton new_button begin
        Npoints = 200
        eos = Symbol(Select_eos)

        try
            model = @eval $eos([$species])
        catch
            notify(__model__, "Species $species is not available in $Select_eos.", :warning)
            throw(TypeError("Species $species is not available in $Select_eos."))
        end

        model = @eval $eos([$species])

        (Tc,Pc,vc) = crit_pure(model)
        ρc = 1/vc
        uc = Clapeyron.VT_internal_energy.(model,vc,Tc)
        hc = Clapeyron.VT_enthalpy.(model,vc,Tc)
        sc = Clapeyron.VT_entropy.(model,vc,Tc)
        cp_c = Clapeyron.VT_isobaric_heat_capacity.(model,vc,Tc)
        cv_c = Clapeyron.VT_isochoric_heat_capacity.(model,vc,Tc)
        w_c = Clapeyron.VT_speed_of_sound.(model,vc,Tc)
        jt_c = Clapeyron.VT_joule_thomson_coefficient.(model,vc,Tc)
        ks_c = Clapeyron.VT_isentropic_compressibility.(model,vc,Tc)
        kt_c = Clapeyron.VT_isothermal_compressibility.(model,vc,Tc)
        β_c = Clapeyron.VT_isobaric_expansivity.(model,vc,Tc)

        T = LinRange(0.3*Tc,Tc,Npoints)
        sat = saturation_pressure.(model,T)

        psat = [sat[i][1] for i in 1:Npoints]
        vl_sat = [sat[i][2] for i in 1:Npoints]
        vv_sat = [sat[i][3] for i in 1:Npoints]
        T = T[.!isnan.(psat)]
        vl_sat = vl_sat[.!isnan.(psat)]
        vv_sat = vv_sat[.!isnan.(psat)]
        psat = psat[.!isnan.(psat)]

        ρl_sat = [1/vl for vl in vl_sat]
        ρv_sat = [1/vv for vv in vv_sat]
        ul_sat = Clapeyron.VT_internal_energy.(model,vl_sat,T)
        uv_sat = Clapeyron.VT_internal_energy.(model,vv_sat,T)
        hl_sat = Clapeyron.VT_enthalpy.(model,vl_sat,T)
        hv_sat = Clapeyron.VT_enthalpy.(model,vv_sat,T)
        sl_sat = Clapeyron.VT_entropy.(model,vl_sat,T)
        sv_sat = Clapeyron.VT_entropy.(model,vv_sat,T)
        cp_l_sat = Clapeyron.VT_isobaric_heat_capacity.(model,vl_sat,T)
        cp_v_sat = Clapeyron.VT_isobaric_heat_capacity.(model,vv_sat,T)
        cv_l_sat = Clapeyron.VT_isochoric_heat_capacity.(model,vl_sat,T)
        cv_v_sat = Clapeyron.VT_isochoric_heat_capacity.(model,vv_sat,T)
        w_l_sat = Clapeyron.VT_speed_of_sound.(model,vl_sat,T)
        w_v_sat = Clapeyron.VT_speed_of_sound.(model,vv_sat,T)
        jt_l_sat = Clapeyron.VT_joule_thomson_coefficient.(model,vl_sat,T)
        jt_v_sat = Clapeyron.VT_joule_thomson_coefficient.(model,vv_sat,T)
        ks_l_sat = Clapeyron.VT_isentropic_compressibility.(model,vl_sat,T)
        ks_v_sat = Clapeyron.VT_isentropic_compressibility.(model,vv_sat,T)
        kt_l_sat = Clapeyron.VT_isothermal_compressibility.(model,vl_sat,T)
        kt_v_sat = Clapeyron.VT_isothermal_compressibility.(model,vv_sat,T)
        β_l_sat = Clapeyron.VT_isobaric_expansivity.(model,vl_sat,T)
        β_v_sat = Clapeyron.VT_isobaric_expansivity.(model,vv_sat,T)

        if x_axis == "Temperature"
            x = vcat(T,reverse(T))
            x_crit = Tc
            x_label = "Temperature  / K"
        elseif x_axis == "Pressure"
            x = vcat(psat,reverse(psat))./1e5
            x_crit = Pc/1e5
            x_label = "Pressure / bar"
        elseif x_axis == "Density"
            x = vcat(ρl_sat,reverse(ρv_sat)).*1e-3
            x_crit = 1/vc*1e-3
            x_label = "Density / (mol/dm³)"
        elseif x_axis == "Volume"
            x = vcat(vl_sat,reverse(vv_sat)).*1e3
            x_crit = vc*1e3
            x_label = "Volume / (dm³/mol)"
        elseif x_axis == "Internal Energy"
            x = vcat(ul_sat,reverse(uv_sat))
            x_crit = uc
            x_label = "Internal Energy / (J/mol)"
        elseif x_axis == "Enthalpy"
            x = vcat(hl_sat,reverse(hv_sat))
            x_crit = hc
            x_label = "Enthalpy / (J/mol)"
        elseif x_axis == "Entropy"
            x = vcat(sl_sat,reverse(sv_sat))
            x_crit = sc
            x_label = "Entropy / (J/mol K)"
        elseif x_axis == "Isobaric Heat Capacity"
            x = vcat(cp_l_sat,reverse(cp_v_sat))
            x_crit = cp_c
            x_label = "Isobaric Heat Capacity / (J/mol K)"
        elseif x_axis == "Isochoric Heat Capacity"
            x = vcat(cv_l_sat,reverse(cv_v_sat))
            x_crit = cv_c
            x_label = "Isochoric Heat Capacity / (J/mol K)"
        elseif x_axis == "Speed of Sound"
            x = vcat(w_l_sat,reverse(w_v_sat))
            x_crit = w_c
            x_label = "Speed of Sound / (m/s)"
        elseif x_axis == "Joule-Thomson Coefficient"
            x = vcat(jt_l_sat,reverse(jt_v_sat)).*1e5
            x_crit = jt_c.*1e5
            x_label = "Joule-Thomson Coefficient / (K/bar)"
        elseif x_axis == "Isentropic Compressibility"
            x = vcat(ks_l_sat,reverse(ks_v_sat))./1e5
            x_crit = ks_c/1e5
            x_label = "Isentropic Compressibility / bar"
        elseif x_axis == "Isothermal Compressibility"
            x = vcat(kt_l_sat,reverse(kt_v_sat))./1e5
            x_crit = kt_c/1e5
            x_label = "Isothermal Compressibility / bar"
        elseif x_axis == "Isobaric Expansivity"
            x = vcat(β_l_sat,reverse(β_v_sat))
            x_crit =β_c
            x_label = "Isobaric Expansivity / K⁻¹"
        end

        if y_axis == "Temperature"
            y = vcat(T,reverse(T))
            y_crit = Tc
            y_label = "Temperature  / K"
        elseif y_axis == "Pressure"
            y = vcat(psat,reverse(psat))./1e5
            y_crit = Pc/1e5
            y_label = "Pressure / bar"
        elseif y_axis == "Density"
            y = vcat(ρl_sat,reverse(ρv_sat)).*1e-3
            y_crit = 1/vc*1e-3
            y_label = "Density / (mol/dm³)"
        elseif y_axis == "Volume"
            y = vcat(vl_sat,reverse(vv_sat)).*1e3
            y_crit = vc*1e3
            y_label = "Volume / (dm³/mol)"
        elseif y_axis == "Internal Energy"
            y = vcat(ul_sat,reverse(uv_sat))
            y_crit = uc
            y_label = "Internal Energy / (J/mol)"
        elseif y_axis == "Enthalpy"
            y = vcat(hl_sat,reverse(hv_sat))
            y_crit = hc
            y_label = "Enthalpy / (J/mol)"
        elseif y_axis == "Entropy"
            y = vcat(sl_sat,reverse(sv_sat))
            y_crit = sc
            y_label = "Entropy / (J/mol K)"
        elseif y_axis == "Isobaric Heat Capacity"
            y = vcat(cp_l_sat,reverse(cp_v_sat))
            y_crit = cp_c
            y_label = "Isobaric Heat Capacity / (J/mol K)"
        elseif y_axis == "Isochoric Heat Capacity"
            y = vcat(cv_l_sat,reverse(cv_v_sat))
            y_crit = cv_c
            y_label = "Isochoric Heat Capacity / (J/mol K)"
        elseif y_axis == "Speed of Sound"
            y = vcat(w_l_sat,reverse(w_v_sat))
            y_crit = w_c
            y_label = "Speed of Sound / (m/s)"
        elseif y_axis == "Joule-Thomson Coefficient"
            y = vcat(jt_l_sat,reverse(jt_v_sat)).*1e5
            y_crit = jt_c
            y_label = "Joule-Thomson Coefficient / (K/bar)"
        elseif y_axis == "Isentropic Compressibility"
            y = vcat(ks_l_sat,reverse(ks_v_sat))./1e5
            y_crit = ks_c/1e5
            y_label = "Isentropic Compressibility / bar"
        elseif y_axis == "Isothermal Compressibility"
            y = vcat(kt_l_sat,reverse(kt_v_sat))./1e5
            y_crit = kt_c/1e5
            y_label = "Isothermal Compressibility / bar"
        elseif y_axis == "Isobaric Expansivity"
            y = vcat(β_l_sat,reverse(β_v_sat))
            y_crit =β_c
            y_label = "Isobaric Expansivity / K⁻¹"
        end

        if log_x
            type_x = "log"
        else
            type_x = "linear"
        end

        if log_y
            type_y = "log"
        else
            type_y = "linear"
        end

        trace = [PlotlyBase.scatter(x=x,y=y,mode="lines",line=PlotlyBase.attr(color="red", dash="solid", width=3),name=""),
                 PlotlyBase.scatter(x=[x_crit],y=[y_crit],mode="markers",marker=PlotlyBase.attr(color="red", size=6))]
        layout = PlotlyBase.Layout(autosize=false,width=700,height=470,
        xaxis = PlotlyBase.attr(title = x_label, font_size=12, showgrid=false,            
                                          ticks="inside",mirror=true,showline=true,linecolor="black",type=type_x),
                             yaxis = PlotlyBase.attr(title = y_label, font_size=12, showgrid=false,       
                                          ticks="inside",mirror=true,showline=true,linecolor="black",type=type_y),
                             showlegend=false, plot_bgcolor="white")
    end
    @onchange x_axis begin
        if x_axis == "Temperature"
            x = vcat(T,reverse(T))
            x_crit = Tc
            x_label = "Temperature  / K"
        elseif x_axis == "Pressure"
            x = vcat(psat,reverse(psat))./1e5
            x_crit = Pc/1e5
            x_label = "Pressure / bar"
        elseif x_axis == "Density"
            x = vcat(ρl_sat,reverse(ρv_sat)).*1e-3
            x_crit = 1/vc*1e-3
            x_label = "Density / (mol/dm³)"
        elseif x_axis == "Volume"
            x = vcat(vl_sat,reverse(vv_sat)).*1e3
            x_crit = vc*1e3
            x_label = "Volume / (dm³/mol)"
        elseif x_axis == "Internal Energy"
            x = vcat(ul_sat,reverse(uv_sat))
            x_crit = uc
            x_label = "Internal Energy / (J/mol)"
        elseif x_axis == "Enthalpy"
            x = vcat(hl_sat,reverse(hv_sat))
            x_crit = hc
            x_label = "Enthalpy / (J/mol)"
        elseif x_axis == "Entropy"
            x = vcat(sl_sat,reverse(sv_sat))
            x_crit = sc
            x_label = "Entropy / (J/mol K)"
        elseif x_axis == "Isobaric Heat Capacity"
            x = vcat(cp_l_sat,reverse(cp_v_sat))
            x_crit = cp_c
            x_label = "Isobaric Heat Capacity / (J/mol K)"
        elseif x_axis == "Isochoric Heat Capacity"
            x = vcat(cv_l_sat,reverse(cv_v_sat))
            x_crit = cv_c
            x_label = "Isochoric Heat Capacity / (J/mol K)"
        elseif x_axis == "Speed of Sound"
            x = vcat(w_l_sat,reverse(w_v_sat))
            x_crit = w_c
            x_label = "Speed of Sound / (m/s)"
        elseif x_axis == "Joule-Thomson Coefficient"
            x = vcat(jt_l_sat,reverse(jt_v_sat)).*1e5
            x_crit = jt_c
            x_label = "Joule-Thomson Coefficient / (K/bar)"
        elseif x_axis == "Isentropic Compressibility"
            x = vcat(ks_l_sat,reverse(ks_v_sat))./1e5
            x_crit = ks_c/1e5
            x_label = "Isentropic Compressibility / bar"
        elseif x_axis == "Isothermal Compressibility"
            x = vcat(kt_l_sat,reverse(kt_v_sat))./1e5
            x_crit = kt_c/1e5
            x_label = "Isothermal Compressibility / bar"
        elseif x_axis == "Isobaric Expansivity"
            x = vcat(β_l_sat,reverse(β_v_sat))
            x_crit =β_c
            x_label = "Isobaric Expansivity / K⁻¹"
        end
        if log_x
            type_x = "log"
        else
            type_x = "linear"
        end

        if log_y
            type_y = "log"
        else
            type_y = "linear"
        end

        trace = [PlotlyBase.scatter(x=x,y=y,mode="lines",line=PlotlyBase.attr(color="red", dash="solid", width=3),name=""),
                 PlotlyBase.scatter(x=[x_crit],y=[y_crit],mode="markers",marker=PlotlyBase.attr(color="red", size=6))]
        layout = PlotlyBase.Layout(autosize=false,width=700,height=470,
        xaxis = PlotlyBase.attr(title = x_label, font_size=12, showgrid=false,            
                                          ticks="inside",mirror=true,showline=true,linecolor="black",type=type_x),
                             yaxis = PlotlyBase.attr(title = y_label, font_size=12, showgrid=false,       
                                          ticks="inside",mirror=true,showline=true,linecolor="black",type=type_y),
                             showlegend=false, plot_bgcolor="white")
    end
    @onchange y_axis begin
        if y_axis == "Temperature"
            y = vcat(T,reverse(T))
            y_crit = Tc
            y_label = "Temperature  / K"
        elseif y_axis == "Pressure"
            y = vcat(psat,reverse(psat))./1e5
            y_crit = Pc/1e5
            y_label = "Pressure / bar"
        elseif y_axis == "Density"
            y = vcat(ρl_sat,reverse(ρv_sat)).*1e-3
            y_crit = 1/vc*1e-3
            y_label = "Density / (mol/dm³)"
        elseif y_axis == "Volume"
            y = vcat(vl_sat,reverse(vv_sat)).*1e3
            y_crit = vc*1e3
            y_label = "Volume / (dm³/mol)"
        elseif y_axis == "Internal Energy"
            y = vcat(ul_sat,reverse(uv_sat))
            y_crit = uc
            y_label = "Internal Energy / (J/mol)"
        elseif y_axis == "Enthalpy"
            y = vcat(hl_sat,reverse(hv_sat))
            y_crit = hc
            y_label = "Enthalpy / (J/mol)"
        elseif y_axis == "Entropy"
            y = vcat(sl_sat,reverse(sv_sat))
            y_crit = sc
            y_label = "Entropy / (J/mol K)"
        elseif y_axis == "Isobaric Heat Capacity"
            y = vcat(cp_l_sat,reverse(cp_v_sat))
            y_crit = cp_c
            y_label = "Isobaric Heat Capacity / (J/mol K)"
        elseif y_axis == "Isochoric Heat Capacity"
            y = vcat(cv_l_sat,reverse(cv_v_sat))
            y_crit = cv_c
            y_label = "Isochoric Heat Capacity / (J/mol K)"
        elseif y_axis == "Speed of Sound"
            y = vcat(w_l_sat,reverse(w_v_sat))
            y_crit = w_c
            y_label = "Speed of Sound / (m/s)"
        elseif y_axis == "Joule-Thomson Coefficient"
            y = vcat(jt_l_sat,reverse(jt_v_sat)).*1e5
            y_crit = jt_c
            y_label = "Joule-Thomson Coefficient / (K/bar)"
        elseif y_axis == "Isentropic Compressibility"
            y = vcat(ks_l_sat,reverse(ks_v_sat))./1e5
            y_crit = ks_c/1e5
            y_label = "Isentropic Compressibility / bar"
        elseif y_axis == "Isothermal Compressibility"
            y = vcat(kt_l_sat,reverse(kt_v_sat))./1e5
            y_crit = kt_c/1e5
            y_label = "Isothermal Compressibility / bar"
        elseif y_axis == "Isobaric Expansivity"
            y = vcat(β_l_sat,reverse(β_v_sat))
            y_crit = β_c
            y_label = "Isobaric Expansivity / K⁻¹"
        end

        if log_x
            type_x = "log"
        else
            type_x = "linear"
        end

        if log_y
            type_y = "log"
        else
            type_y = "linear"
        end

        trace = [PlotlyBase.scatter(x=x,y=y,mode="lines",line=PlotlyBase.attr(color="red", dash="solid", width=3),name=""),
                 PlotlyBase.scatter(x=[x_crit],y=[y_crit],mode="markers",marker=PlotlyBase.attr(color="red", size=6))]
        layout = PlotlyBase.Layout(autosize=false,width=700,height=470,
        xaxis = PlotlyBase.attr(title = x_label, font_size=12, showgrid=false,            
                                          ticks="inside",mirror=true,showline=true,linecolor="black",type=type_x),
                             yaxis = PlotlyBase.attr(title = y_label, font_size=12, showgrid=false,       
                                          ticks="inside",mirror=true,showline=true,linecolor="black",type=type_y),
                             showlegend=false, plot_bgcolor="white")
    end

    @onchange log_x, log_y begin
        if log_x
            type_x = "log"
        else
            type_x = "linear"
        end

        if log_y
            type_y = "log"
        else
            type_y = "linear"
        end

        layout = PlotlyBase.Layout(autosize=false,width=700,height=470,
        xaxis = PlotlyBase.attr(title = x_label, font_size=12, showgrid=false,            
                                          ticks="inside",mirror=true,showline=true,linecolor="black",type=type_x),
                             yaxis = PlotlyBase.attr(title = y_label, font_size=12, showgrid=false,       
                                          ticks="inside",mirror=true,showline=true,linecolor="black",type=type_y),
                             showlegend=false, plot_bgcolor="white")

        trace = deepcopy(trace)
    end
end
end