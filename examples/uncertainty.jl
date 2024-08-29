module UNCERTAINTY
using GenieFramework
using Clapeyron, Main.ThermoPlots
using Measurements, ForwardDiffOverMeasurements
import PlotlyBase, PlotlyKaleido
@genietools
@app begin
    @in tab_selected = "pT"

    @in Slider_epsilon = 0.
    @in Slider_sigma = 0.
    @in Slider_segment = 0.
    @in Slider_pre = 1.

    @in species = "carbon dioxide"
    @in Selected_property = "Density"
    @in new_button = false 

    @out Select_property = ["Density","Volume",
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
    @out trace_rho = []
    @out layout_rho = PlotlyBase.Layout(autosize=false,width=700,height=470,
    xaxis = PlotlyBase.attr(title = "Density / (mol/dm³)", font_size=12, showgrid=false,            
                                      ticks="inside",mirror=true,showline=true,linecolor="black",type="linear"),
                         yaxis = PlotlyBase.attr(title = "Temperature  / K", font_size=12, showgrid=false,       
                                      ticks="inside",mirror=true,showline=true,linecolor="black"),
                         showlegend=false, plot_bgcolor="white")
    @out trace_bulk = []
    @out layout_bulk = PlotlyBase.Layout(autosize=false,width=700,height=470,
    xaxis = PlotlyBase.attr(title = "Temperature  / K", font_size=12, showgrid=false,            
                                      ticks="inside",mirror=true,showline=true,linecolor="black",type="linear"),
                         yaxis = PlotlyBase.attr(title = "Density / (mol/dm³)", font_size=12, showgrid=false,       
                                      ticks="inside",mirror=true,showline=true,linecolor="black"),
                         showlegend=false, plot_bgcolor="white")

    @onbutton new_button begin
        Npoints = 200

        try
            model = PCSAFT([species])
        catch
            notify(__model__, "Species $species is not available in PCSAFT.", :warning)
            throw(TypeError("Species $species is not available in PCSAFT."))
        end

        model = PCSAFT([species])

        model_err = PCSAFT([species];userlocations = (;
                            Mw = [model.params.Mw[1]±0],
                            epsilon = [model.params.epsilon[1]±Slider_epsilon],
                            sigma = [model.params.sigma[1]/1e-10±Slider_sigma],
                            segment = [model.params.segment[1]±Slider_segment],
                            k = [0.0;;], #matrix
                            n_H = [0],
                            n_e = [0],
                            epsilon_assoc = nothing,
                            bondvol = nothing))

        (Tc,Pc,vc) = crit_pure(model)
        ρc = 1/vc

        T = LinRange(0.6*Tc,Tc,Npoints)
        sat = saturation_pressure.(model_err,T)

        psat = [sat[i][1].val for i in 1:Npoints]
        psat_err = [sat[i][1].err for i in 1:Npoints]

        vl_sat = [sat[i][2] for i in 1:Npoints]
        vv_sat = [sat[i][3] for i in 1:Npoints]
        T = T[.!isnan.(psat)]
        vl_sat = vl_sat[.!isnan.(psat)]
        vv_sat = vv_sat[.!isnan.(psat)]
        psat_err = psat_err[.!isnan.(psat)]
        psat = psat[.!isnan.(psat)]

        ρl_sat = [1/vl for vl in vl_sat]
        ρv_sat = [1/vv for vv in vv_sat]
        ρl_sat_err = [ρl_sat[i].err for i in 1:length(psat)]
        ρv_sat_err = [ρv_sat[i].err for i in 1:length(psat)]
        ρl_sat = [ρl_sat[i].val for i in 1:length(psat)]
        ρv_sat = [ρv_sat[i].val for i in 1:length(psat)]

        trace = [PlotlyBase.scatter(x=T,y=(psat.-psat_err)./1e5,mode="lines",line=PlotlyBase.attr(color="blue", dash="solid", width=0),name=""),
            PlotlyBase.scatter(x=T,y=psat./1e5,mode="lines",
            fill="tonexty",fillcolor="rgba(0,0,255,0.5)",line=PlotlyBase.attr(color="blue", dash="solid", width=3),name=""),
            PlotlyBase.scatter(x=T,y=(psat.+psat_err)./1e5,mode="lines",fill="tonexty",fillcolor="rgba(0,0,255,0.5)",line=PlotlyBase.attr(color="blue", dash="solid", width=0),name=""),
                 PlotlyBase.scatter(x=[Tc],y=[Pc/1e5],mode="markers",marker=PlotlyBase.attr(color="blue", size=6))]
        trace_rho = [PlotlyBase.scatter(y=T,x=(ρl_sat.-ρl_sat_err).*1e-3,mode="lines",line=PlotlyBase.attr(color="blue", dash="solid", width=0),name=""),
            PlotlyBase.scatter(y=T,x=ρl_sat.*1e-3,mode="lines",fill="tonextx",fillcolor="rgba(0,0,255,0.5)",line=PlotlyBase.attr(color="blue", dash="solid", width=3),name=""),
            PlotlyBase.scatter(y=T,x=(ρl_sat.+ρl_sat_err).*1e-3,mode="lines",fill="tonextx",fillcolor="rgba(0,0,255,0.5)",line=PlotlyBase.attr(color="blue", dash="solid", width=0),name=""),
            PlotlyBase.scatter(y=T,x=(ρv_sat.+ρv_sat_err).*1e-3,mode="lines",line=PlotlyBase.attr(color="blue", dash="solid", width=0),name=""),
            PlotlyBase.scatter(y=T,x=ρv_sat.*1e-3,mode="lines",fill="tonextx",fillcolor="rgba(0,0,255,0.5)",line=PlotlyBase.attr(color="blue", dash="solid", width=3),name=""),
            PlotlyBase.scatter(y=T,x=(ρv_sat.-ρv_sat_err).*1e-3,mode="lines",fill="tonextx",fillcolor="rgba(0,0,255,0.5)",line=PlotlyBase.attr(color="blue", dash="solid", width=0),name=""),
            PlotlyBase.scatter(y=[Tc],x=[ρc*1e-3],mode="markers",marker=PlotlyBase.attr(color="blue", size=6))]

            pre = Slider_pre*1e5
            T = LinRange(250,500,200)

            if Selected_property == "Volume"
                y = volume.(model_err,pre,T)
                y_err = [y[i].err for i in 1:200].*1e3
                y = [y[i].val for i in 1:200].*1e3
                y_label = "Volume / (dm³/mol)"
            elseif Selected_property == "Density"
                y = 1e-3./volume.(model_err,pre,T)
                y_err = [y[i].err for i in 1:200]
                y = [y[i].val for i in 1:200]
                y_label = "Density / (mol/dm³)"
            elseif Selected_property == "Internal Energy"
                y = internal_energy.(model_err,pre,T)
                y_err = [y[i].err for i in 1:200]
                y = [y[i].val for i in 1:200]
                y_label = "Internal Energy / (J/mol)"
            elseif Selected_property == "Enthalpy"
                y = enthalpy.(model_err,pre,T)
                y_err = [y[i].err for i in 1:200]
                y = [y[i].val for i in 1:200]
                y_label = "Enthalpy / (J/mol)"
            elseif Selected_property == "Entropy"
                y = entropy.(model_err,pre,T)
                y_err = [y[i].err for i in 1:200]
                y = [y[i].val for i in 1:200]
                y_label = "Entropy / (J/mol K)"
            elseif Selected_property == "Isobaric Heat Capacity"
                y = isobaric_heat_capacity.(model_err,pre,T)
                y_err = [y[i].err for i in 1:200]
                y = [y[i].val for i in 1:200]
                y_label = "Isobaric Heat Capacity / (J/mol K)"
            elseif Selected_property == "Isochoric Heat Capacity"
                y = isochoric_heat_capacity.(model_err,pre,T)
                y_err = [y[i].err for i in 1:200]
                y = [y[i].val for i in 1:200]
                y_label = "Isochoric Heat Capacity / (J/mol K)"
            elseif Selected_property == "Speed of Sound"
                y = speed_of_sound.(model_err,pre,T)
                y_err = [y[i].err for i in 1:200]
                y = [y[i].val for i in 1:200]
                y_label = "Speed of Sound / (m/s)"
            elseif Selected_property == "Joule-Thomson Coefficient"
                y = joule_thomson_coefficient.(model_err,pre,T)
                y_err = [y[i].err for i in 1:200]
                y = [y[i].val for i in 1:200]
                y_label = "Joule-Thomson Coefficient / (K/bar)"
            elseif Selected_property == "Isentropic Compressibility"
                y = isentropic_compressibility.(model_err,pre,T)
                y_err = [y[i].err for i in 1:200]
                y = [y[i].val for i in 1:200]
                y_label = "Isentropic Compressibility / (1/bar)"
            elseif Selected_property == "Isothermal Compressibility"
                y = isothermal_compressibility.(model_err,pre,T)
                y_err = [y[i].err for i in 1:200]
                y = [y[i].val for i in 1:200]
                y_label = "Isothermal Compressibility / (1/bar)"
            elseif Selected_property == "Isobaric Expansivity"
                y = isobaric_expansivity.(model_err,pre,T)
                y_err = [y[i].err for i in 1:200]
                y = [y[i].val for i in 1:200]
                y_label = "Isobaric Expansivity / (1/K)"
            end

            trace_bulk = [PlotlyBase.scatter(x=T,y=(y.-y_err),mode="lines",line=PlotlyBase.attr(color="blue", dash="solid", width=0),name=""),
                PlotlyBase.scatter(x=T,y=y,mode="lines",fill="tonexty",fillcolor="rgba(0,0,255,0.5)",line=PlotlyBase.attr(color="blue", dash="solid", width=3),name=""),
                PlotlyBase.scatter(x=T,y=(y.+y_err),mode="lines",fill="tonexty",fillcolor="rgba(0,0,255,0.5)",line=PlotlyBase.attr(color="blue", dash="solid", width=0),name="")]
    end

    @onchange Slider_sigma, Slider_epsilon, Slider_segment begin
        Npoints = 200
        try
            model = PCSAFT([species])
        catch
            notify(__model__, "Species $species is not available in PCSAFT.", :warning)
            throw(TypeError("Species $species is not available in PCSAFT."))
        end

        model = PCSAFT([species])

        model_err = PCSAFT([species];userlocations = (;
                            Mw = [model.params.Mw[1]±0],
                            epsilon = [model.params.epsilon[1]±Slider_epsilon],
                            sigma = [model.params.sigma[1]/1e-10±Slider_sigma],
                            segment = [model.params.segment[1]±Slider_segment],
                            k = [0.0;;], #matrix
                            n_H = [0],
                            n_e = [0],
                            epsilon_assoc = nothing,
                            bondvol = nothing))

        (Tc,Pc,vc) = crit_pure(model)
        ρc = 1/vc

        T = LinRange(0.6*Tc,Tc,Npoints)
        sat = saturation_pressure.(model_err,T)

        psat = [sat[i][1].val for i in 1:Npoints]
        psat_err = [sat[i][1].err for i in 1:Npoints]

        vl_sat = [sat[i][2] for i in 1:Npoints]
        vv_sat = [sat[i][3] for i in 1:Npoints]
        T = T[.!isnan.(psat)]
        vl_sat = vl_sat[.!isnan.(psat)]
        vv_sat = vv_sat[.!isnan.(psat)]
        psat_err = psat_err[.!isnan.(psat)]
        psat = psat[.!isnan.(psat)]

        ρl_sat = [1/vl for vl in vl_sat]
        ρv_sat = [1/vv for vv in vv_sat]
        ρl_sat_err = [ρl_sat[i].err for i in 1:length(psat)]
        ρv_sat_err = [ρv_sat[i].err for i in 1:length(psat)]
        ρl_sat = [ρl_sat[i].val for i in 1:length(psat)]
        ρv_sat = [ρv_sat[i].val for i in 1:length(psat)]

        trace = [PlotlyBase.scatter(x=T,y=(psat.-psat_err)./1e5,mode="lines",line=PlotlyBase.attr(color="blue", dash="solid", width=0),name=""),
            PlotlyBase.scatter(x=T,y=psat./1e5,mode="lines",fill="tonexty",fillcolor="rgba(0,0,255,0.5)",line=PlotlyBase.attr(color="blue", dash="solid", width=3),name=""),
            PlotlyBase.scatter(x=T,y=(psat.+psat_err)./1e5,mode="lines",fill="tonexty",fillcolor="rgba(0,0,255,0.5)",line=PlotlyBase.attr(color="blue", dash="solid", width=0),name=""),
                 PlotlyBase.scatter(x=[Tc],y=[Pc/1e5],mode="markers",marker=PlotlyBase.attr(color="blue", size=6))]

        trace_rho = [PlotlyBase.scatter(y=T,x=(ρl_sat.-ρl_sat_err).*1e-3,mode="lines",line=PlotlyBase.attr(color="blue", dash="solid", width=0),name=""),
            PlotlyBase.scatter(y=T,x=ρl_sat.*1e-3,mode="lines",fill="tonextx",fillcolor="rgba(0,0,255,0.5)",line=PlotlyBase.attr(color="blue", dash="solid", width=3),name=""),
            PlotlyBase.scatter(y=T,x=(ρl_sat.+ρl_sat_err).*1e-3,mode="lines",fill="tonextx",fillcolor="rgba(0,0,255,0.5)",line=PlotlyBase.attr(color="blue", dash="solid", width=0),name=""),
            PlotlyBase.scatter(y=T,x=(ρv_sat.+ρv_sat_err).*1e-3,mode="lines",line=PlotlyBase.attr(color="blue", dash="solid", width=0),name=""),
            PlotlyBase.scatter(y=T,x=ρv_sat.*1e-3,mode="lines",fill="tonextx",fillcolor="rgba(0,0,255,0.5)",line=PlotlyBase.attr(color="blue", dash="solid", width=3),name=""),
            PlotlyBase.scatter(y=T,x=(ρv_sat.-ρv_sat_err).*1e-3,mode="lines",fill="tonextx",fillcolor="rgba(0,0,255,0.5)",line=PlotlyBase.attr(color="blue", dash="solid", width=0),name=""),
            PlotlyBase.scatter(y=[Tc],x=[ρc*1e-3],mode="markers",marker=PlotlyBase.attr(color="blue", size=6))]
            
            pre = Slider_pre*1e5
            T = LinRange(250,500,200)

            if Selected_property == "Volume"
                y = volume.(model_err,pre,T)
                y_err = [y[i].err for i in 1:200].*1e3
                y = [y[i].val for i in 1:200].*1e3
                y_label = "Volume / (dm³/mol)"
            elseif Selected_property == "Density"
                y = 1e-3./volume.(model_err,pre,T)
                y_err = [y[i].err for i in 1:200]
                y = [y[i].val for i in 1:200]
                y_label = "Density / (mol/dm³)"
            elseif Selected_property == "Internal Energy"
                y = internal_energy.(model_err,pre,T)
                y_err = [y[i].err for i in 1:200]
                y = [y[i].val for i in 1:200]
                y_label = "Internal Energy / (J/mol)"
            elseif Selected_property == "Enthalpy"
                y = enthalpy.(model_err,pre,T)
                y_err = [y[i].err for i in 1:200]
                y = [y[i].val for i in 1:200]
                y_label = "Enthalpy / (J/mol)"
            elseif Selected_property == "Entropy"
                y = entropy.(model_err,pre,T)
                y_err = [y[i].err for i in 1:200]
                y = [y[i].val for i in 1:200]
                y_label = "Entropy / (J/mol K)"
            elseif Selected_property == "Isobaric Heat Capacity"
                y = isobaric_heat_capacity.(model_err,pre,T)
                y_err = [y[i].err for i in 1:200]
                y = [y[i].val for i in 1:200]
                y_label = "Isobaric Heat Capacity / (J/mol K)"
            elseif Selected_property == "Isochoric Heat Capacity"
                y = isochoric_heat_capacity.(model_err,pre,T)
                y_err = [y[i].err for i in 1:200]
                y = [y[i].val for i in 1:200]
                y_label = "Isochoric Heat Capacity / (J/mol K)"
            elseif Selected_property == "Speed of Sound"
                y = speed_of_sound.(model_err,pre,T)
                y_err = [y[i].err for i in 1:200]
                y = [y[i].val for i in 1:200]
                y_label = "Speed of Sound / (m/s)"
            elseif Selected_property == "Joule-Thomson Coefficient"
                y = joule_thomson_coefficient.(model_err,pre,T)
                y_err = [y[i].err for i in 1:200]
                y = [y[i].val for i in 1:200]
                y_label = "Joule-Thomson Coefficient / (K/bar)"
            elseif Selected_property == "Isentropic Compressibility"
                y = isentropic_compressibility.(model_err,pre,T)
                y_err = [y[i].err for i in 1:200]
                y = [y[i].val for i in 1:200]
                y_label = "Isentropic Compressibility / (1/bar)"
            elseif Selected_property == "Isothermal Compressibility"
                y = isothermal_compressibility.(model_err,pre,T)
                y_err = [y[i].err for i in 1:200]
                y = [y[i].val for i in 1:200]
                y_label = "Isothermal Compressibility / (1/bar)"
            elseif Selected_property == "Isobaric Expansivity"
                y = isobaric_expansivity.(model_err,pre,T)
                y_err = [y[i].err for i in 1:200]
                y = [y[i].val for i in 1:200]
                y_label = "Isobaric Expansivity / (1/K)"
            end

            trace_bulk = [PlotlyBase.scatter(x=T,y=(y.-y_err),mode="lines",line=PlotlyBase.attr(color="blue", dash="solid", width=0),name=""),
                PlotlyBase.scatter(x=T,y=y,mode="lines",fill="tonexty",fillcolor="rgba(0,0,255,0.5)",line=PlotlyBase.attr(color="blue", dash="solid", width=3),name=""),
                PlotlyBase.scatter(x=T,y=(y.+y_err),mode="lines",fill="tonexty",fillcolor="rgba(0,0,255,0.5)",line=PlotlyBase.attr(color="blue", dash="solid", width=0),name="")]
            layout_bulk = PlotlyBase.Layout(autosize=false,width=700,height=470,
            xaxis = PlotlyBase.attr(title = "Temperature  / K", font_size=12, showgrid=false,            
                                          ticks="inside",mirror=true,showline=true,linecolor="black",type="linear"),
                                 yaxis = PlotlyBase.attr(title = y_label, font_size=12, showgrid=false,       
                                          ticks="inside",mirror=true,showline=true,linecolor="black"),
                                 showlegend=false, plot_bgcolor="white")
    end

    @onchange Slider_pre, Selected_property begin
        try
            model = PCSAFT([species])
        catch
            notify(__model__, "Species $species is not available in PCSAFT.", :warning)
            throw(TypeError("Species $species is not available in PCSAFT."))
        end

        pre = Slider_pre*1e5

        if pre < 0
            notify(__model__, "Pressure must be positive.", :warning)
            throw(TypeError("Pressure must be positive."))
        end

        model = PCSAFT([species])

        model_err = PCSAFT([species];userlocations = (;
                        Mw = [model.params.Mw[1]±0],
                        epsilon = [model.params.epsilon[1]±Slider_epsilon],
                        sigma = [model.params.sigma[1]/1e-10±Slider_sigma],
                        segment = [model.params.segment[1]±Slider_segment],
                        k = [0.0;;], #matrix
                        n_H = [0],
                        n_e = [0],
                        epsilon_assoc = nothing,
                        bondvol = nothing))

        T = LinRange(250,500,200)
        if Selected_property == "Volume"
            y = volume.(model_err,pre,T)
            y_err = [y[i].err for i in 1:200].*1e3
            y = [y[i].val for i in 1:200].*1e3
            y_label = "Volume / (dm³/mol)"
        elseif Selected_property == "Density"
            y = 1e-3./volume.(model_err,pre,T)
            y_err = [y[i].err for i in 1:200]
            y = [y[i].val for i in 1:200]
            y_label = "Density / (mol/dm³)"
        elseif Selected_property == "Internal Energy"
            y = internal_energy.(model_err,pre,T)
            y_err = [y[i].err for i in 1:200]
            y = [y[i].val for i in 1:200]
            y_label = "Internal Energy / (J/mol)"
        elseif Selected_property == "Enthalpy"
            y = enthalpy.(model_err,pre,T)
            y_err = [y[i].err for i in 1:200]
            y = [y[i].val for i in 1:200]
            y_label = "Enthalpy / (J/mol)"
        elseif Selected_property == "Entropy"
            y = entropy.(model_err,pre,T)
            y_err = [y[i].err for i in 1:200]
            y = [y[i].val for i in 1:200]
            y_label = "Entropy / (J/mol K)"
        elseif Selected_property == "Isobaric Heat Capacity"
            y = isobaric_heat_capacity.(model_err,pre,T)
            y_err = [y[i].err for i in 1:200]
            y = [y[i].val for i in 1:200]
            y_label = "Isobaric Heat Capacity / (J/mol K)"
        elseif Selected_property == "Isochoric Heat Capacity"
            y = isochoric_heat_capacity.(model_err,pre,T)
            y_err = [y[i].err for i in 1:200]
            y = [y[i].val for i in 1:200]
            y_label = "Isochoric Heat Capacity / (J/mol K)"
        elseif Selected_property == "Speed of Sound"
            y = speed_of_sound.(model_err,pre,T)
            y_err = [y[i].err for i in 1:200]
            y = [y[i].val for i in 1:200]
            y_label = "Speed of Sound / (m/s)"
        elseif Selected_property == "Joule-Thomson Coefficient"
            y = joule_thomson_coefficient.(model_err,pre,T)
            y_err = [y[i].err for i in 1:200]
            y = [y[i].val for i in 1:200]
            y_label = "Joule-Thomson Coefficient / (K/bar)"
        elseif Selected_property == "Isentropic Compressibility"
            y = isentropic_compressibility.(model_err,pre,T)
            y_err = [y[i].err for i in 1:200]
            y = [y[i].val for i in 1:200]
            y_label = "Isentropic Compressibility / (1/bar)"
        elseif Selected_property == "Isothermal Compressibility"
            y = isothermal_compressibility.(model_err,pre,T)
            y_err = [y[i].err for i in 1:200]
            y = [y[i].val for i in 1:200]
            y_label = "Isothermal Compressibility / (1/bar)"
        elseif Selected_property == "Isobaric Expansivity"
            y = isobaric_expansivity.(model_err,pre,T)
            y_err = [y[i].err for i in 1:200]
            y = [y[i].val for i in 1:200]
            y_label = "Isobaric Expansivity / (1/K)"
        end

        trace_bulk = [PlotlyBase.scatter(x=T,y=(y.-y_err),mode="lines",line=PlotlyBase.attr(color="blue", dash="solid", width=0),name=""),
            PlotlyBase.scatter(x=T,y=y,mode="lines",fill="tonexty",fillcolor="rgba(0,0,255,0.5)",line=PlotlyBase.attr(color="blue", dash="solid", width=3),name=""),
            PlotlyBase.scatter(x=T,y=(y.+y_err),mode="lines",fill="tonexty",fillcolor="rgba(0,0,255,0.5)",line=PlotlyBase.attr(color="blue", dash="solid", width=0),name="")]
        layout_bulk = PlotlyBase.Layout(autosize=false,width=700,height=470,
        xaxis = PlotlyBase.attr(title = "Temperature  / K", font_size=12, showgrid=false,            
                                      ticks="inside",mirror=true,showline=true,linecolor="black",type="linear"),
                             yaxis = PlotlyBase.attr(title = y_label, font_size=12, showgrid=false,       
                                      ticks="inside",mirror=true,showline=true,linecolor="black"),
                             showlegend=false, plot_bgcolor="white")
    end
end
end