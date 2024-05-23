module PARAMETER_ESTIMATION
using GenieFramework
using Clapeyron, Main.ThermoPlots
using Metaheuristics, DataFrames, CSV
import PlotlyBase, PlotlyJS, PlotlyKaleido
@genietools

@app begin
    @in tab_selected = "pT"

    @in species = "carbon dioxide"
    
    @in fit_button = false
    @in data_button = false 

    @out epsilon = 0.
    @out sigma = 0.
    @out segment = 0.

    @out objective_fun = []
    @out iter = []
    @out T_exp = []
    @out exp_sat_p = []
    @out exp_sat_rho = []
    @out fit_sat_p = []
    @out fit_sat_rho = []

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
    @out trace_obj = []
    @out layout_obj = PlotlyBase.Layout(autosize=false,width=700,height=470,
    xaxis = PlotlyBase.attr(title = "Iteration", font_size=12, showgrid=false,            
                                      ticks="inside",mirror=true,showline=true,linecolor="black",type="linear"),
                         yaxis = PlotlyBase.attr(title = "Objective function", font_size=12, showgrid=false,       
                                      ticks="inside",mirror=true,showline=true,linecolor="black",type="log"),
                         showlegend=false, plot_bgcolor="white")

    @onbutton data_button begin
        Npoints = 20
        model = MultiFluid([species])
        
        (Tc,Pc,vc) = crit_pure(model)
        T_exp = LinRange(0.7*Tc,Tc*0.95,Npoints)

        sat = saturation_pressure.(model,T_exp)
        exp_sat_p = [sat[i][1] for i in 1:Npoints]
        exp_sat_rho = [1e-3 /sat[i][2] for i in 1:Npoints]

        a=Any["[method=saturation_p_rhol]","T"]
        append!(a,T_exp)
        b=Any["","out_pre"]
        append!(b,exp_sat_p)
        c=Any["","out_rho"]
        append!(c,exp_sat_rho)

        df_exp = DataFrame(A=a, B=b, C=c)

        CSV.write("exp_data.csv",df_exp)

        trace = [PlotlyBase.scatter(x=T_exp,y=exp_sat_p./1e5,mode="markers",marker=PlotlyBase.attr(color="blue", size=6),name="")]

        trace_rho = [PlotlyBase.scatter(x=exp_sat_rho,y=T_exp,mode="markers",marker=PlotlyBase.attr(color="blue", size=6),name="")]
    end

    @onbutton fit_button begin
        Npoints = 200
        model = PCSAFT([species])
        toestimate = [
                    Dict(
                        :param => :epsilon,
                        :lower => model.params.epsilon[1]*0.5,
                        :upper => model.params.epsilon[1]*1.5,
                        :guess => model.params.epsilon[1]
                    ),
                    Dict(
                        :param => :sigma,
                        :factor => 1e-10,
                        :lower => model.params.sigma[1]/1e-10*0.5,
                        :upper => model.params.sigma[1]/1e-10*1.5,
                        :guess => model.params.sigma[1]/1e-10
                    )
                    ,
                    Dict(
                        :param => :segment,
                        :lower => model.params.segment[1]*0.5,
                        :upper => model.params.segment[1]*1.5,
                        :guess => model.params.segment[1]
                    )
                    ]
        

        e,objective,initial,upper,lower = Estimation(model,toestimate,["exp_data.csv"])

        T = LinRange(minimum(T_exp),maximum(T_exp),Npoints)

        function logger_new(st)
            # if st.iteration % 10 == 0
                # println("Iteration: ",st.iteration)
                # println("Best solution: ",st.best_sol.x)
                epsilon = round(st.best_sol.x[1],digits=2)
                sigma = round(st.best_sol.x[2],digits=2)
                segment = round(st.best_sol.x[3],digits=2)
                # println("Best value: ",st.best_sol.f)
                append!(iter,st.iteration)
                append!(objective_fun,st.best_sol.f)

                trace_obj = [PlotlyBase.scatter(x=iter,y=objective_fun,mode="lines",line=PlotlyBase.attr(color="blue", dash="solid", width=3),name="")]

                model_new = return_model(e,model,st.best_sol.x)

                sat = saturation_pressure.(model_new,T)
                fit_sat_p = [sat[i][1] for i in 1:Npoints]
                fit_sat_rho = [1e-3 /sat[i][2] for i in 1:Npoints]

                trace = [PlotlyBase.scatter(x=T_exp,y=exp_sat_p./1e5,mode="markers",marker=PlotlyBase.attr(color="blue", size=6),name=""),
                 PlotlyBase.scatter(x=T,y=fit_sat_p./1e5,mode="lines",line=PlotlyBase.attr(color="blue", dash="solid", width=3),name="")]

                trace_rho = [PlotlyBase.scatter(x=exp_sat_rho,y=T_exp,mode="markers",marker=PlotlyBase.attr(color="blue", size=6),name=""),
                     PlotlyBase.scatter(x=fit_sat_rho,y=T,mode="lines",line=PlotlyBase.attr(color="blue", dash="solid", width=3),name="")]
                
            # end
        end

        params, model_new = optimize(objective,e,ECA(options=Options(iterations=200));verbose=true,logger=logger_new)

        # (Tc,Pc,vc) = crit_pure(model_new)
        

        

        epsilon = round(params[1],digits=2)
        sigma = round(params[2],digits=2)
        segment = round(params[3],digits=2)

        iter = []
        objective_fun = []

        trace = [PlotlyBase.scatter(x=T_exp,y=exp_sat_p./1e5,mode="markers",marker=PlotlyBase.attr(color="blue", size=6),name=""),
                 PlotlyBase.scatter(x=T,y=fit_sat_p./1e5,mode="lines",line=PlotlyBase.attr(color="blue", dash="solid", width=3),name="")]

        trace_rho = [PlotlyBase.scatter(x=exp_sat_rho,y=T_exp,mode="markers",marker=PlotlyBase.attr(color="blue", size=6),name=""),
                     PlotlyBase.scatter(x=fit_sat_rho,y=T,mode="lines",line=PlotlyBase.attr(color="blue", dash="solid", width=3),name="")]

    end

end
end