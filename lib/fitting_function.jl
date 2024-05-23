import Clapeyron
function saturation_p_rhol(model,T)
    sat = Clapeyron.saturation_pressure(model,T)
    return sat[1], 1e-3/sat[2]
end

function logger_new(st)
    if st.iteration % 10 == 0
        println("Iteration: ",st.iteration)
        println("Best solution: ",st.best_sol.x)
        println("Best value: ",st.best_sol.f)
    end
end