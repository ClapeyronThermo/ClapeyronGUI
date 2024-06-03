macro timeout(seconds, expr_to_run, expr_when_fails="Timed out")
    quote
        tsk = @task $(esc(expr_to_run))
        schedule(tsk)
        Timer($(esc(seconds))) do timer
            istaskdone(tsk) || Base.throwto(tsk, InterruptException())
        end
        # try
            fetch(tsk)
        # catch _
        #     $(esc(expr_when_fails))
        # end
    end
end