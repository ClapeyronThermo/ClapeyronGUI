module App
using GenieFramework
@genietools
include("home.jl")
include("purept.jl")
include("purerhot.jl")
include("mixpxy.jl")

@page("/home", "home_ui.jl", layout = "layout.jl", model = HOME)
@page("/purept", "purept_ui.jl", layout = "layout.jl", model = PUREPT)
@page("/purerhot", "purerhot_ui.jl", layout = "layout.jl", model = PURERHOT)
@page("/mixpxy", "mixpxy_ui.jl", layout = "layout.jl", model = MIXPXY)
route("/") do
    redirect(:get_home)
end

end