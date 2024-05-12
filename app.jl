module App
using GenieFramework
@genietools
include("home.jl")
include("examples/purept.jl")
include("examples/purerhot.jl")
include("examples/mixpxy.jl")

@page("/home", "home_ui.jl", layout = "layout.jl", model = HOME)
@page("/purept", "examples/purept_ui.jl", layout = "layout.jl", model = PUREPT)
@page("/purerhot", "examples/purerhot_ui.jl", layout = "layout.jl", model = PURERHOT)
@page("/mixpxy", "examples/mixpxy_ui.jl", layout = "layout.jl", model = MIXPXY)
route("/") do
    redirect(:get_home)
end

end