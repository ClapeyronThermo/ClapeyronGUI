module App
using GenieFramework
@genietools
include("home.jl")
include("pure.jl")

@page("/home", "home_ui.jl", layout = "layout.jl", model = HOME)
@page("/pure", "pure_ui.jl", layout = "layout.jl", model = PURE)
route("/") do
    redirect(:get_home)
end

end