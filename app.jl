module App
using GenieFramework
@genietools
include("layout.jl")
include("home.jl")
include("research.jl")
include("people.jl")
include("examples/purept.jl")
include("examples/purerhot.jl")
include("examples/saturation.jl")
include("examples/mixpxy.jl")


@page("/home", "home_ui.jl", model = HOME)
@page("/research", "research_ui.jl", model = Research)
@page("/people", "people_ui.jl", model = People)
@page("/purept", "examples/purept_ui.jl", model = PUREPT)
@page("/purerhot", "examples/purerhot_ui.jl", model = PURERHOT)
@page("/saturation", "examples/saturation_ui.jl", model = SATURATION)
@page("/mixpxy", "examples/mixpxy_ui.jl", model = MIXPXY)
route("/") do
  redirect(:get_home)
end
end