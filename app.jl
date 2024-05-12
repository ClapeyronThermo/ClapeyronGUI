module App
using GenieFramework
@genietools
include("layout.jl")
include("home.jl")
include("research.jl")
include("people.jl")
include("examples/purept.jl")
include("examples/purerhot.jl")
include("examples/mixpxy.jl")

@page("/home", "home_ui.jl", model = HOME)
@page("/research", "research_ui.jl", model = Research)
@page("/people", "people_ui.jl", model = People)
@page("/purept", "examples/purept_ui.jl", model = PUREPT)
@page("/purerhot", "examples/purerhot_ui.jl", model = PURERHOT)
@page("/mixpxy", "examples/mixpxy_ui.jl", model = MIXPXY)
route("/") do
  redirect(:get_home)
end
end