module App
using GenieFramework
@genietools
include("layout.jl")
include("home.jl")
include("research.jl")
include("people.jl")
include("examples/bulk.jl")
include("examples/saturation.jl")
include("examples/excess.jl")
include("examples/mixpxy.jl")
include("examples/ternary.jl")
include("examples/api_sol.jl")
include("examples/uncertainty.jl")


@page("/home", "home_ui.jl", model = HOME)
@page("/research", "research_ui.jl", model = Research)
@page("/people", "people_ui.jl", model = People)
@page("/bulk_properties", "examples/bulk_ui.jl", model = BULK)
@page("/saturation_properties", "examples/saturation_ui.jl", model = SATURATION)
@page("/excess_properties", "examples/excess_ui.jl", model = EXCESS)
@page("/binary_phase_diagrams", "examples/mixpxy_ui.jl", model = MIXPXY)
@page("/ternary_phase_diagrams", "examples/ternary_ui.jl", model = TERNARY_DIAGRAM)
@page("/api_solubility", "examples/api_sol_ui.jl", model = API_SOL)
@page("/uncertainty_propagation", "examples/uncertainty_ui.jl", model = UNCERTAINTY)
route("/") do
  redirect(:get_home)
end
end