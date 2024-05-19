module App
using GenieFramework
@genietools
include("layout.jl")
include("home.jl")
include("research.jl")
include("people.jl")
include("examples/bulk.jl")
include("examples/saturation.jl")
include("examples/mixpxy.jl")
include("examples/api_sol.jl")


@page("/home", "home_ui.jl", model = HOME)
@page("/research", "research_ui.jl", model = Research)
@page("/people", "people_ui.jl", model = People)
@page("/bulk_properties", "examples/bulk_ui.jl", model = BULK)
@page("/saturation_properties", "examples/saturation_ui.jl", model = SATURATION)
@page("/pxy_diagram", "examples/mixpxy_ui.jl", model = MIXPXY)
@page("/api_solubility", "examples/api_sol_ui.jl", model = API_SOL)
route("/") do
  redirect(:get_home)
end
end