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
include("examples/polymer.jl")
include("examples/uncertainty.jl")
include("examples/parameter_estimation.jl")
include("examples/diffusion.jl")
include("examples/cdft.jl")

meta = Dict("og:title" => "Clapeyron", "og:description" => "A versatile, open-source thermodynamics toolkit for all your modelling needs.", "og:image" => "public/logo.ico")
layout = DEFAULT_LAYOUT(meta=meta)

@page("/home", "home_ui.jl", model = HOME, layout = layout)
@page("/research", "research_ui.jl", model = Research, layout = layout)
@page("/people", "people_ui.jl", model = People, layout = layout)
@page("/bulk_properties", "examples/bulk_ui.jl", model = BULK, layout = layout)
@page("/saturation_properties", "examples/saturation_ui.jl", model = SATURATION, layout = layout)
@page("/excess_properties", "examples/excess_ui.jl", model = EXCESS, layout = layout)
@page("/binary_phase_diagrams", "examples/mixpxy_ui.jl", model = MIXPXY, layout = layout)
@page("/ternary_phase_diagrams", "examples/ternary_ui.jl", model = TERNARY_DIAGRAM, layout = layout)
@page("/api_solubility", "examples/api_sol_ui.jl", model = API_SOL, layout = layout)
@page("/polymer_properties", "examples/polymer_ui.jl", model = POLYMER, layout = layout)
@page("/uncertainty_propagation", "examples/uncertainty_ui.jl", model = UNCERTAINTY, layout = layout)
@page("/parameter_estimation", "examples/parameter_estimation_ui.jl", model = PARAMETER_ESTIMATION, layout = layout)
@page("/diffusion", "examples/diffusion_ui.jl", model = DIFFUSION, layout = layout)
@page("/cdft", "examples/cdft_ui.jl", model = CDFT, layout = layout)
route("/") do
  redirect(:get_home)
end
end