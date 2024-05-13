module ThermoPlots

using Clapeyron, PlotlyJS

# include("bulk/bulk.jl")

include("pure/pT_curve.jl")
include("pure/rhoT_envelope.jl")
include("pure/rhop_envelope.jl")

include("binary/pxy_diagram.jl")
# include("binary/Txy_diagram.jl")
end