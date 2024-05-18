include("layout.jl")
@eval page_layout( [
row([
    cell(
        class="st-col", style="display:block;margin-left:auto;margin-right:auto;align-items:center", [imageview(
        src = "https://raw.githubusercontent.com/ClapeyronThermo/Clapeyron.jl/master/docs/Clapeyron_logo.svg",
        height = "auto",
        width = "600px"
    )])])
])
