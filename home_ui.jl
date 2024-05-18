include("layout.jl")
@eval page_layout( [
row([cell(class="st-col col-2"),
    cell(
        class="st-col col-8", style="display:block;margin-left:auto;margin-right:auto;align-items:center", [imageview(
        src = "https://raw.githubusercontent.com/ClapeyronThermo/Clapeyron.jl/master/docs/Clapeyron_logo.svg",
        height = "auto",
        width = "800px"
    )]),
    cell(class="st-col col-2")])
])
