include("layout.jl")
Html.div(style="display: flex; align-items: center", [imageview(
        src = "https://raw.githubusercontent.com/ClapeyronThermo/Clapeyron.jl/master/docs/Clapeyron_logo.svg",
        height = "auto",
        width = "600px"
    )])
a(href="$(Router.link_to(:get_pure))",
btn("Example", style = "background-color: #4063D8; color: #ffffff",class = "q-mr-sm"))