
Html.div(style="display: flex; align-items: center", [imageview(
        src = "https://raw.githubusercontent.com/ClapeyronThermo/Clapeyron.jl/master/docs/Clapeyron_logo.svg",
        height = "auto",
        width = "600px"
    )])
a(href="$(Router.link_to(:get_pure))",
btn("Example", style = "background-color: #4063D8; color: #ffffff",class = "q-mr-sm"))
Html.div(style="display: flex; gap: 20px; font-size: 1.5em;",
        expansionitem(
            label = "Examples",
            dense = true,
            var"dense-toggle" = true,
            var"expand-separator" = true,
            var"header-class" = "bg-blue-1",
            list(bordered = true,
            class = "rounded-borders",
            [item(clickable = true, a(href="$(Router.link_to(:get_purept))", style="text-decoration: none;; color: #000000;","Saturation curve")),
            item(clickable = true, a(href="$(Router.link_to(:get_purerhot))", style="text-decoration: none;; color: #000000;","VLE envelope")),
            item(clickable = true, a(href="$(Router.link_to(:get_mixpxy))", style="text-decoration: none;; color: #000000;","Binary <em>pxy</em> diagram")),
            ]
                )
                ))