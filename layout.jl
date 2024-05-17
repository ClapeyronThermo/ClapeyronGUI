function page_layout( content)
StippleUI.layout(
    [
        quasar(:header, style="background-color: #4063D8;", 
toolbar(class="text-primary", style="background-color: #4063D8; ", [
          a(btn( flat=true, dense=true, icon="img:https://raw.githubusercontent.com/ClapeyronThermo/Clapeyron.jl/master/docs/src/assets/logo.svg"), href="$(Router.link_to(:get_home))"),
          toolbartitle("Clapeyron",style = "color: #ffffff; font-family:ubuntu; font-size: 1.7em; font-weight: bold;"),
          a(btn("Research",style = "color: #ffffff"),href="$(Router.link_to(:get_research))"),
          a(btn("People",style = "color: #ffffff"),href="$(Router.link_to(:get_people))"),
          btn("Examples",style = "color: #ffffff", [StippleUI.menu(
            [item(clickable = true, a(href="$(Router.link_to(:get_purept))", style="text-decoration: none;; color: #000000;","Saturation curve")),
            item(clickable = true, a(href="$(Router.link_to(:get_purerhot))", style="text-decoration: none;; color: #000000;","VLE envelope")),
            item(clickable = true, a(href="$(Router.link_to(:get_mixpxy))", style="text-decoration: none;; color: #000000;","Binary <em>pxy</em> diagram")),
            ]
                )]),
       ])
),

        cell(style="margin-top:40px", content)

    ],
    view="hHh lpR fFf",
    class="window-height",
)

end
