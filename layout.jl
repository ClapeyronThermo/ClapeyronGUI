function page_layout(content)
StippleUI.layout(
    [
        quasar(:header, style="background-color: #4063D8;", 
toolbar(class="text-primary", style="background-color: #4063D8; ", [
          a(btn( flat=true, dense=true, icon="img:https://raw.githubusercontent.com/ClapeyronThermo/Clapeyron.jl/master/docs/src/assets/logo.svg"), href="$(Router.link_to(:get_home))"),
          toolbartitle("Clapeyron",style = "color: #ffffff; font-size: 1.7em; font-weight: bold;"),
          a(btn("Publications",style = "color: #ffffff"),href="$(Router.link_to(:get_research))"),
          a(btn("People",style = "color: #ffffff"),href="$(Router.link_to(:get_people))"),
          btn("Examples",style = "color: #ffffff", [StippleUI.menu(
            [item(clickable = true, a(href="$(Router.link_to(:get_bulk_properties))", style="text-decoration: none;; color: #000000;","Bulk properties")),
            item(clickable = true, a(href="$(Router.link_to(:get_saturation_properties))", style="text-decoration: none;; color: #000000;","Saturation properties")),
            item(clickable = true, a(href="$(Router.link_to(:get_excess_properties))", style="text-decoration: none;; color: #000000;","Excess properties")),
            item(clickable = true, a(href="$(Router.link_to(:get_binary_phase_diagrams))", style="text-decoration: none;; color: #000000;","Binary phase diagrams")),
            item(clickable = true, a(href="$(Router.link_to(:get_ternary_phase_diagrams))", style="text-decoration: none;; color: #000000;","Ternary phase diagrams")),
            item(clickable = true, a(href="$(Router.link_to(:get_parameter_estimation))", style="text-decoration: none;; color: #000000;","Parameter estimation")),
            item(clickable = true, a(href="$(Router.link_to(:get_api_solubility))", style="text-decoration: none;; color: #000000;","API solubility")),
            item(clickable = true, a(href="$(Router.link_to(:get_polymer_properties))", style="text-decoration: none;; color: #000000;","Polymer properties")),
            item(clickable = true, a(href="$(Router.link_to(:get_electrolyte_properties))", style="text-decoration: none;; color: #000000;","Electrolytes")),
            item(clickable = true, a(href="$(Router.link_to(:get_uncertainty_propagation))", style="text-decoration: none;; color: #000000;","Uncertainty propagation")),
            item(clickable = true, a(href="$(Router.link_to(:get_cdft))", style="text-decoration: none;; color: #000000;","Classical DFT")),
            item(clickable = true, a(href="$(Router.link_to(:get_diffusion))", style="text-decoration: none;; color: #000000;","Transient Diffusion"))
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

export page_layout