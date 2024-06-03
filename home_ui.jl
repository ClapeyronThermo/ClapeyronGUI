include("layout.jl")
@eval page_layout( [
    row(style="margin:0px;padding:0px;margin-bottom:2.5cm",[p("")])

    row(style="margin:0px;padding:0px",[cell(class="st-col col-2"),
    cell(
        class="st-col col-6", style="display:block;margin-left:auto;margin-right:auto;align-items:center", [h1("Welcome to")]),
    cell(class="st-col col-2")])
row(style="margin:0px;padding:0px",[cell(class="st-col col-1"),
    cell(
        class="st-col col-8", style="display:block;margin-left:auto;margin-right:auto;align-items:center", [imageview(
        src = "Clapeyron_logo.svg",
        height = "auto",
        width = "800px"
    )]),
    cell(class="st-col col-4")])
    row(style="margin:0px;padding:0px",[cell(class="st-col col-2"),
    cell(
        class="st-col col-6", style="display:block;margin-left:auto;margin-right:auto;align-items:center", [p("<font size=3;text-align=center> A versatile, open-source thermodynamics toolkit for all your modelling needs. </font>")]),
    cell(class="st-col col-2")])
row(style="margin:0px;padding:0px",[cell(class="st-col col-5"),
    cell(
        class="st-col col-1", style="display:block;margin-left:auto;margin-right:auto;align-items:center", [a(btn("Code",icon="img:https://raw.githubusercontent.com/simple-icons/simple-icons/6a137bea1fd652227fe9c52d2aafdabe68d2f4af/icons/github.svg",style = "color: #000000;background-color:white;border-radius: 5px;"),href="https://github.com/ClapeyronThermo/Clapeyron.jl")]),
        cell(
            class="st-col col-1", style="display:block;margin-left:auto;margin-right:auto;align-items:center", [a(btn("Docs",icon="book",style = "color: #ffffffff;background-color:#4063D8"),href="https://clapeyronthermo.github.io/Clapeyron.jl/dev/")]),
    cell(class="st-col col-5")])
row(style="margin:0px;padding:0px;margin-bottom:3cm",[p("")])
    separator(color = "primary")
    row(style="margin:0px;padding:0px",[p("")])

row(style="margin:0px;padding:0px;margin-bottom:0.5%",[cell(class="st-col col-4"),
    cell(
        class="st-col col-2", style="display:block;margin-left:auto;margin-right:auto;align-items:center", [h2(style="text-align:center","Basics")]),
    cell(class="st-col col-4")])
row(style="margin:0px;padding:0px",[cell(class="st-col col-4",style="border: solid; border-color:#4063D8; border-width:2pt;border-radius: 25px;margin-right:0.5%;margin-left:-0.5%;margin-bottom:0.5%",[h4(style="text-align:center","Bulk Properties"),
                                                               a(imageview(class="center",style="display: block; margin-left: auto;margin-right: auto;",src="isobar.png",width="400px"),href="$(Router.link_to(:get_bulk_properties))"),
p(style="text-align:center","Obtain pure-component bulk properties, such as densities and heat capacities, at a given temperature and pressure.")]),
                                    cell(class="st-col col-4",style="border: solid; border-color:#4063D8; border-width:2pt;border-radius: 25px;margin-bottom:0.5%",[h4(style="text-align:center","Saturation Properties"),
                                    a(imageview(class="center",style="display: block; margin-left: auto;margin-right: auto;",src="vle_envelope.png",width="400px"),href="$(Router.link_to(:get_saturation_properties))"),
                                    p(style="text-align:center","Obtain pure-component saturation properties, such as saturated liquid and vapour densities, and saturation pressures.")]),
                                    cell(class="st-col col-4",style="border: solid; border-color:#4063D8; border-width:2pt;border-radius: 25px;margin-right:-0.5%;margin-left:0.5%;margin-bottom:0.5%",[h4(style="text-align:center","Excess Properties"),
                                    a(imageview(class="center",style="display: block; margin-left: auto;margin-right: auto;",src="excess.png",width="400px"),href="$(Router.link_to(:get_excess_properties))"),
                                    p(style="text-align:center","Obtain excess properties of a binary mixture at a given temperature and pressure.")])])
                                    
row(style="margin:0px;padding:0px",[cell(class="st-col col-4",style="border: solid; border-color:#4063D8; border-width:2pt;border-radius: 25px;margin-right:0.5%;margin-left:-0.5%;margin-bottom:0.5%",[h4(style="text-align:center","Binary Phase Diagrams"),
                                    a(imageview(class="center",style="display: block; margin-left: auto;margin-right: auto;",src="txy.png",width="400px"),href="$(Router.link_to(:get_binary_phase_diagrams))"),
                                    p(style="text-align:center","Obtain binary phase diagrams, such as <em>pxy</em>, <em>Txy</em> and <em>pT</em> projections.")]),
                                    cell(class="st-col col-4",style="border: solid; border-color:#4063D8; border-width:2pt;border-radius: 25px;margin-bottom:0.5%",[h4(style="text-align:center;margin-bottom:7%","Ternary Phase Diagrams"),
                                    a(imageview(class="center",style="display: block; margin-left: auto;margin-right: auto;",src="ternary.png",width="400px"),href="$(Router.link_to(:get_ternary_phase_diagrams))"),
                                    p(style="text-align:center","Obtain ternary phase diagrams, for liquid–liquid and vapour–liquid equilibria, at a given temperature and pressure.")]),
                                    cell(class="st-col col-4",style="border: solid; border-color:#4063D8; border-width:2pt;border-radius: 25px;margin-right:-0.5%;margin-left:0.5%;margin-bottom:0.5%",[h4(style="text-align:center","Parameter Estimation"),
                                    a(imageview(class="center",style="display: block; margin-left: auto;margin-right: auto;",src="fitting.png",width="400px"),href="$(Router.link_to(:get_parameter_estimation))"),
                                    p(style="text-align:center","Estimate the PC-SAFT parameters for a pure-component using pseudo-experimental data.")])])
row(style="margin:0px;padding:0px;margin-bottom:0.5%",[cell(class="st-col col-3"),
    cell(
        class="st-col col-4", style="display:block;margin-left:auto;margin-right:auto;align-items:center", [h2(style="text-align:center","Applications")]),
    cell(class="st-col col-3")])
row(style="margin:0px;padding:0px",[cell(class="st-col col-4",style="border: solid; border-color:#4063D8; border-width:2pt;border-radius: 25px;margin-right:0.5%;margin-left:-0.5%;margin-bottom:0.5%",[h4(style="text-align:center","Pharmaceuticals"),
                                                               a(imageview(class="center",style="display: block; margin-left: auto;margin-right: auto;",src="api.png",width="400px"),href="$(Router.link_to(:get_api_solubility))"),
                                                               p(style="text-align:center","Obtain API solubilities for a range of temperatures and mixed solvents.")]),
    cell(class="st-col col-4",style="border: solid; border-color:#4063D8; border-width:2pt;border-radius: 25px;margin-bottom:0.5%",[h4(style="text-align:center","Polymers"),
                               a(imageview(class="center",style="display: block; margin-left: auto;margin-right: auto;",src="polymer.png",width="400px"),href="$(Router.link_to(:get_polymer_properties))"),
                               p(style="text-align:center","Obtain polymer-gas solubilities and polymer blend cloud points using PC-SAFT.")]),
    cell(class="st-col col-4",style="border: solid; border-color:#4063D8; border-width:2pt;border-radius: 25px;margin-right:-0.5%;margin-left:0.5%;margin-bottom:0.5%",[h4(style="text-align:center","Electrolytes")])])

row(style="margin:0px;padding:0px;margin-bottom:0.5%",[cell(class="st-col col-3"),
    cell(
        class="st-col col-4", style="display:block;margin-left:auto;margin-right:auto;align-items:center", [h2(style="text-align:center","Extensions")]),
    cell(class="st-col col-3")])
row(style="margin:0px;padding:0px",[cell(class="st-col col-4",style="border: solid; border-color:#4063D8; border-width:2pt;border-radius: 25px;margin-right:0.5%;margin-left:-0.5%;margin-bottom:0.5%",[h4(style="text-align:center","Uncertainty Propagation"),
a(imageview(class="center",style="display: block; margin-left: auto;margin-right: auto;",src="uncertainty.png",width="400px"),href="$(Router.link_to(:get_uncertainty_propagation))"),
p(style="text-align:center","Observe the effect of uncertainties within the PC-SAFT parameters for a pure-component.")]),
    cell(class="st-col col-4",style="border: solid; border-color:#4063D8; border-width:2pt;border-radius: 25px;margin-bottom:0.5%",[h4(style="text-align:center","Classical DFT")]),
    cell(class="st-col col-4",style="border: solid; border-color:#4063D8; border-width:2pt;border-radius: 25px;margin-right:-0.5%;margin-left:0.5%;margin-bottom:0.5%",[h4(style="text-align:center","Transient Diffusion"),
    a(imageview(class="center",style="display: block; margin-left: auto;margin-right: auto;",src="diffusion.png",width="400px"),href="$(Router.link_to(:get_diffusion))"),
    p(style="text-align:center","Predict the gas-phase diffusion from a point source of fragrance molecules using UNIFAC.")])])
])
