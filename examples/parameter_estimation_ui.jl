include("layout.jl")

@eval page_layout( [
row([cell(class="st-col col-8", [
        h1("Parameter Estimation"),
        ])])
row([p("")])
row([cell(class="st-col col-7", [
        tabgroup(
        :tab_selected,
        inlinelabel = true,
        class = "text-white",
        style="background-color: #4063D8 ;activecolor: #4063C0;border-radius: 25px;",
        [
                tab(name = "pT", icon = "thermostat", label = "Saturation Pressure"),
                tab(name = "rhoT", icon = "opacity", label = "Saturated Liquid Density"),
                tab(name = "obj", icon = "flag", label = "Objective function")
        ],
        )
tabpanels(
        :tab_selected,
        animated = true,
        var"transition-prev" = "scale",
        var"transition-next" = "scale",
        [
                tabpanel(name = "pT", [
                        plot(:trace, layout=:layout)
                ]),
                tabpanel(name = "rhoT", [
                        plot(:trace_rho, layout=:layout_rho)
                ]),
                tabpanel(name = "obj", [
                        
                        plot(:trace_obj, layout=:layout_obj),
                                                    
                ])
],
)
        ]), 
     cell(class="st-col col-1"),  
     cell(class="st-col col-2", [
        p(""),
        textfield("Species:", :species),
        btn("Generate Data",style="background-color: #4063D8; color: #ffffff", @click(:data_button), loading=:data_button),
        p(""),
        btn("Fit",style="background-color: #4063D8; color: #ffffff", @click(:fit_button), loading=:fit_button),
        p(""),
        h6("<b> ε </b>= {{epsilon}} K"),
        h6("<b> σ </b>= {{sigma}} Å"),
        h6("<b> m </b>= {{segment}}")])])])
