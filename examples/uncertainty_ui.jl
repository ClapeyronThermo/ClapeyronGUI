include("layout.jl")

@eval page_layout( [
row([cell(class="st-col col-8", [
        h1("Uncertainty Propagation"),
        ])])
row([p("")])
row([cell(class="st-col col-7", [
        tabgroup(
        :tab_selected,
        inlinelabel = true,
        class = "text-white",
        style="background-color: #4063D8 ;activecolor: #4063C0;border-radius: 25px;",
        [
                tab(name = "pT", icon = "thermostat", label = "Saturation Curve"),
                tab(name = "rhoT", icon = "opacity", label = "VLE Envelope"),
                tab(name = "bulk", icon = "local_drink", label = "Bulk Properties")
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
                tabpanel(name = "bulk", [
                        h6("Pressure: {{Slider_pre}} bar"),
                        cell(class="st-col col-3",slider(LinRange(1,200,200),:Slider_pre,label = true)),
                        plot(:trace_bulk, layout=:layout_bulk),
                       
                        select(:Selected_property, options = :Select_property, label = "Property"),
                             
                ])
],
)
        ]), 
     cell(class="st-col col-1"),  
     cell(class="st-col col-2", [
        p(""),
        textfield("Species:", :species),
        btn("Plot",style="background-color: #4063D8; color: #ffffff", @click(:new_button), loading=:new_button),
        p(""),
        h6("<b> δε </b>= {{Slider_epsilon}} K"),
        slider(LinRange(0,20,101), :Slider_epsilon,label = true),
        h6("<b> δσ </b>= {{Slider_sigma}} Å"),
        slider(LinRange(0,0.5,101), :Slider_sigma,label = true),
        h6("<b> δm </b>= {{Slider_segment}}"),
        slider(LinRange(0,0.2,101), :Slider_segment,label = true)])])])
