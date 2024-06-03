include("layout.jl")

@eval page_layout( [
row([cell(class="st-col col-8", [
        h1("Classical Density Functional Theory"),
        ])]),
row([cell(class="st-col col-5"),
     select(:Select_eos, options = :eos_list, label = "Equation of State"),
     cell(class="st-col col-4")]),
row([p("")]),
tabgroup(
        :tab_selected,
        inlinelabel = true,
        class = "text-white",
        style="background-color: #4063D8 ;activecolor: #4063C0;border-radius: 25px;",
        [
                tab(name = "sft", icon = "water_drop", label = "Surface Tension"),
                tab(name = "ift", icon = "local_drink", label = "Interfacial Tension")
        ],
        ),

tabpanels(
        :tab_selected,
        animated = true,
        var"transition-prev" = "scale",
        var"transition-next" = "scale",
        [
                tabpanel(name = "sft", [
row([cell(class="st-col col-7", [
        plot(:trace_sft, layout=:layout_sft)
        ]), 
     cell(class="st-col col-2", [
        p(""),
        textfield("Species", :species),
        textfield("Temperature (K):", :temp),
        p("<font size=2;text-align=center><b> Surface Tension (mN/m): {{sft}} </b></font>"),
        btn("Run",style="background-color: #4063D8; color: #ffffff", @click(:converge_sft_button), loading=:converge_sft_button)])])
                ]),
tabpanel(name = "ift", [
        row([cell(class="st-col col-7", [
                plot(:trace_ift, layout=:layout_ift)
                ]), 
                cell(class="st-col col-2", [
                p(""),
                textfield("Surfactant:", :surfactant),
                p("<font size=2;text-align=center><b> Interfacial Tension (mN/m): {{ift}} </b></font>"),
                btn("Run",style="background-color: #4063D8; color: #ffffff", @click(:converge_ift_button), loading=:converge_ift_button)])])
                        ]),
        ],
        )
])
