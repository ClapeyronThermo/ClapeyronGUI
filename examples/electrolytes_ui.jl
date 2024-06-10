include("layout.jl")

@eval page_layout( [
row([cell(class="st-col col-8", [
                h1("Electrolyte Properties")
                ])]),
row([cell(class="st-col col-3"),
        select(:Select_eos, options = :Select_eos_list, label = "Equation of State"),
        select(:cation, options = :Select_cation_list, label = "Cation"),
        select(:anion, options = :Select_anion_list, label = "Anion"),
        cell(class="st-col col-5")]),
row([p("")]),
tabgroup(
        :tab_selected,
        inlinelabel = true,
        class = "text-white",
        style="background-color: #4063D8 ;activecolor: #4063C0;border-radius: 25px;",
        [
                tab(name = "miac", icon = "opacity", label = "Mean Ionic Activity Coefficient"),
                tab(name = "oc", icon = "science", label = "Osmotic Coefficient"),
                tab(name = "sat", icon = "thermostat", label = "Saturation Pressure")
        ],
        ),
tabpanels(
        :tab_selected,
        animated = true,
        var"transition-prev" = "scale",
        var"transition-next" = "scale",
        [
                tabpanel(name = "miac", [
                row([cell(class="st-col col-7", [
                        plot(:trace_T, layout=:layout_T),
                        ]),
                    cell(class="st-col col-1", [
                        p(""),
                        textfield("Temperature (K):", :temp),
                        btn("Plot",style="background-color: #4063D8; color: #ffffff", @click(:new_miac_button), loading=:new_miac_button)
                        ]),
                     cell(class="st-col col-1", [
                                
                                ])])
                ]),
        tabpanel(name = "oc", [
                row([cell(class="st-col col-7", [
                        plot(:trace_p, layout=:layout_p),
                        ]),
                    cell(class="st-col col-1", [
                        p(""),
                        textfield("Temperature (K):", :temp),
                        btn("Plot",style="background-color: #4063D8; color: #ffffff", @click(:new_oc_button), loading=:new_oc_button)
                        ]),
                     cell(class="st-col col-1", [
                                
                                ])])
                                ]),
        tabpanel(name = "sat", [
                row([cell(class="st-col col-7", [
                        plot(:trace_pT, layout=:layout_pT),
                        ]),
                        cell(class="st-col col-1", [
                        p(""),
                        textfield("Temperature (K):", :temp),
                        btn("Plot",style="background-color: #4063D8; color: #ffffff", @click(:new_sat_button), loading=:new_sat_button)
                        ]),
                        cell(class="st-col col-1", [
                                
                                ])])
                                ])
],
)

])