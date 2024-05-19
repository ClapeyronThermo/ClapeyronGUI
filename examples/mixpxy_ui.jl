include("layout.jl")

@eval page_layout( [
row([cell(class="st-col col-8", [
                h1("Binary Phase Diagrams")
                ])]),
row([cell(class="st-col col-3"),
        select(:Select_eos, options = :Select_eos_list, label = "Equation of State"),
        textfield("Species 1:", :species1),
        textfield("Species 2:", :species2),
        cell(class="st-col col-5")]),
row([p("")]),
tabgroup(
        :tab_selected,
        inlinelabel = true,
        class = "text-white",
        style="background-color: #4063D8 ;activecolor: #4063C0;border-radius: 25px;",
        [
                tab(name = "pxy", icon = "thermostat", label = "Isothermal"),
                tab(name = "Txy", icon = "compress", label = "Isobaric"),
                tab(name = "pT", icon = "map", label = "pT Projection")
        ],
        ),
tabpanels(
        :tab_selected,
        animated = true,
        var"transition-prev" = "scale",
        var"transition-next" = "scale",
        [
                tabpanel(name = "pxy", [
                row([cell(class="st-col col-7", [
                        plot(:trace_T, layout=:layout_T),
                        ]),
                    cell(class="st-col col-1", [
                        p(""),
                        textfield("Temperature (K):", :temp),
                        btn("Plot",style="background-color: #4063D8; color: #ffffff", @click(:new_T_button), loading=:new_T_button)
                        ]),
                     cell(class="st-col col-1", [
                                
                                ])])
                ]),
        tabpanel(name = "Txy", [
                row([cell(class="st-col col-7", [
                        plot(:trace_p, layout=:layout_p),
                        ]),
                    cell(class="st-col col-1", [
                        p(""),
                        textfield("Pressure (bar):", :pre),
                        expansionitem(label="Options", [
                checkbox("Check LLE",:check_lle),
                p("{{log_x}}"),
                ]),
                        btn("Plot",style="background-color: #4063D8; color: #ffffff", @click(:new_p_button), loading=:new_p_button)
                        ]),
                     cell(class="st-col col-1", [
                                
                                ])])
                                ]),
],
)

])