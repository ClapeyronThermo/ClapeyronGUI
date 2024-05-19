include("layout.jl")

@eval page_layout( [
row([cell(class="st-col col-8", [
        h1("Bulk Properties"),
        ])]),
row([cell(class="st-col col-5"),
     select(:Select_eos, options = :Select_eos_list, label = "Equation of State"),
     textfield("Species:", :species),
     cell(class="st-col col-5")]),
row([p("")]),
tabgroup(
        :tab_selected,
        inlinelabel = true,
        class = "text-white shadow-2",
        style="background-color: #4063D8 ;activecolor: #4063C0;border-radius: 25px;",
        [
                tab(name = "isobaric", icon = "compress", label = "Isobaric"),
                tab(name = "isothermal", icon = "thermostat", label = "Isothermal")
        ],
        ),

tabpanels(
        :tab_selected,
        animated = true,
        var"transition-prev" = "scale",
        var"transition-next" = "scale",
        [
                tabpanel(name = "isothermal", [
row([cell(class="st-col col-7", [
        plot(:trace_T, layout=:layout_T)
        ]), 
     cell(class="st-col col-2", [
        p(""),
        textfield("Temperature (K):", :temp),
        textfield("Starting Pressure (bar):", :pre_start),
        textfield("Stopping Pressure (bar):", :pre_end),
        expansionitem(label="Options", [
                toggle("Log x-axis",:log_x_T),
                toggle("Log y-axis",:log_y_T),
                checkbox("Pseudo-experimental data",:exp_data),
                p("{{log_x}}"),
                ]),
        btn("Plot",style="background-color: #4063D8; color: #ffffff", @click(:new_T_button), loading=:new_T_button)])])
row([cell(class="st-ol col-1"),
        cell(class="st-col col-2", [
        select(:y_T_axis, options = :Select_property, label = "Y-axis")
        ])])
                ]),
                tabpanel(name = "isobaric", [
                        row([cell(class="st-col col-7", [
                                plot(:trace_p, layout=:layout_p)
                                ]), 
                             cell(class="st-col col-2", [
                                p(""),
                                textfield("Pressure (bar):", :pre),
                                textfield("Starting Temperature (K):", :temp_start),
                                textfield("Stopping Temperature (K):", :temp_end),
                                expansionitem(label="Options", [
                                        toggle("Log x-axis",:log_x_p),
                                        toggle("Log y-axis",:log_y_p),
                                        checkbox("Pseudo-experimental data",:exp_data),
                                        ]),
                                btn("Plot",style="background-color: #4063D8; color: #ffffff", @click(:new_p_button), loading=:new_p_button)])])
                        row([cell(class="st-ol col-1"),
                                cell(class="st-col col-2", [
                                select(:y_p_axis, options = :Select_property, label = "Y-axis")
                                ])])
                                        ]),
        ],
        )
# row([cell(class="st-col col-7", [
#         plot(:trace, layout=:layout)
#         ]), 
#      cell(class="st-col col-2", [
#         p(""),
#         select(:Select_eos, options = :Select_eos_list, label = "Equation of State"),
#         textfield("Species:", :species),
#         expansionitem(label="Options", [
#                 toggle("Log x-axis",:log_x),
#                 toggle("Log y-axis",:log_y),
#                 checkbox("Pseudo-experimental data",:exp_data),
#                 p("{{log_x}}"),
#                 ]),
#         btn("Plot",style="background-color: #4063D8; color: #ffffff", @click(:new_button), loading=:new_button)])])
# row([cell(class="st-ol col-1"),
#         cell(class="st-col col-2", [
#         select(:x_axis, options = :Select_property, label = "X-axis")
#         ]),
#     cell(class="st-col col-1"),
#     cell(class="st-col col-2", [
#         select(:y_axis, options = :Select_property, label = "Y-axis")
                # ])])
])
