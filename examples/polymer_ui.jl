include("layout.jl")

@eval page_layout( [
row([cell(class="st-col col-8", [
        h1("Polymer Properties"),
        ])]),
row([p("")]),
tabgroup(
        :tab_selected,
        inlinelabel = true,
        class = "text-white",
        style="background-color: #4063D8 ;activecolor: #4063C0;border-radius: 25px;",
        [
                tab(name = "gas_solubility", icon = "bubble_chart", label = "Gas Solubility"),
                tab(name = "cloud_point", icon = "cloud", label = "Cloud Points")
        ],
        ),
tabpanels(
        :tab_selected,
        animated = true,
        var"transition-prev" = "scale",
        var"transition-next" = "scale",
        [
                tabpanel(name = "gas_solubility", [
row([cell(class="st-col col-7", [
        plot(:trace_T, layout=:layout_T)
        ]), 
     cell(class="st-col col-2", [
        p(""),
        select(:polymer, options = :polymer_list, label = "Polymer"),
        textfield("Molecular Weight (g/mol):", :Mw_poly),
        select(:gas, options = :gas_list, label = "Gas"),
        textfield("Temperature (K):", :temp),
        btn("Plot",style="background-color: #4063D8; color: #ffffff", @click(:new_T_button), loading=:new_T_button)])])
                ]),
        tabpanel(name = "cloud_point", [
                row([cell(class="st-col col-7", [
                        plot(:trace_p, layout=:layout_p)
                        ]), 
                        cell(class="st-col col-2", [
                        p(""),
                        select(:polymer1, options = :polymer_list, label = "Polymer 1"),
                        textfield("Molecular Weight (g/mol):", :Mw_poly1),
                        select(:polymer2, options = :polymer_list, label = "Polymer 2"),
                        textfield("Molecular Weight (g/mol):", :Mw_poly2),
                        btn("Plot",style="background-color: #4063D8; color: #ffffff", @click(:new_p_button), loading=:new_p_button)])])
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
