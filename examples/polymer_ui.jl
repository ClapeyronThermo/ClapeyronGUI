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
        select(:Select_eos_1, options = :Select_eos_list_1, label = "Equation of State"),
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
                        select(:Select_eos_2, options = :Select_eos_list_2, label = "Equation of State"),
                        select(:polymer1, options = :polymer_list, label = "Polymer 1"),
                        textfield("Molecular Weight (g/mol):", :Mw_poly1),
                        select(:polymer2, options = :polymer_list, label = "Polymer 2"),
                        textfield("Molecular Weight (g/mol):", :Mw_poly2),
                        btn("Plot",style="background-color: #4063D8; color: #ffffff", @click(:new_p_button), loading=:new_p_button)])])
                                ]),
],
)
])
