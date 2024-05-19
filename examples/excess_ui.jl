include("layout.jl")

@eval page_layout( [
row([cell(class="st-col col-8", [
        h1("Excess Properties"),
        ])])
row([cell(class="st-col col-7", [
        plot(:trace, layout=:layout)
        ]), 
     cell(class="st-col col-2", [
        p(""),
        select(:Select_eos, options = :Select_eos_list, label = "Equation of State"),
        textfield("Species 1:", :species1),
        textfield("Species 2:", :species2),
        textfield("Temperature (K):", :temp),
        textfield("Pressure (bar):", :pre),
        btn("Plot",style="background-color: #4063D8; color: #ffffff", @click(:new_button), loading=:new_button)])])
row([cell(class="st-ol col-1"),
    cell(class="st-col col-2", [
        select(:y_axis, options = :Select_property, label = "Y-axis")
                ])
])])
