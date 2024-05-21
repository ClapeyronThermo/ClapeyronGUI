include("layout.jl")

@eval page_layout( [
row([cell(class="st-col col-8", [
        h1("Ternary Phase Diagrams"),
        ])])
row([cell(class="st-col col-7", [
        plot(:trace, layout=:layout)
        ]), 
     cell(class="st-col col-2", [
        p(""),
        select(:Select_eos, options = :Select_eos_list, label = "Equation of State"),
        textfield("Species 1:", :species1),
        textfield("Species 2:", :species2),
        textfield("Species 3:", :species3),
        textfield("Temperature (K):", :temp),
        textfield("Pressure (bar):", :pre),
        btn("Plot",style="background-color: #4063D8; color: #ffffff", @click(:new_button), loading=:new_button)])])
])