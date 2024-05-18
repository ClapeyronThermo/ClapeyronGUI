include("layout.jl")

@eval page_layout( [

row([cell(class="st-col col-2", [
        select(:Select_eos, options = :Select_eos_list, label = "Equation of State")
        ]), 
        cell(class="st-col col-2", [
        textfield("Species 1:", :species1)
       ]),
       cell(class="st-col col-2", [
        textfield("Species 2:", :species2)
       ]),
        cell(class="st-col col-2", [
         textfield("Temperature (K):", :temp)
        ]),
        cell(class="st-col col-1", [
                btn("Plot",style="background-color: #4063D8; color: #ffffff", @click(:new_button), loading=:new_button)
               ])
       ])
row([cell(class="st-col col-8", [
        plot(:trace, layout=:layout)
        ])])
])
