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
        btn("New", @click(:new_button), loading=:new_button)
       ]),
       cell(class="st-col col-1", [
        btn("Add", @click(:add_button), loading=:add_button)
       ])
       ])
row([cell(class="st-col col-8", [
        plot(:trace, layout=:layout)
        ])])
])
