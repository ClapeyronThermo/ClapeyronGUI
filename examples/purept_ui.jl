include("layout.jl")

page_layout( [
row([cell(class="st-col col-3", [
        select(:Select_eos, options = :Select_eos_list, label = "Equation of State"),
        ]), 
        cell(class="st-col col-3", [
        textfield("Species:", :species)
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
