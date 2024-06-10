include("layout.jl")

@eval page_layout( [
row([cell(class="st-col col-8", [
        h1("Saturation Properties"),
        ])])
row([cell(class="st-col col-7", [
        plot(:trace, layout=:layout)
        ]), 
     cell(class="st-col col-2", [
        p(""),
        select(:Select_eos, options = :Select_eos_list, label = "Equation of State"),
        textfield("Species:", :species),
        expansionitem(label="Options", [
                toggle("Log x-axis",:log_x,color="blue-14"),
                toggle("Log y-axis",:log_y,color="blue-14"),
                ]),
        btn("Plot",style="background-color: #4063D8; color: #ffffff", @click(:new_button), loading=:new_button)])])
row([cell(class="st-ol col-1"),
        cell(class="st-col col-2", [
        select(:x_axis, options = :Select_property, label = "X-axis")
        ]),
    cell(class="st-col col-1"),
    cell(class="st-col col-2", [
        select(:y_axis, options = :Select_property, label = "Y-axis")
                ])
])])
