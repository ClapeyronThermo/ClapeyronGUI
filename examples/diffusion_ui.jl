include("layout.jl")

@eval page_layout( [
row([cell(class="st-col col-8", [
        h1("Transient Diffusion"),
        ])])
row([cell(class="st-col col-7", [
        plot(:trace, layout=:layout)
        ]), 
     cell(class="st-col col-3", [
        p(""),
        h6("Limonene: {{x1}} (mol/mol)"),
        slider(LinRange(0,0.25,101), :x1,label = true),
        h6("Geraniol: {{x2}} (mol/mol)"),
        slider(LinRange(0,0.25,101), :x2,label = true),
        h6("Vanillin: {{x3}} (mol/mol)"),
        slider(LinRange(0,0.25,101), :x3,label = true),
        btn("Simulate",style="background-color: #4063D8; color: #ffffff", @click(:simulate_button), loading=:simulate_button)])])
])