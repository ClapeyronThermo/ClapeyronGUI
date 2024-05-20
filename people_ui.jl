include("layout.jl")
@eval page_layout( [
h1("Developers")
p("")
row([cell(class="st-col col-3", [
        imageview(src="https://media.licdn.com/dms/image/D5603AQEGMAZzdPepcw/profile-displayphoto-shrink_400_400/0/1695345101523?e=1721260800&v=beta&t=5u4nGrMktND95FVkQ5JIRbRWWzdSh3f_jp8vsAhCI4I",
        height = "350px",
        width = "300px"),
        ]),
      cell(class="st-col col-7", [
        h4("Pierre J Walker"),
        p("<font size=3> Pierre is the co-founder and principal developer of Clapeyron.jl. He obtained his MEng in Chemical Engineering from Imperial College London, with a year abroad at MIT. While at Imperial, he worked with Prof. George Jackson and Dr. Andrew Haslam on the development of new thermodynamic models and ab initio methods. At MIT, he worked with Prof. William Green on the development of machine learning and group-contribution methods for the prediction of solvation properties. He is currently a PhD candidate at California Institute of Technology, working with Prof. Zhen-Gang Wang where his research focuses on everaging computational tools and theory to model thermodynamic and transport properties of ionic liquid and polyelectrolyte systems. </font>  "),
        row([a(btn(icon="img:https://raw.githubusercontent.com/simple-icons/simple-icons/6a137bea1fd652227fe9c52d2aafdabe68d2f4af/icons/github.svg"),href="https://github.com/pw0908"),
        a(btn(icon="img:https://raw.githubusercontent.com/simple-icons/simple-icons/master/icons/googlescholar.svg"),href="https://scholar.google.com/citations?user=TU9sngYAAAAJ&hl=en"),
        a(btn(icon="img:https://raw.githubusercontent.com/simple-icons/simple-icons/6a137bea1fd652227fe9c52d2aafdabe68d2f4af/icons/linkedin.svg"),href="https://www.linkedin.com/in/pierre-walker/")])
        # h6("Awards:"),
        # p("2023: Constantin G. Economou Prize"),
        # p("2021: Hinchley Medal"),
        # p("2020: Associate Fellowship of the Higher Education Academy")
        ])
       ])

row([cell(class="st-col col-7", [
        h4("Andrés Riedemann"),
        p("Text"),
        row([a(btn(icon="img:https://raw.githubusercontent.com/simple-icons/simple-icons/6a137bea1fd652227fe9c52d2aafdabe68d2f4af/icons/github.svg"),href="https://github.com/longemen3000"),
        a(btn(icon="img:https://raw.githubusercontent.com/simple-icons/simple-icons/6a137bea1fd652227fe9c52d2aafdabe68d2f4af/icons/linkedin.svg"),href="https://www.linkedin.com/in/andrés-riedemann-rubilar-b89241a1/")])
        ]),
    cell(class="st-col col-3", [
      imageview(src="https://pbs.twimg.com/profile_images/1040248301554212869/6NMXLMrP_400x400.jpg",
      height = "auto",
      width = "300px"),
        ])
       ])
])
