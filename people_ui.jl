include("layout.jl")
page_layout( [
row([cell(class="st-col col-3", [
        imageview(src="https://media.licdn.com/dms/image/D5603AQEGMAZzdPepcw/profile-displayphoto-shrink_400_400/0/1695345101523?e=1721260800&v=beta&t=5u4nGrMktND95FVkQ5JIRbRWWzdSh3f_jp8vsAhCI4I",
        height = "350px",
        width = "300px"),
        ]),
      cell(class="st-col col-7", [
        h4("Pierre J Walker"),
        ])
       ])

row([cell(class="st-col col-7", [
        h4("Andrés Riedemann"),
        p("Text"),
        ]),
    cell(class="st-col col-3", [
      imageview(src="https://pbs.twimg.com/profile_images/1040248301554212869/6NMXLMrP_400x400.jpg",
      height = "auto",
      width = "300px"),
        ])
       ])
])
