cell(style="display: flex; justify-content: space-between; align-items: center; background-color: #4063D8; padding: 10px 50px; color: #ffffff; top: 0; width: 100%; box-sizing: border-box;", [
    cell(style="font-size: 1.5em; font-weight: bold;",[
        a("Clapeyron.jl",href="$(Router.link_to(:get_home))", style="text-decoration: none; color: #ffffff;"),
        ]
    ),
    
])
page(model, partial=true, [@yield])