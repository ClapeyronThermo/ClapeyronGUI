include("layout.jl")
@eval page_layout( [
# row([cell(class="st-col col-10", [
#         h3("Topics")
#         # list()
#         ])    
#        ])
# row([cell(class="st-col col-3", [
#         select(:Select_eos, options = :Select_eos_list, label = "Equation of State"),
#         ]),
#       cell(class="st-col col-7", [
#         h4("Predictive Thermodynamics"),
#         ])
#        ])

# row([cell(class="st-col col-7", [
#         h4("Computer-Aided Molecular Design"),
#         p("Text"),
#         ]),
#     cell(class="st-col col-3", [
#         select(:Select_eos, options = :Select_eos_list, label = "Equation of State"),
#         ])
#        ])

# row([cell(class="st-col col-3", [
#         select(:Select_eos, options = :Select_eos_list, label = "Equation of State"),
#         ]),
#         cell(class="st-col col-7", [
#             h4("Seamless Multi-scale Modelling"),
#             ])    
#        ])

row([cell(class="st-col col-10", [
        h3("Publication List"),
        h5("Software"),
        p("[1] Walker, P.J., Yew, H-W., Riedeman, A., 2022. ‘Clapeyron.jl: An Extensible, Open-Source Fluid Thermodynamics Toolkit’, Industrial and Engineering Chemistry Research, 61, 20, pp. 7130-7153"),
        p("[2] Walker, P.J., Riedemann, A., Wang, Z-G., 2024. ‘GCIdentifier.jl: A Julia package for identifying molecular fragments from SMILES’, Journal of Open-Source Software, 9, 96, pp. 6453."),
        h5("Applications"),
        p("[3] Walker, P.J., Mueller, S., Smirnova, I., 2023. ‘Confidence Interval and Uncertainty Propagation Analysis of SAFT-type Equations of State’, Journal of Chemical and Engineering Data, 69, 2, pp. 495-508."),
        p("[4] Walker, P.J., 2022. ‘Towards advanced, predictive mixing rules in SAFT equations of state’, Industrial Engineering and Chemistry Research, 61, 49, pp. 18165–18175."),
        p("[5] Ylitalo, A.S., Chao, H., Walker., P.J., Fitzgibbons, T.C., Ginzburg, V.G., Zhou, W., Wang, Z.-G., Di Maio, E., Kornfield, J.A., 2022. ‘Competition between CO2-philicity and mixing entropy leads to CO2 solubility maximum in polyether polyols’, Industrial and Engineering Chemistry Research, 61, 34, pp. 12835–12844."),
        p("[6] Walker, P.J., Zhao, T., Haslam, A.J., Jackson, G., 2022. ‘Ab initio development of generalized Lennard-Jones (Mie) force fields for predictions of thermodynamic properties in advanced molecular-based SAFT equations of state’, Journal of Chemical Physics, 156, pp. 154106."),
        p("[7] Walker, P.J., Liang, X., Kontogeorgis, G.M., 2022. ‘Importance of the Relative Static Permittivity in Electrolyte SAFT-VR Mie Equations of state’, Fluid Phase Equilibria, 551, pp. 113256"),
        p("[8] Inguva, P.K., Walker, P.J., Zhu, K., Yew, H-W., Haslam, A.J., Matar, O.K., 2021. ‘Continuum-scale modelling of polymer blends using the Cahn–Hilliard equation: Transport and thermodynamics’, Soft Matter, 17, pp. 5645-5665."),
        p("[9] Walker, P. J., Haslam, A. J., 2020. ‘A new predictive group-contribution ideal-heat-capacity model, and its influence on second-derivative properties calculated using a free-energy equation of state’, Journal of Chemical and Engineering Data, 65(12), pp. 5809-5829."),
        h5("Education"),
        p("[10] Paoli, L., Inguva, P.K., Haslam, A.J., Walker, P.J., 2024. ‘Introducing students to research codes: An introductory course to computational thermodynamics in Julia’, Education for Chemical Engineers, 48, pp. 1-14.")
        # list()
        ])    
       ])
])
