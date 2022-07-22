using Documenter
using EnergyExpressions
using AtomicLevels

makedocs(
    modules = [EnergyExpressions],
    sitename = "EnergyExpressions",
    pages = [
        "Home" => "index.md",
        "Theory" => [
            "Notation" => "notation.md",
            "Energy Expressions" => "energy_expressions.md",
            "Calculus of Variations" => "calculus_of_variations.md"
        ],
        "Implementation" => [
            "Conjugate orbitals" => "conjugate_orbitals.md",
            "Slater determinants" => "slater_determinants.md",
            "N-body operators" => "nbody_operators.md",
            "N-body matrix elements" => "nbody_matrix_elements.md",
            "Common N-body operators" => "common_operators.md",
            "N-body equations" => "equations.md",
            "Variation" => "variations.md",
            "System of equations" => "system_of_equations.md"
        ]
    ],
    format = Documenter.HTML(
        mathengine = MathJax2(Dict(:TeX => Dict(
            :equationNumbers => Dict(:autoNumber => "AMS"),
            :Macros => Dict(
                :defd => "≝",
                :ket => ["|#1\\rangle",1],
                :bra => ["\\langle#1|",1],
                :braket => ["\\langle#1|#2\\rangle",2],
                :matrixel => ["\\langle#1|#2|#3\\rangle",3],
                :vec => ["\\mathbf{#1}",1],
                :mat => ["\\mathsf{#1}",1],
                :conj => ["#1^*",1],
                :im => "\\mathrm{i}",
                :operator => ["\\mathfrak{#1}",1],
                :Hamiltonian => "\\operator{H}",
                :hamiltonian => "\\operator{h}",
                :Lagrangian => "\\operator{L}",
                :fock => "\\operator{f}",
                :lagrange => ["\\epsilon_{#1}",1],
                :vary => ["\\delta_{#1}",1],
                :onebody => ["(#1|#2)",2],
                :twobody => ["[#1|#2]",2],
                :twobodydx => ["[#1||#2]",2],
                :direct => ["{\\operator{J}_{#1}}",1],
                :exchange => ["{\\operator{K}_{#1}}",1],
                :diff => ["\\mathrm{d}#1\\,",1],
            ),
        ))),
    ),
    doctest = false
)

deploydocs(
    repo = "github.com/JuliaAtoms/EnergyExpressions.jl.git",
    push_preview = true,
)
