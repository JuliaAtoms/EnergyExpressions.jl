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
        "Slater determinants" => "slater_determinants.md",
        "N-body operators" => "nbody_operators.md",
        "N-body matrix elements" => "nbody_matrix_elements.md",
        "Common N-body operators" => "common_operators.md",
    ],
    assets = ["assets/latex.js"],
)

deploydocs(repo = "github.com/JuliaAtoms/EnergyExpressions.jl.git")
