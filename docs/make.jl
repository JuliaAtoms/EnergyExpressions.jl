using Documenter
import EnergyExpressions
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
        "Slater determinants" => "slater_determinants.md"
    ],
    assets = ["assets/latex.js"],
)

deploydocs(repo = "github.com/JuliaAtoms/EnergyExpressions.jl.git")
