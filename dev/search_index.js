var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "#EnergyExpressions.jl-1",
    "page": "Home",
    "title": "EnergyExpressions.jl",
    "category": "section",
    "text": "EnergyExpressions.jl provides the basis for deriving equations of motion for many-body quantum systems, including, but not limited to, atoms and molecules. Please see the theory section for an overview of the mathematics and physics involved."
},

{
    "location": "notation/#",
    "page": "Notation",
    "title": "Notation",
    "category": "page",
    "text": "Throughout, Einstein notation is employed, meaning that indices appearing twice on only one side of an equation are summed over, e.g.:c = braketii iff c equiv sum_i braketiiwhereasc_ka = braketakimplies no summation."
},

{
    "location": "notation/#Spin-orbitals-1",
    "page": "Notation",
    "title": "Spin-orbitals",
    "category": "section",
    "text": "beginequation\nchi_a(tau_i) defd\npsi_a(vecr_i)sigma_a(s_i)\nendequationwhere sigma(s) is either alpha (spin up) or beta (spin down).The overlap between two spin-orbitals is given bybeginequation\nbraketij defd intdifftau conjchi_i(tau) chi_j(tau)\nendequationwhich, the case of orthogonal bases (this is a matter of choice), simplifies tobraketij = delta_ij"
},

{
    "location": "notation/#One-body-matrix-elements-1",
    "page": "Notation",
    "title": "One-body matrix elements",
    "category": "section",
    "text": "beginequation\nI(ab) equiv onebodyab defd matrixelahamiltonianb\nendequationwherebeginequation\nhamiltonian defd T + V\nendequationis the (possibly time-dependent) one-body Hamiltonian."
},

{
    "location": "notation/#Two-body-matrix-elements-1",
    "page": "Notation",
    "title": "Two-body matrix elements",
    "category": "section",
    "text": "beginequation\ntwobodydxabcd defd\ntwobodyabcd - twobodyabdc\nendequationwherebeginequation\nbeginaligned\ntwobodyabcd defd\nintdifftau_1difftau_2\nconjchi_a(tau_1)\nconjchi_b(tau_2)\nfrac1r_12\nchi_c(tau_1)\nchi_d(tau_2) \n=\ndelta(sigma_asigma_c)\ndelta(sigma_bsigma_d)\nintdiffvecr_1diffvecr_2\nconjchi_a(vecr_1)\nconjchi_b(vecr_2)\nfrac1r_12\nchi_c(vecr_1)\nchi_d(vecr_2)\nendaligned\nendequationThe special casebeginequation\nF(ab) defd twobodyabab\nendequationis called the direct interaction (gives rise to the screening potential), and the other special casebeginequation\nG(ab) defd twobodyabba\nendequationis called the exchange interaction (gives rise to the non-local potential).NBbeginequation\ntwobodydxiiii = 0\nendequationsuch thatbeginequation\nfrac12twobodydxijij equiv\nsum_ijitwobodydxijij\nendequationi.e. we sum over ij, divide by two to avoid double-counting and avoid automatically the case i=j."
},

{
    "location": "notation/#Repulsion-potential-1",
    "page": "Notation",
    "title": "Repulsion potential",
    "category": "section",
    "text": "We also have the repulsion potential formed from two orbitals ab:beginequation\ntwobodyab defd\nintdifftau_1conjchi_a(tau_1)frac1r_12chi_b(tau_1)=\ndelta(sigma_asigma_b)intdiffvecr_1conjchi_a(vecr_1)frac1r_12chi_b(vecr_1)\nendequation"
},

{
    "location": "energy_expressions/#",
    "page": "Energy Expressions",
    "title": "Energy Expressions",
    "category": "page",
    "text": ""
},

{
    "location": "energy_expressions/#Energy-expressions-1",
    "page": "Energy Expressions",
    "title": "Energy expressions",
    "category": "section",
    "text": "The average energy of the system is given bybeginequation\nE_textrmav defd matrixelPsiHamiltonianPsi\nendequationwhere Psi is the (multi-electron) wavefunction of the system and Hamiltonian is the full Hamiltonian. A common approach is to approximate the wavefunction as a linear combination of Slater determinants:beginequation\nPsi approx sum_K D_KPhi_K\nendequationwhere each Phi_K constitutes an anti-symmetrized product of one-particle spin-orbitals chi_i. Depending on whether these spin-orbitals are chosen to be orthogonal or not, with respect to each other, the energy expression takes different forms.Other choices for the expansion are possible, for instance Configuration state functions. The algebra for generating energy expressions from such objects is more involved, and is not implemented in this library.The examples below are for atomic configurations, but the library is not fundamentally limited to this case."
},

{
    "location": "energy_expressions/#Orthogonal-case,-Slater–Condon-rules-1",
    "page": "Energy Expressions",
    "title": "Orthogonal case, Slater–Condon rules",
    "category": "section",
    "text": "DocTestSetup = quote\n    using AtomicLevels\n    using EnergyExpressions\nendFirst, we find all configurations possible for two configurations of helium:julia> he = spin_configurations(c\"1s2\")\n1-element Array{Configuration{SpinOrbital{Orbital{Int64}}},1}:\n 1s²\n\njulia> he_exc = spin_configurations(c\"1s 2p\")\n12-element Array{Configuration{SpinOrbital{Orbital{Int64}}},1}:\n 1s₀α 2p₋₁α\n 1s₀α 2p₋₁β\n 1s₀α 2p₀α\n 1s₀α 2p₀β\n 1s₀α 2p₁α\n 1s₀α 2p₁β\n 1s₀β 2p₋₁α\n 1s₀β 2p₋₁β\n 1s₀β 2p₀α\n 1s₀β 2p₀β\n 1s₀β 2p₁α\n 1s₀β 2p₁β"
},

{
    "location": "energy_expressions/#One-body-energy-term-1",
    "page": "Energy Expressions",
    "title": "One-body energy term",
    "category": "section",
    "text": "The Slater–Condon rules state that the one-body energy expression between two configurations isin case the configurations are identical, onebodyPhi_APhi_B = onebodyii:\njulia> OneBodyEnergyExpression(he[1],he[1])\nI(1s₀α) + I(1s₀β)\nin case the configurations differ by one orbital, onebodyPhi_APhi_B = onebodyij:\njulia> OneBodyEnergyExpression(he[1],he_exc[1])\nI(1s₀β, 2p₋₁α)\nin case the configurations differ by more than one orbital, onebodyPhi_APhi_B = 0:\njulia> OneBodyEnergyExpression(he_exc[1],he_exc[9])\n0We can easily generate the full one-body matrix:julia> one_body_hamiltonian_matrix(SpinOrbital, he_exc)\n12×12 Array{OneBodyEnergyExpression{SpinOrbital,SpinOrbital},2}:\n I(1s₀α) + I(2p₋₁α)  I(2p₋₁α, 2p₋₁β)     I(2p₋₁α, 2p₀α)     I(2p₋₁α, 2p₀β)     I(2p₋₁α, 2p₁α)     I(2p₋₁α, 2p₁β)     I(1s₀α, 1s₀β)       0                   0                  0                  0                  0\n I(2p₋₁β, 2p₋₁α)     I(1s₀α) + I(2p₋₁β)  I(2p₋₁β, 2p₀α)     I(2p₋₁β, 2p₀β)     I(2p₋₁β, 2p₁α)     I(2p₋₁β, 2p₁β)     0                   I(1s₀α, 1s₀β)       0                  0                  0                  0\n I(2p₀α, 2p₋₁α)      I(2p₀α, 2p₋₁β)      I(1s₀α) + I(2p₀α)  I(2p₀α, 2p₀β)      I(2p₀α, 2p₁α)      I(2p₀α, 2p₁β)      0                   0                   I(1s₀α, 1s₀β)      0                  0                  0\n I(2p₀β, 2p₋₁α)      I(2p₀β, 2p₋₁β)      I(2p₀β, 2p₀α)      I(1s₀α) + I(2p₀β)  I(2p₀β, 2p₁α)      I(2p₀β, 2p₁β)      0                   0                   0                  I(1s₀α, 1s₀β)      0                  0\n I(2p₁α, 2p₋₁α)      I(2p₁α, 2p₋₁β)      I(2p₁α, 2p₀α)      I(2p₁α, 2p₀β)      I(1s₀α) + I(2p₁α)  I(2p₁α, 2p₁β)      0                   0                   0                  0                  I(1s₀α, 1s₀β)      0\n I(2p₁β, 2p₋₁α)      I(2p₁β, 2p₋₁β)      I(2p₁β, 2p₀α)      I(2p₁β, 2p₀β)      I(2p₁β, 2p₁α)      I(1s₀α) + I(2p₁β)  0                   0                   0                  0                  0                  I(1s₀α, 1s₀β)\n I(1s₀β, 1s₀α)       0                   0                  0                  0                  0                  I(1s₀β) + I(2p₋₁α)  I(2p₋₁α, 2p₋₁β)     I(2p₋₁α, 2p₀α)     I(2p₋₁α, 2p₀β)     I(2p₋₁α, 2p₁α)     I(2p₋₁α, 2p₁β)\n 0                   I(1s₀β, 1s₀α)       0                  0                  0                  0                  I(2p₋₁β, 2p₋₁α)     I(1s₀β) + I(2p₋₁β)  I(2p₋₁β, 2p₀α)     I(2p₋₁β, 2p₀β)     I(2p₋₁β, 2p₁α)     I(2p₋₁β, 2p₁β)\n 0                   0                   I(1s₀β, 1s₀α)      0                  0                  0                  I(2p₀α, 2p₋₁α)      I(2p₀α, 2p₋₁β)      I(1s₀β) + I(2p₀α)  I(2p₀α, 2p₀β)      I(2p₀α, 2p₁α)      I(2p₀α, 2p₁β)\n 0                   0                   0                  I(1s₀β, 1s₀α)      0                  0                  I(2p₀β, 2p₋₁α)      I(2p₀β, 2p₋₁β)      I(2p₀β, 2p₀α)      I(1s₀β) + I(2p₀β)  I(2p₀β, 2p₁α)      I(2p₀β, 2p₁β)\n 0                   0                   0                  0                  I(1s₀β, 1s₀α)      0                  I(2p₁α, 2p₋₁α)      I(2p₁α, 2p₋₁β)      I(2p₁α, 2p₀α)      I(2p₁α, 2p₀β)      I(1s₀β) + I(2p₁α)  I(2p₁α, 2p₁β)\n 0                   0                   0                  0                  0                  I(1s₀β, 1s₀α)      I(2p₁β, 2p₋₁α)      I(2p₁β, 2p₋₁β)      I(2p₁β, 2p₀α)      I(2p₁β, 2p₀β)      I(2p₁β, 2p₁α)      I(1s₀β) + I(2p₁β)"
},

{
    "location": "energy_expressions/#Two-body-energy-term-1",
    "page": "Energy Expressions",
    "title": "Two-body energy term",
    "category": "section",
    "text": "Similar considerations apply for the two-body energy terms between two configurations. To make it more interesting, we consider lithium which has three electrons:julia> li = spin_configurations(c\"1s2 2s\")\n2-element Array{Configuration{SpinOrbital{Orbital{Int64}}},1}:\n 1s² 2s₀α\n 1s² 2s₀β\n\njulia> li_exc = spin_configurations(c\"1s 2s 2p\")\n24-element Array{Configuration{SpinOrbital{Orbital{Int64}}},1}:\n 1s₀α 2s₀α 2p₋₁α\n 1s₀α 2s₀α 2p₋₁β\n 1s₀α 2s₀α 2p₀α\n 1s₀α 2s₀α 2p₀β\n 1s₀α 2s₀α 2p₁α\n 1s₀α 2s₀α 2p₁β\n 1s₀α 2s₀β 2p₋₁α\n 1s₀α 2s₀β 2p₋₁β\n 1s₀α 2s₀β 2p₀α\n 1s₀α 2s₀β 2p₀β\n 1s₀α 2s₀β 2p₁α\n 1s₀α 2s₀β 2p₁β\n 1s₀β 2s₀α 2p₋₁α\n 1s₀β 2s₀α 2p₋₁β\n 1s₀β 2s₀α 2p₀α\n 1s₀β 2s₀α 2p₀β\n 1s₀β 2s₀α 2p₁α\n 1s₀β 2s₀α 2p₁β\n 1s₀β 2s₀β 2p₋₁α\n 1s₀β 2s₀β 2p₋₁β\n 1s₀β 2s₀β 2p₀α\n 1s₀β 2s₀β 2p₀β\n 1s₀β 2s₀β 2p₁α\n 1s₀β 2s₀β 2p₁βin case the configurations are identical, twobodydxPhi_APhi_B = twobodydxijij:\njulia> TwoBodyEnergyExpression(li[1],li[1])\n[1s₀β 1s₀α||1s₀β 1s₀α] + [2s₀α 1s₀α||2s₀α 1s₀α] + [2s₀α 1s₀β||2s₀α 1s₀β]\nin case the configurations differ by one orbital, twobodydxPhi_APhi_B = twobodydxikjk:\njulia> TwoBodyEnergyExpression(li[1],li[2])\n[2s₀α 1s₀α||2s₀β 1s₀α] + [2s₀α 1s₀β||2s₀β 1s₀β]\nin case the configurations differ by two orbitals, twobodydxPhi_APhi_B = twobodydxijkl:\njulia> TwoBodyEnergyExpression(li[1],li_exc[end])\n[1s₀α 2s₀α||2s₀β 2p₁β]\nin case the configurations differ by more than two orbital, twobodydxPhi_APhi_B = 0:\njulia> TwoBodyEnergyExpression(li_exc[1],li_exc[end])\n0Again, we can generate the full two-body matrix:julia> two_body_hamiltonian_matrix(SpinOrbital, li_exc)\n24×24 Array{TwoBodyEnergyExpression{SpinOrbital,SpinOrbital,SpinOrbital,SpinOrbital},2}:\n [2s₀α 1s₀α||2s₀α 1s₀α] + … + [2p₋₁α 2s₀α||2p₋₁α 2s₀α]  [2p₋₁α 1s₀α||2p₋₁β 1s₀α] + [2p₋₁α 2s₀α||2p₋₁β 2s₀α]    …  0\n [2p₋₁β 1s₀α||2p₋₁α 1s₀α] + [2p₋₁β 2s₀α||2p₋₁α 2s₀α]    [2s₀α 1s₀α||2s₀α 1s₀α] + … + [2p₋₁β 2s₀α||2p₋₁β 2s₀α]     0\n [2p₀α 1s₀α||2p₋₁α 1s₀α] + [2p₀α 2s₀α||2p₋₁α 2s₀α]      [2p₀α 1s₀α||2p₋₁β 1s₀α] + [2p₀α 2s₀α||2p₋₁β 2s₀α]         0\n [2p₀β 1s₀α||2p₋₁α 1s₀α] + [2p₀β 2s₀α||2p₋₁α 2s₀α]      [2p₀β 1s₀α||2p₋₁β 1s₀α] + [2p₀β 2s₀α||2p₋₁β 2s₀α]         0\n [2p₁α 1s₀α||2p₋₁α 1s₀α] + [2p₁α 2s₀α||2p₋₁α 2s₀α]      [2p₁α 1s₀α||2p₋₁β 1s₀α] + [2p₁α 2s₀α||2p₋₁β 2s₀α]         0\n [2p₁β 1s₀α||2p₋₁α 1s₀α] + [2p₁β 2s₀α||2p₋₁α 2s₀α]      [2p₁β 1s₀α||2p₋₁β 1s₀α] + [2p₁β 2s₀α||2p₋₁β 2s₀α]      …  [1s₀α 2s₀α||1s₀β 2s₀β]\n [2s₀β 1s₀α||2s₀α 1s₀α] + [2s₀β 2p₋₁α||2s₀α 2p₋₁α]      [2s₀β 2p₋₁α||2s₀α 2p₋₁β]                                  [1s₀α 2p₋₁α||1s₀β 2p₁β]\n [2s₀β 2p₋₁β||2s₀α 2p₋₁α]                               [2s₀β 1s₀α||2s₀α 1s₀α] + [2s₀β 2p₋₁β||2s₀α 2p₋₁β]         [1s₀α 2p₋₁β||1s₀β 2p₁β]\n [2s₀β 2p₀α||2s₀α 2p₋₁α]                                [2s₀β 2p₀α||2s₀α 2p₋₁β]                                   [1s₀α 2p₀α||1s₀β 2p₁β]\n [2s₀β 2p₀β||2s₀α 2p₋₁α]                                [2s₀β 2p₀β||2s₀α 2p₋₁β]                                   [1s₀α 2p₀β||1s₀β 2p₁β]\n [2s₀β 2p₁α||2s₀α 2p₋₁α]                                [2s₀β 2p₁α||2s₀α 2p₋₁β]                                …  [1s₀α 2p₁α||1s₀β 2p₁β]\n [2s₀β 2p₁β||2s₀α 2p₋₁α]                                [2s₀β 2p₁β||2s₀α 2p₋₁β]                                   [1s₀α 2s₀β||1s₀β 2s₀β] + [1s₀α 2p₁β||1s₀β 2p₁β]\n [1s₀β 2s₀α||1s₀α 2s₀α] + [1s₀β 2p₋₁α||1s₀α 2p₋₁α]      [1s₀β 2p₋₁α||1s₀α 2p₋₁β]                                  [2s₀α 2p₋₁α||2s₀β 2p₁β]\n [1s₀β 2p₋₁β||1s₀α 2p₋₁α]                               [1s₀β 2s₀α||1s₀α 2s₀α] + [1s₀β 2p₋₁β||1s₀α 2p₋₁β]         [2s₀α 2p₋₁β||2s₀β 2p₁β]\n [1s₀β 2p₀α||1s₀α 2p₋₁α]                                [1s₀β 2p₀α||1s₀α 2p₋₁β]                                   [2s₀α 2p₀α||2s₀β 2p₁β]\n [1s₀β 2p₀β||1s₀α 2p₋₁α]                                [1s₀β 2p₀β||1s₀α 2p₋₁β]                                …  [2s₀α 2p₀β||2s₀β 2p₁β]\n [1s₀β 2p₁α||1s₀α 2p₋₁α]                                [1s₀β 2p₁α||1s₀α 2p₋₁β]                                   [2s₀α 2p₁α||2s₀β 2p₁β]\n [1s₀β 2p₁β||1s₀α 2p₋₁α]                                [1s₀β 2p₁β||1s₀α 2p₋₁β]                                   [2s₀α 1s₀β||2s₀β 1s₀β] + [2s₀α 2p₁β||2s₀β 2p₁β]\n [1s₀β 2s₀β||1s₀α 2s₀α]                                 0                                                         [2p₋₁α 1s₀β||2p₁β 1s₀β] + [2p₋₁α 2s₀β||2p₁β 2s₀β]\n 0                                                      [1s₀β 2s₀β||1s₀α 2s₀α]                                    [2p₋₁β 1s₀β||2p₁β 1s₀β] + [2p₋₁β 2s₀β||2p₁β 2s₀β]\n 0                                                      0                                                      …  [2p₀α 1s₀β||2p₁β 1s₀β] + [2p₀α 2s₀β||2p₁β 2s₀β]\n 0                                                      0                                                         [2p₀β 1s₀β||2p₁β 1s₀β] + [2p₀β 2s₀β||2p₁β 2s₀β]\n 0                                                      0                                                         [2p₁α 1s₀β||2p₁β 1s₀β] + [2p₁α 2s₀β||2p₁β 2s₀β]\n 0                                                      0\n [2s₀β 1s₀β||2s₀β 1s₀β] + … + [2p₁β 2s₀β||2p₁β 2s₀β]"
},

{
    "location": "energy_expressions/#Non-orthogonal-case,-Löwdin-rules-1",
    "page": "Energy Expressions",
    "title": "Non-orthogonal case, Löwdin rules",
    "category": "section",
    "text": "This case is more complex and not yet fully implemented in EnergyExpressions.jl"
},

{
    "location": "energy_expressions/#One-body-energy-term-2",
    "page": "Energy Expressions",
    "title": "One-body energy term",
    "category": "section",
    "text": "onebodyPhi_APhi_B = (-)^i+jonebodyijD^AB_jiwhere D^AB_ji is the determinant minor.The implementation is only correct for two-electron configurations at the moment:julia> OneBodyEnergyExpression{SpinOrbital,SpinOrbital}(he[1],he[1],orthogonal=false)\nI(1s₀α) - I(1s₀α, 1s₀β) - I(1s₀β, 1s₀α) + I(1s₀β)"
},

{
    "location": "energy_expressions/#References-1",
    "page": "Energy Expressions",
    "title": "References",
    "category": "section",
    "text": "Per-Olov Löwdin (1955). Quantum Theory of Many-Particle Systems. I. Physical Interpretations by Means of Density Matrices, Natural Spin-Orbitals, and Convergence Problems in the Method of Configurational Interaction. Physical Review, 97(6), 1474–1489. 10.1103/physrev.97.1474\nPer-Olov Löwdin (1955). Quantum Theory of Many-Particle Systems. II. Study of the Ordinary Hartree-Fock Approximation. Physical Review, 97(6), 1490–1508. 10.1103/physrev.97.1490\nPer-Olov Löwdin (1955). Quantum Theory of Many-Particle Systems. III. Extension of the Hartree-Fock Scheme to Include Degenerate Systems and Correlation Effects. Physical Review, 97(6), 1509–1520. 10.1103/physrev.97.1509DocTestSetup = nothing"
},

{
    "location": "calculus_of_variations/#",
    "page": "Calculus of Variations",
    "title": "Calculus of Variations",
    "category": "page",
    "text": ""
},

{
    "location": "calculus_of_variations/#Calculus-of-Variations-1",
    "page": "Calculus of Variations",
    "title": "Calculus of Variations",
    "category": "section",
    "text": "Calculus of Variations is a mathematical technique for deriving differential equations, whose solutions fulfil a certain criterion. For the case at hand, we are interested in solutions whose total energy matches with the one set up by the energy expression for all the electrons in a system. The following variational rules are useful:beginequation\nvarybraabraketab = ketb\nendequationwhere the notation varybraabraketab means vary braketab with respect to braa.beginequation\nbegincases\nvarybralonebodylk =\nvarybralmatrixellhamiltoniank =\nhamiltonianketk\nvaryketkonebodylk =\nbralhamiltonian\nendcases\nendequationbeginequation\nbeginaligned\nvarybraatwobodydxabcd =\nvarybraatwobodyabcd-twobodyabdc =\nvarybraamatrixelabr_12^-1cd-matrixelabr_12^-1dc \n=\nmatrixelbr_12^-1dketc-\nmatrixelbr_12^-1cketd equiv\ntwobodybdketc-\ntwobodybcketddefd\n(directbd-exchangebd)ketc\nendaligned\nendequation"
},

{
    "location": "calculus_of_variations/#Helium-1",
    "page": "Calculus of Variations",
    "title": "Helium",
    "category": "section",
    "text": "Returning to our helium example (now only considering the ground state):DocTestSetup = quote\n    using AtomicLevels\n    using EnergyExpressions\nendjulia> he = spin_configurations(c\"1s2\")[1]\n1s²\n\njulia> a,b = he.orbitals\n2-element Array{SpinOrbital{Orbital{Int64}},1}:\n 1s₀α\n 1s₀βFirst we find the one- and two-body energy expressions:julia> h = OneBodyEnergyExpression(he,he)\nI(1s₀α) + I(1s₀β)\n\njulia> HC = TwoBodyEnergyExpression(he,he)\n[1s₀β 1s₀α||1s₀β 1s₀α]We can now vary these expressions with respect to the different spin-orbitals; by convention, we vary the expressions with respect to the conjugate orbitals, to derive equations for the unconjugated ones (this is important for complex orbitals, which is the case when studying time-dependent problems):julia> h.integrals\n2-element Array{OneBodyIntegral{SpinOrbital{Orbital{Int64}},SpinOrbital{Orbital{Int64}}},1}:\n I(1s₀α)\n I(1s₀β)\n\njulia> diff.(h.integrals, Ref(conj(a)))\n2-element Array{OneBodyHamiltonian,1}:\n ĥ1s₀α\n ĥ0This result means that the variation of the first integral in the one-body expression yields the one-body Hamiltonian hamiltonian acting on the orbital 1s₀α, whereas the second integral, not containing the first orbital, varies to yield zero.If we instead vary with respect to the second orbital, unconjugated this time, we getjulia> diff.(h.integrals, Ref(b))\n2-element Array{OneBodyHamiltonian,1}:\n ĥ0\n 1s₀β†ĥSimilarly, for the two-body energy expression:julia> diff.(HC.integrals, Ref(conj(a)))\n1-element Array{DirectExchangePotentials{SpinOrbital{Orbital{Int64}},SpinOrbital{Orbital{Int64}},SpinOrbital{Orbital{Int64}}},1}:\n [1s₀β||1s₀β]1s₀αDocTestSetup = nothing"
},

{
    "location": "slater_determinants/#",
    "page": "Slater determinants",
    "title": "Slater determinants",
    "category": "page",
    "text": ""
},

{
    "location": "slater_determinants/#Slater-determinants-1",
    "page": "Slater determinants",
    "title": "Slater determinants",
    "category": "section",
    "text": "Slater determinants are wavefunctions constructed from anti-symmetrized one-particle states.CurrentModule = EnergyExpressions\nDocTestSetup = quote\n    using EnergyExpressions\n    using AtomicLevels\nend"
},

{
    "location": "slater_determinants/#EnergyExpressions.SlaterDeterminant",
    "page": "Slater determinants",
    "title": "EnergyExpressions.SlaterDeterminant",
    "category": "type",
    "text": "SlaterDeterminant(orbitals::Vector{O})\n\nConstructs a Slater determinant from a set of spin-orbitals.\n\nExamples\n\njulia> SlaterDeterminant([:a, :b])\na(1)b(2) - a(2)b(1)\n\njulia> SlaterDeterminant([:a, :b, :c])\na(1)b(2)c(3) - a(1)b(3)c(2) - a(2)b(1)c(3) + a(2)b(3)c(1) + a(3)b(1)c(2) - a(3)b(2)c(1)\n\n\n\n\n\n"
},

{
    "location": "slater_determinants/#EnergyExpressions.SlaterDeterminant-Tuple{Configuration{#s1} where #s1<:SpinOrbital}",
    "page": "Slater determinants",
    "title": "EnergyExpressions.SlaterDeterminant",
    "category": "method",
    "text": "SlaterDeterminant(configuration::Configuration{<:SpinOrbital})\n\nConstructs a Slater determinant from the spin-orbitals of the spin-configuration configuration.\n\nExamples\n\njulia> SlaterDeterminant(spin_configurations(c\"1s2\")[1])\n1s₀α(1)1s₀β(2) - 1s₀α(2)1s₀β(1)\n\n\n\n\n\n"
},

{
    "location": "slater_determinants/#Base.length-Tuple{SlaterDeterminant}",
    "page": "Slater determinants",
    "title": "Base.length",
    "category": "method",
    "text": "length(slater_determinant)\n\nReturn the number of spin-orbitals in the Slater determinant.\n\n\n\n\n\n"
},

{
    "location": "slater_determinants/#EnergyExpressions.AdjointSlaterDeterminant",
    "page": "Slater determinants",
    "title": "EnergyExpressions.AdjointSlaterDeterminant",
    "category": "type",
    "text": "AdjointSlaterDeterminant(slater_determinant)\n\nRepresentation of the Hermitian conjugate (dual vector) of a Slater determinant. Constructed via the usual adjoint operator.\n\n\n\n\n\n"
},

{
    "location": "slater_determinants/#Base.adjoint-Tuple{SlaterDeterminant}",
    "page": "Slater determinants",
    "title": "Base.adjoint",
    "category": "method",
    "text": "adjoint(slater_determinant)\n\nConstruct the adjoint of slater_determinant\n\nExamples\n\njulia> SlaterDeterminant([:a, :b])\'\n[a(1)b(2) - a(2)b(1)]†\n\n\n\n\n\n"
},

{
    "location": "slater_determinants/#Base.length-Tuple{EnergyExpressions.AdjointSlaterDeterminant}",
    "page": "Slater determinants",
    "title": "Base.length",
    "category": "method",
    "text": "length(adjoint_slater_determinant)\n\nReturn the number of spin-orbitals in the adjoint Slater determinant.\n\n\n\n\n\n"
},

{
    "location": "slater_determinants/#Construction-of-Slater-determinants-1",
    "page": "Slater determinants",
    "title": "Construction of Slater determinants",
    "category": "section",
    "text": "SlaterDeterminant\nSlaterDeterminant(::Configuration{<:SpinOrbital})\nlength(::SlaterDeterminant)\nAdjointSlaterDeterminant\nadjoint(::SlaterDeterminant)\nlength(::AdjointSlaterDeterminant)"
},

{
    "location": "slater_determinants/#Example-usage-1",
    "page": "Slater determinants",
    "title": "Example usage",
    "category": "section",
    "text": "julia> sa = SlaterDeterminant([:l, :a])\nl(1)a(2) - l(2)a(1)\n\njulia> sb = SlaterDeterminant([:k, :b])\nk(1)b(2) - k(2)b(1)\n\njulia> using AtomicLevels\n\njulia> SlaterDeterminant(spin_configurations(c\"1s2\")[1])\n1s₀α(1)1s₀β(2) - 1s₀α(2)1s₀β(1) DocTestSetup = nothing"
},

{
    "location": "nbody_operators/#",
    "page": "N-body operators",
    "title": "N-body operators",
    "category": "page",
    "text": ""
},

{
    "location": "nbody_operators/#EnergyExpressions.NBodyOperator",
    "page": "N-body operators",
    "title": "EnergyExpressions.NBodyOperator",
    "category": "type",
    "text": "NBodyOperator{N}\n\nAbstract N-body operator coupling N bodies each between two Slater determinants\n\n\n\n\n\n"
},

{
    "location": "nbody_operators/#EnergyExpressions.LinearCombinationOperator",
    "page": "N-body operators",
    "title": "EnergyExpressions.LinearCombinationOperator",
    "category": "type",
    "text": "LinearCombinationOperator(operators)\n\nRepresents a linear combination of NBodyOperators.\n\n\n\n\n\n"
},

{
    "location": "nbody_operators/#N-body-operators-1",
    "page": "N-body operators",
    "title": "N-body operators",
    "category": "section",
    "text": "CurrentModule = EnergyExpressions\nDocTestSetup = quote\n    using EnergyExpressions\n    using AtomicLevels\nendNBodyOperator\nLinearCombinationOperator DocTestSetup = nothing"
},

{
    "location": "nbody_matrix_elements/#",
    "page": "N-body matrix elements",
    "title": "N-body matrix elements",
    "category": "page",
    "text": ""
},

{
    "location": "nbody_matrix_elements/#EnergyExpressions.NBodyTermFactor",
    "page": "N-body matrix elements",
    "title": "EnergyExpressions.NBodyTermFactor",
    "category": "type",
    "text": "NBodyTermFactor\n\nAbstract type for a factor in a term in a N-body matrix element expansion\n\n\n\n\n\n"
},

{
    "location": "nbody_matrix_elements/#EnergyExpressions.OrbitalOverlap",
    "page": "N-body matrix elements",
    "title": "EnergyExpressions.OrbitalOverlap",
    "category": "type",
    "text": "OrbitalOverlap(a,b)\n\nRepresents the overlap between the orbitals a and b in a N-body matrix element expansion.\n\nExamples\n\njulia> EnergyExpressions.OrbitalOverlap(:a,:b)\n⟨a|b⟩\n\n\n\n\n\n"
},

{
    "location": "nbody_matrix_elements/#EnergyExpressions.OrbitalMatrixElement",
    "page": "N-body matrix elements",
    "title": "EnergyExpressions.OrbitalMatrixElement",
    "category": "type",
    "text": "OrbitalMatrixElement(a,o,b)\n\nRepresents the N-body matrix element between the sets of orbitals a and b.\n\nExamples\n\njulia> struct MyTwoBodyOperator <: TwoBodyOperator end\n\njulia> EnergyExpressions.OrbitalMatrixElement((:a,:b), MyTwoBodyOperator(), (:c,:d))\n⟨a b|MyTwoBodyOperator()|c d⟩\n\n\n\n\n\n"
},

{
    "location": "nbody_matrix_elements/#EnergyExpressions.NBodyTerm",
    "page": "N-body matrix elements",
    "title": "EnergyExpressions.NBodyTerm",
    "category": "type",
    "text": "NBodyTerm(factors, coeff)\n\nStructure representing one term in the expansion of a N-body matrix element.\n\n\n\n\n\n"
},

{
    "location": "nbody_matrix_elements/#EnergyExpressions.NBodyMatrixElement",
    "page": "N-body matrix elements",
    "title": "EnergyExpressions.NBodyMatrixElement",
    "category": "type",
    "text": "NBodyMatrixElement(terms)\n\nStructure representing the expansion of a N-body matrix element.\n\n\n\n\n\n"
},

{
    "location": "nbody_matrix_elements/#EnergyExpressions.isabovediagonal",
    "page": "N-body matrix elements",
    "title": "EnergyExpressions.isabovediagonal",
    "category": "function",
    "text": "isabovediagonal(i::CartesianIndex{N})\n\nReturns true if the CartesianIndex i is above the hyper-diagonal.\n\n\n\n\n\n"
},

{
    "location": "nbody_matrix_elements/#EnergyExpressions.isdiagonal",
    "page": "N-body matrix elements",
    "title": "EnergyExpressions.isdiagonal",
    "category": "function",
    "text": "isdiagonal(i::CartesianIndex)\n\nReturns true if the CartesianIndex i is diagonal in any dimensions, i.e. if not all coordinates are unique. The density matrices vanish in this case, due to the Pauli principle; this is also known as the “Fermi hole”.\n\n\n\n\n\n"
},

{
    "location": "nbody_matrix_elements/#EnergyExpressions.detaxis",
    "page": "N-body matrix elements",
    "title": "EnergyExpressions.detaxis",
    "category": "function",
    "text": "detaxis(i::CartesianIndex{N})\n\nGenerate the axis index vector for the determinant minor, whose rows or columns represented by the CartesianIndex i should be omitted.\n\n\n\n\n\n"
},

{
    "location": "nbody_matrix_elements/#EnergyExpressions.detminor",
    "page": "N-body matrix elements",
    "title": "EnergyExpressions.detminor",
    "category": "function",
    "text": "detminor(k, l, A)\n\nCalculate the determinant minor of A, where the rows k and the columns l have been stricken out.\n\n\n\n\n\n"
},

{
    "location": "nbody_matrix_elements/#EnergyExpressions.cofactor",
    "page": "N-body matrix elements",
    "title": "EnergyExpressions.cofactor",
    "category": "function",
    "text": "cofactor(k, l, A)\n\nCalculate the cofactor of A, where the rows k and the columns l have been stricken out.\n\n\n\n\n\n"
},

{
    "location": "nbody_matrix_elements/#LinearAlgebra.det",
    "page": "N-body matrix elements",
    "title": "LinearAlgebra.det",
    "category": "function",
    "text": "det(A)\n\nCalculate the determinant of the matrix A whose elements are of the NBodyTerm type, by expanding the determinant along the first column. This is an expensive operation, and should only be done with relatively sparse matrices.\n\n\n\n\n\n"
},

{
    "location": "nbody_matrix_elements/#EnergyExpressions.permutation_sign",
    "page": "N-body matrix elements",
    "title": "EnergyExpressions.permutation_sign",
    "category": "function",
    "text": "permutation_sign(p)\n\nCalculate the sign of the permutation p, 1 if iseven(p), -1 otherwise.\n\n\n\n\n\n"
},

{
    "location": "nbody_matrix_elements/#EnergyExpressions.overlap_matrix",
    "page": "N-body matrix elements",
    "title": "EnergyExpressions.overlap_matrix",
    "category": "function",
    "text": "overlap_matrix(a::SlaterDeterminant, b::SlaterDeterminant[, overlaps=[]])\n\nGenerate the single-particle orbital overlap matrix, between the orbitals in the Slater determinants a and b. All orbitals are assumed to be orthogonal, except for those which are given in overlaps.\n\nExamples\n\nFirst we define two Slater determinants that have some orbitals in common:\n\njulia> sa = SlaterDeterminant([:i, :j, :l,:k̃])\ni(1)j(2)l(3)k̃(4) - i(1)j(2)l(4)k̃(3) - i(1)j(3)l(2)k̃(4) + i(1)j(3)l(4)k̃(2) + …  + i(4)j(1)l(3)k̃(2) + i(4)j(2)l(1)k̃(3) - i(4)j(2)l(3)k̃(1) - i(4)j(3)l(1)k̃(2) + i(4)j(3)l(2)k̃(1)\n\njulia> sb = SlaterDeterminant([:i, :j, :k, :l̃])\ni(1)j(2)k(3)l̃(4) - i(1)j(2)k(4)l̃(3) - i(1)j(3)k(2)l̃(4) + i(1)j(3)k(4)l̃(2) + …  + i(4)j(1)k(3)l̃(2) + i(4)j(2)k(1)l̃(3) - i(4)j(2)k(3)l̃(1) - i(4)j(3)k(1)l̃(2) + i(4)j(3)k(2)l̃(1)\n\nThe orbital overlap matrix by default is\n\njulia> overlap_matrix(sa, sb)\n4×4 SparseArrays.SparseMatrixCSC{EnergyExpressions.NBodyTerm,Int64} with 2 stored entries:\n  [1, 1]  =  1\n  [2, 2]  =  1\n\nwhich has only two non-zero entries, since only two of the orbitals are common between the Slater determinants sa and sb.\n\nWe can then define that the orbitals k̃ and l̃ are non-orthogonal:\n\njulia> overlap_matrix(sa, sb, [OrbitalOverlap(:k̃,:l̃)])\n4×4 SparseArrays.SparseMatrixCSC{EnergyExpressions.NBodyTerm,Int64} with 3 stored entries:\n  [1, 1]  =  1\n  [2, 2]  =  1\n  [4, 4]  =  ⟨k̃|l̃⟩\n\nWe can even specify that the orbital k̃ is non-orthogonal to itself (this can be useful when the k̃ is a linear combination of orthogonal orbital):\n\njulia> overlap_matrix(sa, sa, [OrbitalOverlap(:k̃,:k̃)])\n4×4 SparseArrays.SparseMatrixCSC{EnergyExpressions.NBodyTerm,Int64} with 4 stored entries:\n  [1, 1]  =  1\n  [2, 2]  =  1\n  [3, 3]  =  1\n  [4, 4]  =  ⟨k̃|k̃⟩\n\nNotice that this overlap matrix was calculated between the Slater determinant sa and itself. \n\n\n\n\n\n"
},

{
    "location": "nbody_matrix_elements/#Base.Matrix",
    "page": "N-body matrix elements",
    "title": "Base.Matrix",
    "category": "type",
    "text": "Matrix(op::QuantumOperator, slater_determinants[, overlaps])\n\nGenerate the matrix corresponding to the quantum operator op, between the different slater_determinants. It is possible to specify non-orthogonalities between single-particle orbitals in overlaps.\n\n\n\n\n\n"
},

{
    "location": "nbody_matrix_elements/#N-body-matrix-elements-1",
    "page": "N-body matrix elements",
    "title": "N-body matrix elements",
    "category": "section",
    "text": "CurrentModule = EnergyExpressions\nDocTestSetup = quote\n    using EnergyExpressions\n    using AtomicLevels\nendThe matrix element of an N-body operator between two Slater determinants may be expanded according to the Löwdin rules (which reduce to the Slater–Condon rules if all single-particle orbitals are orthogonal):beginequation\nmatrixelPhi_AOmega_nPhi_B =\nsum_p (-)^p\nmatrixelk_1k_2k_nOmega_nl_1l_2l_n\nD^AB(k_1k_2k_nl_1l_2l_n)\nendequationwhere D^AB(k_1k_2k_nl_1l_2l_n) is the determinant minor of the orbital overlap determinant D^AB with the rows k_1k_2k_n and columns l_1l_2l_n stricken out, and p runs over all permutations.In general, a term in the expansion is thus of the formbeginequation\ncmatrixelk_1k_2k_nOmega_nl_1l_2l_nbraketabbraketcddotsbraketyz\nendequationwhere c is a scalar. This is represented by NBodyTerm type.NBodyTermFactor\nOrbitalOverlap\nOrbitalMatrixElement\nNBodyTerm\nNBodyMatrixElement\nisabovediagonal\nisdiagonal\ndetaxis\ndetminor\n# indexsum\ncofactor\ndet\npermutation_sign\noverlap_matrix\nMatrix DocTestSetup = nothing"
},

{
    "location": "common_operators/#",
    "page": "Common N-body operators",
    "title": "Common N-body operators",
    "category": "page",
    "text": ""
},

{
    "location": "common_operators/#EnergyExpressions.OneBodyHamiltonian",
    "page": "Common N-body operators",
    "title": "EnergyExpressions.OneBodyHamiltonian",
    "category": "type",
    "text": "OneBodyHamiltonian\n\nThe one-body Hamiltonian, may include external fields. Is diagonal in spin, i.e.\n\n\n\n\n\n"
},

{
    "location": "common_operators/#EnergyExpressions.FieldFreeOneBodyHamiltonian",
    "page": "Common N-body operators",
    "title": "EnergyExpressions.FieldFreeOneBodyHamiltonian",
    "category": "type",
    "text": "FieldFreeOneBodyHamiltonian\n\nThe one-body Hamiltonian, with no external fields. Is diagonal in the orbitals, i.e. does not couple unequal orbitals.\n\nExamples\n\njulia> EnergyExpressions.OrbitalMatrixElement((:a,:b), CoulombInteraction(), (:c,:d))\n[a b|c d]\n\n\n\n\n\n"
},

{
    "location": "common_operators/#EnergyExpressions.CoulombInteraction",
    "page": "Common N-body operators",
    "title": "EnergyExpressions.CoulombInteraction",
    "category": "type",
    "text": "CoulombInteraction\n\nTwo-body Hamiltonian, representing the mutual Coulombic repulsion between two electrons. Is diagonal in spin, i.e. the spin of the orbitals associated with the same coordinate must be the same.\n\nExamples\n\njulia> EnergyExpressions.OrbitalMatrixElement((:a,:b), CoulombInteraction(), (:c,:d))\n[a b|c d]\n\n\n\n\n\n"
},

{
    "location": "common_operators/#Base.iszero",
    "page": "Common N-body operators",
    "title": "Base.iszero",
    "category": "function",
    "text": "iszero(me::EnergyExpressions.OrbitalMatrixElement{1,<:SpinOrbital,OneBodyHamiltonian,<:SpinOrbital})\n\nThe matrix element vanishes if the spin-orbitals do not have the same spin.\n\n\n\n\n\niszero(me::EnergyExpressions.OrbitalMatrixElement{2,<:SpinOrbital,CoulombInteraction,<:SpinOrbital})\n\nThe matrix element vanishes if the spin-orbitals associated with the same coordinate do not have the same spin.\n\n\n\n\n\n"
},

{
    "location": "common_operators/#Common-N-body-operators-1",
    "page": "Common N-body operators",
    "title": "Common N-body operators",
    "category": "section",
    "text": "DocTestSetup = quote\n    using EnergyExpressions\n    using AtomicLevels\nendA few N-body operators that are common in quantum mechanics. It is possible to form a composite Hamiltonian thus:julia> H = OneBodyHamiltonian() + CoulombInteraction()\nĥ + ĝIf then define a set of Slater determinants, we can easily form the energy expression:julia> slaters = SlaterDeterminant.([[:a, :b], [:c, :d], [:a, :c]])\n3-element Array{SlaterDeterminant{Symbol},1}:\n a(1)b(2) - a(2)b(1)\n c(1)d(2) - c(2)d(1)\n a(1)c(2) - a(2)c(1)\n\njulia> Matrix(H, slaters)\n3×3 Array{EnergyExpressions.NBodyMatrixElement,2}:\n (a|a) + (b|b) - [a b|b a] + [a b|a b]  - [a b|d c] + [a b|c d]                (b|c) - [a b|c a] + [a b|a c]\n - [c d|b a] + [c d|a b]                (c|c) + (d|d) - [c d|d c] + [c d|c d]  - (d|a) - [c d|c a] + [c d|a c]\n (c|b) - [a c|b a] + [a c|a b]          - (a|d) - [a c|d c] + [a c|c d]        (a|a) + (c|c) - [a c|c a] + [a c|a c]We can also specify that e.g. the orbitals b and c are non-orthogonal, and thus derive a slightly different energy expression, that takes this into account:julia> Matrix(H, slaters, [OrbitalOverlap(:b,:c)])\n3×3 Array{EnergyExpressions.NBodyMatrixElement,2}:\n (a|a) + (b|b) - [a b|b a] + [a b|a b]       - ⟨b|c⟩(a|d) - [a b|d c] + [a b|c d]   ⟨b|c⟩(a|a) + (b|c) - [a b|c a] + [a b|a c]\n - ⟨c|b⟩(d|a) - [c d|b a] + [c d|a b]        (c|c) + (d|d) - [c d|d c] + [c d|c d]  - (d|a) - [c d|c a] + [c d|a c]\n ⟨c|b⟩(a|a) + (c|b) - [a c|b a] + [a c|a b]  - (a|d) - [a c|d c] + [a c|c d]        (a|a) + (c|c) - [a c|c a] + [a c|a c]Beware that the computational complexity grows factorially with the amount of non-orthogonal orbitals!OneBodyHamiltonian\nFieldFreeOneBodyHamiltonian\nCoulombInteraction\niszero DocTestSetup = nothing"
},

]}
