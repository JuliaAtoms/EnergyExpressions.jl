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

]}
