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
    "text": "beginequation\ntwobodydxabcd defd\ntwobodyabcd - twobodyabdc\nendequationwherebeginequation\nbeginaligned\ntwobodyabcd defd\nintdifftau_1difftau_2\nconjchi_a(tau_1)\nconjchi_b(tau_2)\nfrac1r_12\nchi_c(tau_1)\nchi_d(tau_2) \n=\ndelta(sigma_asigma_c)\ndelta(sigma_bsigma_d)\nintdiffvecr_1diffvecr_2\nconjchi_a(vecr_1)\nconjchi_b(vecr_2)\nfrac1r_12\nchi_c(vecr_1)\nchi_d(vecr_2)\nendaligned\nendequationThe special casebeginequation\nF(ab) defd twobodyabab\nendequationis called the direct interaction (gives rise to the screening potential), and the other special casebeginequation\nG(ab) defd twobodyabba\nendequationis called the exchange interaction (gives rise to the non-local potential).note: Note\nSincebeginequation\ntwobodydxiiii = 0\nendequationwe havebeginequation\nfrac12twobodydxijij equiv\nsum_ijitwobodydxijij\nendequationi.e. we sum over ij, divide by two to avoid double-counting and avoid automatically the case i=j."
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
    "text": "The average energy of the system is given bybeginequation\nE_textrmav defd matrixelPsiHamiltonianPsi\nendequationwhere Psi is the (multi-electron) wavefunction of the system and Hamiltonian is the full Hamiltonian. A common approach is to approximate the wavefunction as a linear combination of Slater determinants:beginequation\nPsi approx sum_K D_KPhi_K\nendequationwhere each Phi_K constitutes an anti-symmetrized product of one-particle spin-orbitals chi_i. Depending on whether these spin-orbitals are chosen to be orthogonal or not, with respect to each other, the energy expression takes different forms.Other choices for the expansion are possible, for instance Configuration state functions. The algebra for generating energy expressions from such objects is more involved, and is not implemented in this library.The examples below are for atomic configurations, but the library is not limited to this case, rather, it could work with any kind of configuration constructed from single-particle spin-orbitals. For example, it possible to work with normal Julia Symbols denoting the spin-orbitals. This can be used to derive equations of motion purely symbolically."
},

{
    "location": "energy_expressions/#Orthogonal-case,-Slater–Condon-rules-1",
    "page": "Energy Expressions",
    "title": "Orthogonal case, Slater–Condon rules",
    "category": "section",
    "text": "DocTestSetup = quote\n    using AtomicLevels\n    using EnergyExpressions\nendFirst, we find all configurations possible for two configurations of helium, and convert them into SlaterDeterminants.julia> he = SlaterDeterminant.(spin_configurations(c\"1s2\"))\n1-element Array{SlaterDeterminant{SpinOrbital{Orbital{Int64}}},1}:\n |1s₀α 1s₀β|\n\njulia> he_exc = SlaterDeterminant.(spin_configurations(c\"1s 2p\"))\n12-element Array{SlaterDeterminant{SpinOrbital{Orbital{Int64}}},1}:\n |1s₀α 2p₋₁α|\n |1s₀α 2p₋₁β|\n |1s₀α 2p₀α|\n |1s₀α 2p₀β|\n |1s₀α 2p₁α|\n |1s₀α 2p₁β|\n |1s₀β 2p₋₁α|\n |1s₀β 2p₋₁β|\n |1s₀β 2p₀α|\n |1s₀β 2p₀β|\n |1s₀β 2p₁α|\n |1s₀β 2p₁β|"
},

{
    "location": "energy_expressions/#One-body-energy-term-1",
    "page": "Energy Expressions",
    "title": "One-body energy term",
    "category": "section",
    "text": "First we create a OneBodyHamiltonian, the one-body operator that corresponds to the physical quantity energy:julia> h = OneBodyHamiltonian()\nĥWe can then evaluate matrix elements between Slater determinants corresponding to different configurations.The Slater–Condon rules state that the one-body energy expression between two configurations isin case the configurations are identical, onebodyPhi_APhi_B = onebodyii:\njulia> Matrix(h, he)\n1×1 Array{EnergyExpressions.NBodyMatrixElement,2}:\n (1s₀α|1s₀α) + (1s₀β|1s₀β)\nin case the configurations differ by one orbital, onebodyPhi_APhi_B = onebodyij:\njulia> Matrix(h, [he[1],he_exc[2]])\n2×2 Array{EnergyExpressions.NBodyMatrixElement,2}:\n (1s₀α|1s₀α) + (1s₀β|1s₀β)  (1s₀β|2p₋₁β)\n (2p₋₁β|1s₀β)               (1s₀α|1s₀α) + (2p₋₁β|2p₋₁β)\nWe had to choose the second excited state, since the OneBodyHamiltonian does not couple orbitals of differing spin; therefore the corresponding matrix elements would be zero. The diagonal matrix elements correspond to the first case.\nin case the configurations differ by more than one orbital, onebodyPhi_APhi_B = 0:\njulia> Matrix(h, [he[1],he_exc[10]])\n2×2 Array{EnergyExpressions.NBodyMatrixElement,2}:\n (1s₀α|1s₀α) + (1s₀β|1s₀β)  0\n 0                          (1s₀β|1s₀β) + (2p₀β|2p₀β)We can easily generate the full one-body matrix:julia> Matrix(h, vcat(he,he_exc))\n13×13 Array{EnergyExpressions.NBodyMatrixElement,2}:\n (1s₀α|1s₀α) + (1s₀β|1s₀β)  0                            (1s₀β|2p₋₁β)                 …  - (1s₀α|2p₁α)              0\n 0                          (1s₀α|1s₀α) + (2p₋₁α|2p₋₁α)  0                               0                          0\n (2p₋₁β|1s₀β)               0                            (1s₀α|1s₀α) + (2p₋₁β|2p₋₁β)     0                          0\n 0                          (2p₀α|2p₋₁α)                 0                               0                          0\n (2p₀β|1s₀β)                0                            (2p₀β|2p₋₁β)                    0                          0\n 0                          (2p₁α|2p₋₁α)                 0                            …  0                          0\n (2p₁β|1s₀β)                0                            (2p₁β|2p₋₁β)                    0                          0\n - (2p₋₁α|1s₀α)             0                            0                               (2p₋₁α|2p₁α)               0\n 0                          0                            0                               0                          (2p₋₁β|2p₁β)\n - (2p₀α|1s₀α)              0                            0                               (2p₀α|2p₁α)                0\n 0                          0                            0                            …  0                          (2p₀β|2p₁β)\n - (2p₁α|1s₀α)              0                            0                               (1s₀β|1s₀β) + (2p₁α|2p₁α)  0\n 0                          0                            0                               0                          (1s₀β|1s₀β) + (2p₁β|2p₁β)"
},

{
    "location": "energy_expressions/#Two-body-energy-term-1",
    "page": "Energy Expressions",
    "title": "Two-body energy term",
    "category": "section",
    "text": "Similar considerations apply for the two-body energy terms between two configurations. To make it more interesting, we consider lithium which has three electrons:julia> li = SlaterDeterminant.(spin_configurations(c\"1s2 2s\"))\n2-element Array{SlaterDeterminant{SpinOrbital{Orbital{Int64}}},1}:\n |1s₀α 1s₀β 2s₀α|\n |1s₀α 1s₀β 2s₀β|\n\njulia> li_exc = SlaterDeterminant.(spin_configurations(c\"1s 2s 2p\"))\n24-element Array{SlaterDeterminant{SpinOrbital{Orbital{Int64}}},1}:\n |1s₀α 2s₀α 2p₋₁α|\n |1s₀α 2s₀α 2p₋₁β|\n |1s₀α 2s₀α 2p₀α|\n |1s₀α 2s₀α 2p₀β|\n |1s₀α 2s₀α 2p₁α|\n |1s₀α 2s₀α 2p₁β|\n |1s₀α 2s₀β 2p₋₁α|\n |1s₀α 2s₀β 2p₋₁β|\n |1s₀α 2s₀β 2p₀α|\n |1s₀α 2s₀β 2p₀β|\n |1s₀α 2s₀β 2p₁α|\n |1s₀α 2s₀β 2p₁β|\n |1s₀β 2s₀α 2p₋₁α|\n |1s₀β 2s₀α 2p₋₁β|\n |1s₀β 2s₀α 2p₀α|\n |1s₀β 2s₀α 2p₀β|\n |1s₀β 2s₀α 2p₁α|\n |1s₀β 2s₀α 2p₁β|\n |1s₀β 2s₀β 2p₋₁α|\n |1s₀β 2s₀β 2p₋₁β|\n |1s₀β 2s₀β 2p₀α|\n |1s₀β 2s₀β 2p₀β|\n |1s₀β 2s₀β 2p₁α|\n |1s₀β 2s₀β 2p₁β|The operator we choose is the Coulomb repulsion, implemented by CoulombInteraction:julia> H = CoulombInteraction()\nĝin case the configurations are identical, twobodydxPhi_APhi_B = twobodydxijij:\njulia> Matrix(H, li)[1,1]\nF(1s₀α,1s₀β) - G(1s₀α,2s₀α) + F(1s₀α,2s₀α) + F(1s₀β,2s₀α)\nNB that some terms in the sum vanish due to spin-conservation in the two-body integral.\nin case the configurations differ by one orbital, twobodydxPhi_APhi_B = twobodydxikjk:\njulia> Matrix(H, [li[1],li_exc[6]])[2,1]\n- [1s₀α 2p₁β|1s₀α 1s₀β] - [2s₀α 2p₁β|2s₀α 1s₀β]\nin case the configurations differ by two orbitals, twobodydxPhi_APhi_B = twobodydxijkl:\njulia> Matrix(H, [li[1],li_exc[11]])[1,2]\n[1s₀β 2s₀α|2s₀β 2p₁α]\nin case the configurations differ by more than two orbital, twobodydxPhi_APhi_B = 0:\njulia> Matrix(H, [li[1],li_exc[23]])[1,2]\n0\nIn this particular case, the matrix element vanishes because of spin-conservation, as well.Again, we can generate the full two-body matrix:julia> Matrix(H, vcat(li,li_exc))\n26×26 Array{EnergyExpressions.NBodyMatrixElement,2}:\n F(1s₀α,1s₀β) - G(1s₀α,2s₀α) + F(1s₀α,2s₀α) + F(1s₀β,2s₀α)                 0                                                                           …  0\n 0                                                                         F(1s₀α,1s₀β) + F(1s₀α,2s₀β) - G(1s₀β,2s₀β) + F(1s₀β,2s₀β)                      0\n 0                                                                         0                                                                              0\n - [1s₀α 2p₋₁β|1s₀α 1s₀β] - [2s₀α 2p₋₁β|2s₀α 1s₀β]                         0                                                                              0\n 0                                                                         0                                                                              0\n - [1s₀α 2p₀β|1s₀α 1s₀β] - [2s₀α 2p₀β|2s₀α 1s₀β]                           0                                                                           …  0\n 0                                                                         0                                                                              0\n - [1s₀α 2p₁β|1s₀α 1s₀β] - [2s₀α 2p₁β|2s₀α 1s₀β]                           0                                                                              0\n [2s₀β 2p₋₁α|1s₀β 2s₀α]                                                    0                                                                              0\n 0                                                                         - [1s₀α 2p₋₁β|1s₀α 1s₀β] - [2s₀β 2p₋₁β|2s₀β 1s₀β] + [2s₀β 2p₋₁β|1s₀β 2s₀β]     0\n [2s₀β 2p₀α|1s₀β 2s₀α]                                                     0                                                                           …  0\n 0                                                                         - [1s₀α 2p₀β|1s₀α 1s₀β] - [2s₀β 2p₀β|2s₀β 1s₀β] + [2s₀β 2p₀β|1s₀β 2s₀β]        0\n [2s₀β 2p₁α|1s₀β 2s₀α]                                                     0                                                                              0\n 0                                                                         - [1s₀α 2p₁β|1s₀α 1s₀β] - [2s₀β 2p₁β|2s₀β 1s₀β] + [2s₀β 2p₁β|1s₀β 2s₀β]        0\n [1s₀β 2p₋₁α|1s₀β 1s₀α] + [2s₀α 2p₋₁α|2s₀α 1s₀α] - [2s₀α 2p₋₁α|1s₀α 2s₀α]  0                                                                              0\n 0                                                                         - [2s₀α 2p₋₁β|1s₀α 2s₀β]                                                    …  0\n [1s₀β 2p₀α|1s₀β 1s₀α] + [2s₀α 2p₀α|2s₀α 1s₀α] - [2s₀α 2p₀α|1s₀α 2s₀α]     0                                                                              0\n 0                                                                         - [2s₀α 2p₀β|1s₀α 2s₀β]                                                        0\n [1s₀β 2p₁α|1s₀β 1s₀α] + [2s₀α 2p₁α|2s₀α 1s₀α] - [2s₀α 2p₁α|1s₀α 2s₀α]     0                                                                              0\n 0                                                                         - [2s₀α 2p₁β|1s₀α 2s₀β]                                                        0\n 0                                                                         [1s₀β 2p₋₁α|1s₀β 1s₀α] + [2s₀β 2p₋₁α|2s₀β 1s₀α]                             …  0\n 0                                                                         0                                                                              - [1s₀β 2p₋₁β|2p₁β 1s₀β] + [1s₀β 2p₋₁β|1s₀β 2p₁β] - [2s₀β 2p₋₁β|2p₁β 2s₀β] + [2s₀β 2p₋₁β|2s₀β 2p₁β]\n 0                                                                         [1s₀β 2p₀α|1s₀β 1s₀α] + [2s₀β 2p₀α|2s₀β 1s₀α]                                  0\n 0                                                                         0                                                                              - [1s₀β 2p₀β|2p₁β 1s₀β] + [1s₀β 2p₀β|1s₀β 2p₁β] - [2s₀β 2p₀β|2p₁β 2s₀β] + [2s₀β 2p₀β|2s₀β 2p₁β]\n 0                                                                         [1s₀β 2p₁α|1s₀β 1s₀α] + [2s₀β 2p₁α|2s₀β 1s₀α]                                  0\n 0                                                                         0                                                                           …  - G(1s₀β,2s₀β) + F(1s₀β,2s₀β) - G(1s₀β,2p₁β) + F(1s₀β,2p₁β) - G(2s₀β,2p₁β) + F(2s₀β,2p₁β)"
},

{
    "location": "energy_expressions/#Linear-combination-of-N-body-operators-1",
    "page": "Energy Expressions",
    "title": "Linear combination of N-body operators",
    "category": "section",
    "text": "Finally, we may form a linear combination of different many-body operators:julia> H = OneBodyHamiltonian() + CoulombInteraction()which allows to generate the full energy-expression matrix simultaneously:julia> Matrix(H, vcat(he,he_exc[1:2]))\n3×3 Array{EnergyExpressions.NBodyMatrixElement,2}:\n (1s₀α|1s₀α) + (1s₀β|1s₀β) + F(1s₀α,1s₀β)  0                                                            (1s₀β|2p₋₁β) + [1s₀α 1s₀β|1s₀α 2p₋₁β]\n 0                                         (1s₀α|1s₀α) + (2p₋₁α|2p₋₁α) - G(1s₀α,2p₋₁α) + F(1s₀α,2p₋₁α)  0\n (2p₋₁β|1s₀β) + [1s₀α 2p₋₁β|1s₀α 1s₀β]     0                                                            (1s₀α|1s₀α) + (2p₋₁β|2p₋₁β) + F(1s₀α,2p₋₁β)"
},

{
    "location": "energy_expressions/#Non-orthogonal-case,-Löwdin-rules-1",
    "page": "Energy Expressions",
    "title": "Non-orthogonal case, Löwdin rules",
    "category": "section",
    "text": "This case is more complex and is invoked by providing a list of OrbitalOverlaps designating those pairs of orbitals which are chosen to be non-orthogonal. Beware that the computation complexity increases factorially with the amount of non-orthogonalities! It is thus not a good idea to choose all orbitals to non-orthogonal to one another.For simplicity, we now consider a set of symbolic orbitals: a,b,c, where specify that b and c are non-orthogonal:julia> cfgs = SlaterDeterminant.([[:a, :b, :c], [:b, :d, :e]])\n2-element Array{SlaterDeterminant{Symbol},1}:\n |a b c|\n |b d e|The pairwise non-orthogonality can be specified simply asjulia> overlaps = [OrbitalOverlap(:b, :c)]\n1-element Array{OrbitalOverlap{Symbol,Symbol},1}:\n ⟨b|c⟩"
},

{
    "location": "energy_expressions/#One-body-energy-term-2",
    "page": "Energy Expressions",
    "title": "One-body energy term",
    "category": "section",
    "text": "onebodyPhi_APhi_B = (-)^i+jonebodyijD^AB(ij)where D^AB(ij) is the determinant minor, where the row i and column j are stricken out.julia> Matrix(OneBodyHamiltonian(), cfgs, overlaps)\n2×2 Array{EnergyExpressions.NBodyMatrixElement,2}:\n (a|a) - (a|a)⟨c|b⟩⟨b|c⟩ + (b|b) - (b|c)⟨c|b⟩ - (c|b)⟨b|c⟩ + (c|c)  0\n 0                                                                  (b|b) + (d|d) + (e|e)"
},

{
    "location": "energy_expressions/#Two-body-energy-term-2",
    "page": "Energy Expressions",
    "title": "Two-body energy term",
    "category": "section",
    "text": "Similarly, we havematrixelPhi_AOmega_2Phi_B = (-)^i+j+k+lmatrixelijOmega_2klD^AB(ijkl)where D^AB(ijkl) is the determinant minor, where the rows i,j and columns k,l are stricken out.The expressions become rather lengthy, so we only look at one particular matrix element:julia> Matrix(CoulombInteraction(), cfgs, overlaps)[1,2]\n- ⟨c|b⟩[a b|e d] + ⟨c|b⟩[a b|d e] + [a c|e d] - [a c|d e]"
},

{
    "location": "energy_expressions/#Higher-order-terms-1",
    "page": "Energy Expressions",
    "title": "Higher-order terms",
    "category": "section",
    "text": "For details of the implementation and the general N-body case, see N-body matrix elements."
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
    "text": "Returning to our helium example (now only considering the ground state):DocTestSetup = quote\n    using AtomicLevels\n    using EnergyExpressions\nendjulia> he = SlaterDeterminant.(spin_configurations(c\"1s2\"))\n1-element Array{SlaterDeterminant{SpinOrbital{Orbital{Int64}}},1}:\n |1s₀α 1s₀β|\n\njulia> a,b = he[1].orbitals\n2-element Array{SpinOrbital{Orbital{Int64}},1}:\n 1s₀α\n 1s₀βFirst we find the energy expression:julia> E = Matrix(OneBodyHamiltonian() + CoulombInteraction(), he)\n1×1 Array{EnergyExpressions.NBodyMatrixElement,2}:\n (1s₀α|1s₀α) + (1s₀β|1s₀β) + F(1s₀α,1s₀β)We can now vary this expression with respect to the different spin-orbitals; by convention, we vary the expression with respect to the conjugate orbitals, to derive equations for the unconjugated ones (this is important for complex orbitals, which is the case when studying time-dependent problems):julia> diff(E, Conjugate(a))\n1×1 SparseArrays.SparseMatrixCSC{LinearCombinationEquation,Int64} with 1 stored entry:\n  [1, 1]  =  ĥ|1s₀α⟩ + [1s₀β|1s₀β]|1s₀α⟩\n\njulia> diff(E, Conjugate(b))\n1×1 SparseArrays.SparseMatrixCSC{LinearCombinationEquation,Int64} with 1 stored entry:\n  [1, 1]  =  ĥ|1s₀β⟩ + [1s₀α|1s₀α]|1s₀β⟩This result means that the variation of the first integral in the one-body part of the energy expression yields the one-body Hamiltonian hamiltonian acting on the orbital 1s₀α, whereas the second integral, not containing the first orbital, varies to yield zero. The two-body energy F(1s₀α,1s₀β) varies to yield the two-body potential [1s₀β|1s₀β], again acting on the orbital 1s₀α.If we instead vary with respect to the unconjugated orbitals, we getjulia> diff(E, a)\n1×1 SparseArrays.SparseMatrixCSC{LinearCombinationEquation,Int64} with 1 stored entry:\n  [1, 1]  =  ⟨1s₀α|ĥ + ⟨1s₀α|[1s₀β|1s₀β]\n\njulia> diff(E, b)\n1×1 SparseArrays.SparseMatrixCSC{LinearCombinationEquation,Int64} with 1 stored entry:\n  [1, 1]  =  ⟨1s₀β|ĥ + ⟨1s₀β|[1s₀α|1s₀α]that is, we instead derived equations for the conjugated orbitals.DocTestSetup = nothing"
},

{
    "location": "conjugate_orbitals/#",
    "page": "Conjugate orbitals",
    "title": "Conjugate orbitals",
    "category": "page",
    "text": ""
},

{
    "location": "conjugate_orbitals/#EnergyExpressions.Conjugate",
    "page": "Conjugate orbitals",
    "title": "EnergyExpressions.Conjugate",
    "category": "type",
    "text": "Conjugate(orbital)\n\nType representing the conjugation of an orbital.\n\nExamples\n\njulia> Conjugate(:a)\n:a†\n\n\n\n\n\n"
},

{
    "location": "conjugate_orbitals/#Base.conj",
    "page": "Conjugate orbitals",
    "title": "Base.conj",
    "category": "function",
    "text": "conj(o::AbstractOrbital)\n\nConvenience function to conjugate an AbstractOrbital.\n\nExamples\n\njulia> conj(o\"1s\")\n1s†\n\n\n\n\n\nconj(o::Conjugate)\n\nConvenience function to unconjugate a conjugated orbital.\n\nExamples\n\njulia> conj(Conjugate(:a))\n:a\n\n\n\n\n\n"
},

{
    "location": "conjugate_orbitals/#Conjugate-orbitals-1",
    "page": "Conjugate orbitals",
    "title": "Conjugate orbitals",
    "category": "section",
    "text": "DocTestSetup = quote\n    using EnergyExpressions\n    using AtomicLevels\nendA conjugated orbital, conjchi (often written brachi) is the dual to the unconjugated orbital chi (often written ketchi). In the code, conjugation of orbitals is denoted with a dagger (†) to avoid confusion with the multiplication operator *.Conjugate\nconj DocTestSetup = nothing"
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
    "location": "nbody_operators/#EnergyExpressions.IdentityOperator",
    "page": "N-body operators",
    "title": "EnergyExpressions.IdentityOperator",
    "category": "type",
    "text": "IdentityOperator{N}\n\nThe N-body identity operator. Leaves the orbital(s) acted upon unchanged.\n\n\n\n\n\n"
},

{
    "location": "nbody_operators/#EnergyExpressions.ContractedOperator",
    "page": "N-body operators",
    "title": "EnergyExpressions.ContractedOperator",
    "category": "type",
    "text": "ContractedOperator(a, o, b)\n\nAn NBodyOperator representing the contraction of the operator o over the orbital sets a and b. The lengths of a and b have to equal, and they cannot exceed the dimension of o.\n\n\n\n\n\n"
},

{
    "location": "nbody_operators/#Base.in-Union{Tuple{O}, Tuple{O,ContractedOperator}} where O",
    "page": "N-body operators",
    "title": "Base.in",
    "category": "method",
    "text": "in(orbital, co::ContractedOperator)\n\nTest if orbital is among the right set of orbitals of the ContractedOperator co. Useful to test if co is an integral operator with respect to orbital.\n\n\n\n\n\n"
},

{
    "location": "nbody_operators/#Base.in-Union{Tuple{O}, Tuple{Conjugate{O},ContractedOperator}} where O",
    "page": "N-body operators",
    "title": "Base.in",
    "category": "method",
    "text": "in(corbital::Conjugate, co::ContractedOperator)\n\nTest if corbital is among the left set of orbitals of the ContractedOperator co. Useful to test if co is an integral operator with respect to corbital.\n\n\n\n\n\n"
},

{
    "location": "nbody_operators/#EnergyExpressions.contract",
    "page": "N-body operators",
    "title": "EnergyExpressions.contract",
    "category": "function",
    "text": "contract(orbital_matrix_element, i...)\n\nContract the orbital_matrix_element over all coordinates i....\n\n\n\n\n\ncontract(ome::OrbitalMatrixElement{N}, i...)\n\nContract ome over all coordinates i.... length(i) cannot be larger than N.\n\n\n\n\n\n"
},

{
    "location": "nbody_operators/#EnergyExpressions.complement",
    "page": "N-body operators",
    "title": "EnergyExpressions.complement",
    "category": "function",
    "text": "complement(N, i...)\n\nGenerate the complement to i... in the set 1:N. Useful for contracting OrbitalMatrixElements over all coordinates except i....\n\n\n\n\n\n"
},

{
    "location": "nbody_operators/#N-body-operators-1",
    "page": "N-body operators",
    "title": "N-body operators",
    "category": "section",
    "text": "CurrentModule = EnergyExpressions\nDocTestSetup = quote\n    using EnergyExpressions\n    using AtomicLevels\nendNBodyOperator\nLinearCombinationOperator\nIdentityOperator\nContractedOperator\nBase.in(orbital::O, co::ContractedOperator) where O\nBase.in(corbital::Conjugate{O}, co::ContractedOperator) where O\ncontract\ncomplement DocTestSetup = nothing"
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
    "location": "nbody_matrix_elements/#EnergyExpressions.numbodies",
    "page": "N-body matrix elements",
    "title": "EnergyExpressions.numbodies",
    "category": "function",
    "text": "numbodies(::NBodyOperator{N})\n\nReturns the number of bodies coupled by the N-body operator, i.e. N.\n\n\n\n\n\nnumbodies(lco::LinearCombinationOperator)\n\nReturns the maximum number of bodies coupled by any of the N-body operators in the LinearCombinationOperator.\n\n\n\n\n\nnumbodies(::OrbitalOverlap)\n\nReturns the number of bodies coupled by the zero-body operator in the orbital overlap, i.e. 0.\n\n\n\n\n\nnumbodies(::OrbitalMatrixElement{N})\n\nReturns the number of bodies coupled by the operator, i.e. N.\n\n\n\n\n\n"
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
    "location": "nbody_matrix_elements/#EnergyExpressions.isdependent",
    "page": "N-body matrix elements",
    "title": "EnergyExpressions.isdependent",
    "category": "function",
    "text": "isdependent(o::OrbitalOverlap, orbital)\n\nReturns true if the OrbitalOverlap o depends on orbital.\n\nExamples\n\njulia> isdependent(OrbitalOverlap(:a,:b), :a)\nfalse\n\njulia> isdependent(OrbitalOverlap(:a,:b), Conjugate(:a))\ntrue\n\njulia> isdependent(OrbitalOverlap(:a,:b), :b)\ntrue\n\n\n\n\n\nisdependent(o::OrbitalMatrixElement, orbital)\n\nReturns true if the OrbitalMatrixElement o depends on orbital.\n\nExamples\n\njulia> isdependent(EnergyExpressions.OrbitalMatrixElement((:a,), OneBodyHamiltonian(), (:b,)), :a)\nfalse\n\njulia> isdependent(EnergyExpressions.OrbitalMatrixElement((:a,), OneBodyHamiltonian(), (:b,)), Conjugate(:a))\ntrue\n\njulia> isdependent(EnergyExpressions.OrbitalMatrixElement((:a,), OneBodyHamiltonian(), (:b,)), :b)\ntrue\n\njulia> isdependent(EnergyExpressions.OrbitalMatrixElement((:a,:b,), CoulombInteraction(), (:c,:d)), :c)\ntrue\n\n\n\n\n\nisdependent(nbt::NBodyTerm, o)\n\nReturns true if any of the factors comprising nbt is dependent on the orbital o. Not that the result is dependent on whether o is conjugated or not.\n\n\n\n\n\n"
},

{
    "location": "nbody_matrix_elements/#EnergyExpressions.transform",
    "page": "N-body matrix elements",
    "title": "EnergyExpressions.transform",
    "category": "function",
    "text": "transform(f::Function, nbt::NBodyTerm)\n\nTransform integrals of the the N-body matrix element expansion term nbt according to the function f, which should accept a single NBodyTermFactor as its argument.\n\n\n\n\n\ntransform(f::Function, nbme::NBodyMatrixElement)\n\nTransform integrals of the the N-body matrix element nbme according to the function f, which should accept a single NBodyTermFactor as its argument, and return a NBodyMatrixElement. This is useful for adapting energy expressions to specific symmetries of the system under consideration.\n\n\n\n\n\n"
},

{
    "location": "nbody_matrix_elements/#EnergyExpressions.overlap_matrix",
    "page": "N-body matrix elements",
    "title": "EnergyExpressions.overlap_matrix",
    "category": "function",
    "text": "overlap_matrix(a::SlaterDeterminant, b::SlaterDeterminant[, overlaps=[]])\n\nGenerate the single-particle orbital overlap matrix, between the orbitals in the Slater determinants a and b. All orbitals are assumed to be orthogonal, except for those which are given in overlaps.\n\nExamples\n\nFirst we define two Slater determinants that have some orbitals in common:\n\njulia> sa = SlaterDeterminant([:i, :j, :l,:k̃])\ni(1)j(2)l(3)k̃(4) - i(1)j(2)l(4)k̃(3) - i(1)j(3)l(2)k̃(4) + i(1)j(3)l(4)k̃(2) + …  + i(4)j(1)l(3)k̃(2) + i(4)j(2)l(1)k̃(3) - i(4)j(2)l(3)k̃(1) - i(4)j(3)l(1)k̃(2) + i(4)j(3)l(2)k̃(1)\n\njulia> sb = SlaterDeterminant([:i, :j, :k, :l̃])\ni(1)j(2)k(3)l̃(4) - i(1)j(2)k(4)l̃(3) - i(1)j(3)k(2)l̃(4) + i(1)j(3)k(4)l̃(2) + …  + i(4)j(1)k(3)l̃(2) + i(4)j(2)k(1)l̃(3) - i(4)j(2)k(3)l̃(1) - i(4)j(3)k(1)l̃(2) + i(4)j(3)k(2)l̃(1)\n\nThe orbital overlap matrix by default is\n\njulia> overlap_matrix(sa, sb)\n4×4 SparseArrays.SparseMatrixCSC{EnergyExpressions.NBodyTerm,Int64} with 2 stored entries:\n  [1, 1]  =  1\n  [2, 2]  =  1\n\nwhich has only two non-zero entries, since only two of the orbitals are common between the Slater determinants sa and sb.\n\nWe can then define that the orbitals k̃ and l̃ are non-orthogonal:\n\njulia> overlap_matrix(sa, sb, [OrbitalOverlap(:k̃,:l̃)])\n4×4 SparseArrays.SparseMatrixCSC{EnergyExpressions.NBodyTerm,Int64} with 3 stored entries:\n  [1, 1]  =  1\n  [2, 2]  =  1\n  [4, 4]  =  ⟨k̃|l̃⟩\n\nWe can even specify that the orbital k̃ is non-orthogonal to itself (this can be useful when the k̃ is a linear combination of orthogonal orbitals):\n\njulia> overlap_matrix(sa, sa, [OrbitalOverlap(:k̃,:k̃)])\n4×4 SparseArrays.SparseMatrixCSC{EnergyExpressions.NBodyTerm,Int64} with 4 stored entries:\n  [1, 1]  =  1\n  [2, 2]  =  1\n  [3, 3]  =  1\n  [4, 4]  =  ⟨k̃|k̃⟩\n\nNotice that this overlap matrix was calculated between the Slater determinant sa and itself.\n\n\n\n\n\n"
},

{
    "location": "nbody_matrix_elements/#EnergyExpressions.EnergyExpression",
    "page": "N-body matrix elements",
    "title": "EnergyExpressions.EnergyExpression",
    "category": "type",
    "text": "EnergyExpression\n\nAn energy expression is given by an energy matrix, or interaction matrix, sandwiched between a vector of mixing coefficients: E = c\'H*c, where c are the mixing coefficients and H the energy matrix.\n\n\n\n\n\n"
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
    "text": "CurrentModule = EnergyExpressions\nDocTestSetup = quote\n    using EnergyExpressions\n    using AtomicLevels\nendThe matrix element of an N-body operator between two Slater determinants may be expanded according to the Löwdin rules (which reduce to the Slater–Condon rules if all single-particle orbitals are orthogonal):beginequation\nlabeleqnmatrix-element-expansion\nmatrixelPhi_AOmega_nPhi_B =\nfrac1nsum_p (-)^p\nmatrixelk_1k_2k_nOmega_nl_1l_2l_n\nD^AB(k_1k_2k_nl_1l_2l_n)\nendequationwhere D^AB(k_1k_2k_nl_1l_2l_n) is the determinant minor of the orbital overlap determinant D^AB with the rows k_1k_2k_n and columns l_1l_2l_n stricken out, and p runs over all permutations.In general, a term in the expansion is thus of the formbeginequation\nalphamatrixelk_1k_2k_nOmega_nl_1l_2l_nbraketabbraketcddotsbraketyz\nendequationwhere alpha is a scalar. This is represented by NBodyTerm type.NBodyTermFactor\nOrbitalOverlap\nOrbitalMatrixElement\nnumbodies\nNBodyTerm\nNBodyMatrixElement\nisdependent\ntransform\noverlap_matrix\nEnergyExpression\nMatrix"
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
    "text": "detaxis(i::CartesianIndex{N})\n\nGenerate the axis index vector for the determinant minor, whose rows or columns represented by the CartesianIndex i should be omitted. Implemented via complement.\n\n\n\n\n\n"
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
    "location": "nbody_matrix_elements/#Calculation-of-determinants-1",
    "page": "N-body matrix elements",
    "title": "Calculation of determinants",
    "category": "section",
    "text": "Actually computing the matrix element expansion eqrefeqnmatrix-element-expansion is a combinatorial problem, that grows factorially with the amount of non-orthogonal orbital pairs. Furthermore, of the (n)^2 terms generated from the expansion, only n are distinct, due to the integrals being symmetric with respect to interchange of the coordinates [hence the normalization factor (n)^-1]. Thankfully, there are few symmetries that can be employed, to generate only the distinct permutations.We use Julia\'s built in CartesianIndex iterators to span the space of all possible choices of orbitals for the overlap determinant. If two or more indices in the CartesianIndex are the same, the overlap is trivially zero (the “Fermi hole”). To avoid double-counting, we also only consider those indices that are above the hyper-diagonal.isabovediagonal\nisdiagonal\ndetaxis\ndetminor\n# indexsum\ncofactor\ndet\npermutation_sign DocTestSetup = nothing"
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
    "text": "OneBodyHamiltonian\n\nThe one-body Hamiltonian, may include external fields. It is diagonal in spin, i.e. it does not couple orbitals of opposite spin.\n\n\n\n\n\n"
},

{
    "location": "common_operators/#EnergyExpressions.FieldFreeOneBodyHamiltonian",
    "page": "Common N-body operators",
    "title": "EnergyExpressions.FieldFreeOneBodyHamiltonian",
    "category": "type",
    "text": "FieldFreeOneBodyHamiltonian\n\nThe one-body Hamiltonian, with no external fields. It is diagonal in the orbitals, i.e. does not couple unequal orbitals.\n\nExamples\n\njulia> EnergyExpressions.OrbitalMatrixElement((:a,), FieldFreeOneBodyHamiltonian(), (:a,))\n⟨a|ĥ₀|a⟩\n\njulia> iszero(EnergyExpressions.OrbitalMatrixElement((:a,), FieldFreeOneBodyHamiltonian(), (:a,)))\nfalse\n\njulia> EnergyExpressions.OrbitalMatrixElement((:a,), FieldFreeOneBodyHamiltonian(), (:b,))\n⟨a|ĥ₀|b⟩\n\njulia> iszero(EnergyExpressions.OrbitalMatrixElement((:a,), FieldFreeOneBodyHamiltonian(), (:b,)))\ntrue\n\n\n\n\n\n"
},

{
    "location": "common_operators/#EnergyExpressions.CoulombInteraction",
    "page": "Common N-body operators",
    "title": "EnergyExpressions.CoulombInteraction",
    "category": "type",
    "text": "CoulombInteraction\n\nTwo-body Hamiltonian, representing the mutual Coulombic repulsion between two electrons. Is diagonal in spin, i.e. the spin of the orbitals associated with the same coordinate must be the same.\n\nExamples\n\njulia> EnergyExpressions.OrbitalMatrixElement((:a,:b), CoulombInteraction(), (:c,:d))\n[a b|c d]\n\njulia> EnergyExpressions.OrbitalMatrixElement((:a,:b), CoulombInteraction(), (:b,:a))\nG(a,b)\n\n\n\n\n\n"
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
    "text": "DocTestSetup = quote\n    using EnergyExpressions\n    using AtomicLevels\nendA few N-body operators that are common in quantum mechanics. It is possible to form a composite Hamiltonian thus:julia> H = OneBodyHamiltonian() + CoulombInteraction()\nĥ + ĝIf then define a set of Slater determinants, we can easily form the energy expression:julia> slaters = SlaterDeterminant.([[:a, :b], [:c, :d], [:a, :c]])\n3-element Array{SlaterDeterminant{Symbol},1}:\n |a b|\n |c d|\n |a c|\n\njulia> Matrix(H, slaters)\n3×3 Array{EnergyExpressions.NBodyMatrixElement,2}:\n (a|a) + (b|b) - G(a,b) + F(a,b)  - [a b|d c] + [a b|c d]          (b|c) - [a b|c a] + [a b|a c]\n - [c d|b a] + [c d|a b]          (c|c) + (d|d) - G(c,d) + F(c,d)  - (d|a) - [c d|c a] + [c d|a c]\n (c|b) - [a c|b a] + [a c|a b]    - (a|d) - [a c|d c] + [a c|c d]  (a|a) + (c|c) - G(a,c) + F(a,c)An energy expression like this can then be used to derive the multi-configurational Hartree–Fock equations for the orbitals a,b,c,d.We can also specify that e.g. the orbitals b and c are non-orthogonal, and thus derive a slightly different energy expression, that takes this into account:julia> Matrix(H, slaters, [OrbitalOverlap(:b,:c)])\n3×3 Array{EnergyExpressions.NBodyMatrixElement,2}:\n (a|a) + (b|b) - G(a,b) + F(a,b)             - ⟨b|c⟩(a|d) - [a b|d c] + [a b|c d]  ⟨b|c⟩(a|a) + (b|c) - [a b|c a] + [a b|a c]\n - ⟨c|b⟩(d|a) - [c d|b a] + [c d|a b]        (c|c) + (d|d) - G(c,d) + F(c,d)       - (d|a) - [c d|c a] + [c d|a c]\n ⟨c|b⟩(a|a) + (c|b) - [a c|b a] + [a c|a b]  - (a|d) - [a c|d c] + [a c|c d]       (a|a) + (c|c) - G(a,c) + F(a,c)Beware that the computational complexity grows factorially with the amount of non-orthogonal orbitals!OneBodyHamiltonian\nFieldFreeOneBodyHamiltonian\nCoulombInteraction\niszero DocTestSetup = nothing"
},

{
    "location": "equations/#",
    "page": "N-body equations",
    "title": "N-body equations",
    "category": "page",
    "text": ""
},

{
    "location": "equations/#EnergyExpressions.NBodyEquation",
    "page": "N-body equations",
    "title": "EnergyExpressions.NBodyEquation",
    "category": "type",
    "text": "NBodyEquation{N,O}(orbital, operator::NBodyOperator[, factor::NBodyTerm])\n\nEquation for an orbital, acted upon by an operator, which may be a single-particle operator, or an N-body operator, contracted over all coordinates but one, and optionally multiplied by an NBodyTerm, corresponding to overlaps/matrix elements between other orbitals.\n\n\n\n\n\n"
},

{
    "location": "equations/#EnergyExpressions.LinearCombinationEquation",
    "page": "N-body equations",
    "title": "EnergyExpressions.LinearCombinationEquation",
    "category": "type",
    "text": "LinearCombinationEquation(equations)\n\nA type representing a linear combination of NBodyEquations. Typically arises when varying a multi-term energy expression.\n\n\n\n\n\n"
},

{
    "location": "equations/#N-body-equations-1",
    "page": "N-body equations",
    "title": "N-body equations",
    "category": "section",
    "text": "DocTestSetup = quote\n    using EnergyExpressions\n    using AtomicLevels\nendThese types are used to represent the integro-differential equations that result from the variation of the energy expression with respect to an orbital. They can, in general, be written on the formbeginequation\n0 = Omega_1ketchibraketab\nendequationwhen varying with respect to a conjugated orbital, andbeginequation\n0 = brachiOmega_1braketab\nendequationwhen varying with respect to an orbital. In both cases, Omega_1 is a one-body operator, either in itself, or resulting from a contraction over all coordinates but one of a many-body operator (see ContractedOperator).NBodyEquation\nLinearCombinationEquation DocTestSetup = nothing"
},

{
    "location": "variations/#",
    "page": "Variation",
    "title": "Variation",
    "category": "page",
    "text": ""
},

{
    "location": "variations/#Base.diff",
    "page": "Variation",
    "title": "Base.diff",
    "category": "function",
    "text": "diff(ab::OrbitalOverlap, o::O)\n\nVary the orbital overlap ⟨a|b⟩ with respect to |o⟩.\n\n\n\n\n\ndiff(ab::OrbitalOverlap, o::Conjugate{O})\n\nVary the orbital overlap ⟨a|b⟩ with respect to ⟨o|.\n\n\n\n\n\ndiff(ome::OrbitalMatrixElement, o::O)\n\nVary the orbital overlap ⟨abc...|Ω|xyz...⟩ with respect to |o⟩.\n\n\n\n\n\ndiff(ome::OrbitalMatrixElement, o::Conjugate{O})\n\nVary the orbital overlap ⟨abc...|Ω|xyz...⟩ with respect to ⟨o|.\n\n\n\n\n\ndiff(me::NBodyMatrixElement, o::O)\n\nVary the NBodyMatrixElement me with respect to the orbital o.\n\n\n\n\n\ndiff(E::Matrix{NBodyMatrixElement}, o::O)\n\nVary the matrix of NBodyMatrixElements with respect to the orbital o.\n\nExamples\n\njulia> E = Matrix(OneBodyHamiltonian()+CoulombInteraction(),\n                  SlaterDeterminant.([[:a, :b], [:c, :d]]))\n2×2 Array{EnergyExpressions.NBodyMatrixElement,2}:\n (a|a) + (b|b) - G(a,b) + F(a,b)  - [a b|d c] + [a b|c d]\n - [c d|b a] + [c d|a b]          (c|c) + (d|d) - G(c,d) + F(c,d)\n\njulia> diff(E, :a)\n2×2 SparseArrays.SparseMatrixCSC{LinearCombinationEquation,Int64} with 2 stored entries:\n  [1, 1]  =  ⟨a|ĥ + -⟨b|[a|b] + ⟨a|[b|b]\n  [2, 1]  =  -⟨d|[c|b] + ⟨c|[d|b]\n\njulia> diff(E, Conjugate(:b))\n2×2 SparseArrays.SparseMatrixCSC{LinearCombinationEquation,Int64} with 2 stored entries:\n  [1, 1]  =  ĥ|b⟩ + -[a|b]|a⟩ + [a|a]|b⟩\n  [1, 2]  =  -[a|d]|c⟩ + [a|c]|d⟩\n\n\n\n\n\ndiff(energy_expression, orbitals)\n\nDerive the integro-differential equations for all orbitals, from energy_expression. Returns a MCEquationSystem, that gathers information on which integrals are common to all equations, for efficient equation solving.\n\n\n\n\n\n"
},

{
    "location": "variations/#Calculus-of-Variations-–-Implementation-1",
    "page": "Variation",
    "title": "Calculus of Variations – Implementation",
    "category": "section",
    "text": "DocTestSetup = quote\n    using EnergyExpressions\n    using AtomicLevels\nenddiff DocTestSetup = nothing"
},

{
    "location": "system_of_equations/#",
    "page": "System of equations",
    "title": "System of equations",
    "category": "page",
    "text": ""
},

{
    "location": "system_of_equations/#System-of-equations-1",
    "page": "System of equations",
    "title": "System of equations",
    "category": "section",
    "text": "When dealing with large energy expressions resulting from many configurations, when deriving the orbital equations (see N-body equations and Calculus of Variations – Implementation), many of the integrals involved will be shared between the different terms of the equations. It is therefore of interest to gather a list of integrals that can be calculated once at every iteration (in a self-consistent procedure for finding eigenstates, or in time-propagation), and whose values can then be reused when solving the individual equations. The routines described below, although primitive, aid in this effort.The idea is the following: With an energy expression on the formbeginequation\nE(vecPvecc) = fracvecc^H\nmatHveccvecc^Hvecc\nendequationwhere vecP is a set of orbitals, vecc is a vector of mixing coefficients, andbeginequation\nmatH_ij defd matrixelchi_iHamiltonianchi_j\nendequationall terms of the equations can be written on the formbeginequation\nlabeleqnequation-term\noperatorAketchi\nunderbraceleft\nsum_n alpha_n conjc_i_n c_j_n\nprod_k int_n_k\nright_textrmRank 0\nendequation(or on its dual form with conjugated orbitals) where operatorA is some one-body operator, alpha_n a coefficient due to the energy expression, conjc_i_n c_j_n the coefficient due to the multi-configurational expansion, and finally int_n_k all the extra integrals (which are necessarily zero-body operators) that may appear due to non-orthogonal orbitals (and which may be shared between many equations). The operator operatorA can have ranks 0–2; formally, it can only have rank 0 (multiplication by a scalar) or 2 (multiplication by a matrix). What we somewhat sloppily to refer by “rank 1” is an operator of rank 2, but which is diagonal in the underlying coordinate, such as a local potential.Thus, to efficiently perform an iteration, one would first compute all common integrals int_k, and then for every equation term of the form eqrefeqnequation-term, form the coefficient in the brakets, calculate the action of the operator operatorA on the orbital chi, multiplied the coefficient.There are a few important special cases of eqrefeqnequation-term:Field-free, one-body Hamiltonian, i.e. operatorA=hamiltonian_0. This is the contribution from the orbital ketchi to itself. Rank 2.\nOne-body Hamiltonian, including field, i.e. operatorA=hamiltonian. This is the contribution from another orbital ketchi to ketchi via some off-diagonal coupling, such as an external field (e.g. an electro–magnetic field). Rank 2.\nDirect interaction, i.e. operatorA=directkl, where two other orbitals ketchi_k and ketchi_l together form a potential acting on ketchi. “Rank 1”.\nExchange interaction, i.e. operatorA=exchangekl, where another orbital ketchi_k and ketchi together form a potential acting on a third orbital ketchi_l. Rank 2.\nSource term, i.e. a contribution that does not involve ketchi in any way. This term arises from other configurations in the multi-configurational expansion. Case 2. is also formulated in this way. Rank 0–2. If for some reason the source orbital and/or the operator acting on it is fixed, it may be possible to precompute the effect of the operator on the source orbital and reduce the computational complexity to a rank 0-equivalent operation.CurrentModule = EnergyExpressions\nDocTestSetup = quote\n    using EnergyExpressions\nend"
},

{
    "location": "system_of_equations/#Example-1",
    "page": "System of equations",
    "title": "Example",
    "category": "section",
    "text": "We start by defining a set of configurations and specifying that a few of the constituent orbitals are non-orthogonal:julia> cfgs = [[:i, :j, :l, :k̃], [:i, :j, :k, :l̃]]\n2-element Array{Array{Symbol,1},1}:\n [:i, :j, :l, :k̃]\n [:i, :j, :k, :l̃]\n\njulia> continua = [:k̃, :l̃]\n2-element Array{Symbol,1}:\n :k̃\n :l̃\n\njulia> overlaps = [OrbitalOverlap(i,j) for i in continua for j in continua]\n4-element Array{OrbitalOverlap{Symbol,Symbol},1}:\n ⟨k̃|k̃⟩\n ⟨k̃|l̃⟩\n ⟨l̃|k̃⟩\n ⟨l̃|l̃⟩We then set up the energy expression as before:julia> H = OneBodyHamiltonian() + CoulombInteraction()\nĥ + ĝ\n\njulia> E = Matrix(H, SlaterDeterminant.(cfgs), overlaps)\n2×2 Array{EnergyExpressions.NBodyMatrixElement,2}:\n (i|i)⟨k̃|k̃⟩ + (j|j)⟨k̃|k̃⟩ + (l|l)⟨k̃|k̃⟩ + (k̃|k̃) - G(i,j)⟨k̃|k̃⟩ + F(i,j)⟨k̃|k̃⟩ + … - G(i,k̃) + F(i,k̃) - G(j,k̃) + F(j,k̃) - G(l,k̃) + F(l,k̃)  …  (l|k)⟨k̃|l̃⟩ - [i l|k i]⟨k̃|l̃⟩ + [i l|i k]⟨k̃|l̃⟩ - [j l|k j]⟨k̃|l̃⟩ + [j l|j k]⟨k̃|l̃⟩ - [l k̃|l̃ k] + [l k̃|k l̃]\n (k|l)⟨l̃|k̃⟩ - [i k|l i]⟨l̃|k̃⟩ + [i k|i l]⟨l̃|k̃⟩ - [j k|l j]⟨l̃|k̃⟩ + [j k|j l]⟨l̃|k̃⟩ - [k l̃|k̃ l] + [k l̃|l k̃]                                     (i|i)⟨l̃|l̃⟩ + (j|j)⟨l̃|l̃⟩ + (k|k)⟨l̃|l̃⟩ + (l̃|l̃) - G(i,j)⟨l̃|l̃⟩ + F(i,j)⟨l̃|l̃⟩ + … - G(i,l̃) + F(i,l̃) - G(j,l̃) + F(j,l̃) - G(k,l̃) + F(k,l̃)Finally, we derive the coupled integro-differential equation system for the continuum orbitals k̃, l̃:julia> eqs = diff(E, Conjugate.(continua))\nEnergyExpressions.MCEquationSystem(EnergyExpressions.OrbitalEquation{Symbol,SparseArrays.SparseMatrixCSC{LinearCombinationEquation,Int64}}[OrbitalEquation(k̃):\n  [1, 1]  =  𝐈₁|k̃⟩(i|i) + 𝐈₁|k̃⟩(j|j) + 𝐈₁|k̃⟩(l|l) + ĥ|k̃⟩ + 𝐈₁|k̃⟩(- G(i,j)) + 𝐈₁|k̃⟩F(i,j) + 𝐈₁|k̃⟩(- G(i,l)) + 𝐈₁|k̃⟩F(i,l) + 𝐈₁|k̃⟩(- G(j,l)) + 𝐈₁|k̃⟩F(j,l) + -[i|k̃]|i⟩ + [i|i]|k̃⟩ + -[j|k̃]|j⟩ + [j|j]|k̃⟩ + -[l|k̃]|l⟩ + [l|l]|k̃⟩\n  [1, 2]  =  𝐈₁|l̃⟩(l|k) + 𝐈₁|l̃⟩(- [i l|k i]) + 𝐈₁|l̃⟩[i l|i k] + 𝐈₁|l̃⟩(- [j l|k j]) + 𝐈₁|l̃⟩[j l|j k] + -[l|l̃]|k⟩ + [l|k]|l̃⟩\n, OrbitalEquation(l̃):\n  [2, 1]  =  𝐈₁|k̃⟩(k|l) + 𝐈₁|k̃⟩(- [i k|l i]) + 𝐈₁|k̃⟩[i k|i l] + 𝐈₁|k̃⟩(- [j k|l j]) + 𝐈₁|k̃⟩[j k|j l] + -[k|k̃]|l⟩ + [k|l]|k̃⟩\n  [2, 2]  =  𝐈₁|l̃⟩(i|i) + 𝐈₁|l̃⟩(j|j) + 𝐈₁|l̃⟩(k|k) + ĥ|l̃⟩ + 𝐈₁|l̃⟩(- G(i,j)) + 𝐈₁|l̃⟩F(i,j) + 𝐈₁|l̃⟩(- G(i,k)) + 𝐈₁|l̃⟩F(i,k) + 𝐈₁|l̃⟩(- G(j,k)) + 𝐈₁|l̃⟩F(j,k) + -[i|l̃]|i⟩ + [i|i]|l̃⟩ + -[j|l̃]|j⟩ + [j|j]|l̃⟩ + -[k|l̃]|k⟩ + [k|k]|l̃⟩\n], Any[(i|i), 𝐈₁, (j|j), (l|l), ĥ, G(i,j), F(i,j), G(i,l), F(i,l), G(j,l)  …  [j k|l j], [j k|j l], [k|k̃], [k|l], (k|k), G(i,k), F(i,k), G(j,k), F(j,k), [k|k]])We can investigate the MCEquationSystem object eqs a bit. It consists of two coupled equations:julia> eqs.equations\n2-element Array{EnergyExpressions.OrbitalEquation{Symbol,SparseArrays.SparseMatrixCSC{LinearCombinationEquation,Int64}},1}:\n OrbitalEquation(k̃):\n  [1, 1]  =  𝐈₁|k̃⟩(i|i) + 𝐈₁|k̃⟩(j|j) + 𝐈₁|k̃⟩(l|l) + ĥ|k̃⟩ + 𝐈₁|k̃⟩(- G(i,j)) + 𝐈₁|k̃⟩F(i,j) + 𝐈₁|k̃⟩(- G(i,l)) + 𝐈₁|k̃⟩F(i,l) + 𝐈₁|k̃⟩(- G(j,l)) + 𝐈₁|k̃⟩F(j,l) + -[i|k̃]|i⟩ + [i|i]|k̃⟩ + -[j|k̃]|j⟩ + [j|j]|k̃⟩ + -[l|k̃]|l⟩ + [l|l]|k̃⟩\n  [1, 2]  =  𝐈₁|l̃⟩(l|k) + 𝐈₁|l̃⟩(- [i l|k i]) + 𝐈₁|l̃⟩[i l|i k] + 𝐈₁|l̃⟩(- [j l|k j]) + 𝐈₁|l̃⟩[j l|j k] + -[l|l̃]|k⟩ + [l|k]|l̃⟩\n\n OrbitalEquation(l̃):\n  [2, 1]  =  𝐈₁|k̃⟩(k|l) + 𝐈₁|k̃⟩(- [i k|l i]) + 𝐈₁|k̃⟩[i k|i l] + 𝐈₁|k̃⟩(- [j k|l j]) + 𝐈₁|k̃⟩[j k|j l] + -[k|k̃]|l⟩ + [k|l]|k̃⟩\n  [2, 2]  =  𝐈₁|l̃⟩(i|i) + 𝐈₁|l̃⟩(j|j) + 𝐈₁|l̃⟩(k|k) + ĥ|l̃⟩ + 𝐈₁|l̃⟩(- G(i,j)) + 𝐈₁|l̃⟩F(i,j) + 𝐈₁|l̃⟩(- G(i,k)) + 𝐈₁|l̃⟩F(i,k) + 𝐈₁|l̃⟩(- G(j,k)) + 𝐈₁|l̃⟩F(j,k) + -[i|l̃]|i⟩ + [i|i]|l̃⟩ + -[j|l̃]|j⟩ + [j|j]|l̃⟩ + -[k|l̃]|k⟩ + [k|k]|l̃⟩The first equation consists of the following terms:julia> eqs.equations[1].terms\n6-element Array{Pair{Int64,Array{EnergyExpressions.MCTerm,1}},1}:\n  0 => [MCTerm{Int64,IdentityOperator{1},Symbol}(1, 1, 1, 𝐈₁, :k̃, [1]), MCTerm{Int64,IdentityOperator{1},Symbol}(1, 1, 1, 𝐈₁, :k̃, [2]), MCTerm{Int64,IdentityOperator{1},Symbol}(1, 1, 1, 𝐈₁, :k̃, [3]), MCTerm{Int64,OneBodyHamiltonian,Symbol}(1, 1, 1, ĥ, :k̃, Int64[]), MCTerm{Int64,IdentityOperator{1},Symbol}(1, 1, -1, 𝐈₁, :k̃, [4]), MCTerm{Int64,IdentityOperator{1},Symbol}(1, 1, 1, 𝐈₁, :k̃, [5]), MCTerm{Int64,IdentityOperator{1},Symbol}(1, 1, -1, 𝐈₁, :k̃, [6]), MCTerm{Int64,IdentityOperator{1},Symbol}(1, 1, 1, 𝐈₁, :k̃, [7]), MCTerm{Int64,IdentityOperator{1},Symbol}(1, 1, -1, 𝐈₁, :k̃, [8]), MCTerm{Int64,IdentityOperator{1},Symbol}(1, 1, 1, 𝐈₁, :k̃, [9]), MCTerm{Int64,ContractedOperator{1,2,1,Symbol,CoulombInteraction,Symbol},Symbol}(1, 1, -1, [i|k̃], :i, Int64[]), MCTerm{Int64,ContractedOperator{1,2,1,Symbol,CoulombInteraction,Symbol},Symbol}(1, 1, -1, [j|k̃], :j, Int64[]), MCTerm{Int64,ContractedOperator{1,2,1,Symbol,CoulombInteraction,Symbol},Symbol}(1, 1, -1, [l|k̃], :l, Int64[]), MCTerm{Int64,IdentityOperator{1},Symbol}(1, 2, 1, 𝐈₁, :l̃, [13]), MCTerm{Int64,IdentityOperator{1},Symbol}(1, 2, -1, 𝐈₁, :l̃, [14]), MCTerm{Int64,IdentityOperator{1},Symbol}(1, 2, 1, 𝐈₁, :l̃, [15]), MCTerm{Int64,IdentityOperator{1},Symbol}(1, 2, -1, 𝐈₁, :l̃, [16]), MCTerm{Int64,IdentityOperator{1},Symbol}(1, 2, 1, 𝐈₁, :l̃, [17])]\n 10 => [MCTerm{Int64,ContractedOperator{1,2,1,Symbol,CoulombInteraction,Symbol},Symbol}(1, 1, 1, [i|i], :k̃, Int64[])]\n 19 => [MCTerm{Int64,ContractedOperator{1,2,1,Symbol,CoulombInteraction,Symbol},Symbol}(1, 2, 1, [l|k], :l̃, Int64[])]\n 11 => [MCTerm{Int64,ContractedOperator{1,2,1,Symbol,CoulombInteraction,Symbol},Symbol}(1, 1, 1, [j|j], :k̃, Int64[])]\n 12 => [MCTerm{Int64,ContractedOperator{1,2,1,Symbol,CoulombInteraction,Symbol},Symbol}(1, 1, 1, [l|l], :k̃, Int64[])]\n 18 => [MCTerm{Int64,ContractedOperator{1,2,1,Symbol,CoulombInteraction,Symbol},Symbol}(1, 2, -1, [l|l̃], :k, Int64[])]where the MCTerm objects indicate which components of the mixing coefficient vector vecc need to be multiplied, and all the integers are pointers to the list of common integrals:julia> eqs.integrals\n32-element Array{Any,1}:\n (i|i)\n (j|j)\n (l|l)\n G(i,j)\n F(i,j)\n G(i,l)\n F(i,l)\n G(j,l)\n F(j,l)\n [i|i]\n [j|j]\n [l|l]\n (l|k)\n [i l|k i]\n [i l|i k]\n [j l|k j]\n [j l|j k]\n [l|l̃]\n [l|k]\n (k|l)\n [i k|l i]\n [i k|i l]\n [j k|l j]\n [j k|j l]\n [k|k̃]\n [k|l]\n (k|k)\n G(i,k)\n F(i,k)\n G(j,k)\n F(j,k)\n [k|k]From this we see that the 𝐈₁ (one-body identity operator) contribution to |k̃⟩ can be written as a linear combination of |k̃⟩ and |l̃⟩, weighted by different components of the mixing coefficient vector vecc and various other integrals. This is all the information necessary to set up an efficient equation solver."
},

{
    "location": "system_of_equations/#EnergyExpressions.MCTerm",
    "page": "System of equations",
    "title": "EnergyExpressions.MCTerm",
    "category": "type",
    "text": "MCTerm(i, j, coeff, operator, source_orbital, integrals=[])\n\nRepresents one term in the multi-configurational expansion. i and j are indices in the mixing-coefficient vector c (which is subject to optimization, and thus has to be referred to), coeff is an additional coefficient, and integrals is a list of indices into the vector of common integrals, the values of which should be multiplied to form the overall coefficient.\n\n\n\n\n\n"
},

{
    "location": "system_of_equations/#EnergyExpressions.OrbitalEquation",
    "page": "System of equations",
    "title": "EnergyExpressions.OrbitalEquation",
    "category": "type",
    "text": "OrbitalEquation(orbital, equation,\n                one_body, direct_terms, exchange_terms, source_terms)\n\nRepresents the integro-differential equation for orbital, expressed as a linear combination of the different terms, with pointers to the list of common integrals that is stored by the encompassing MCEquationSystem object.\n\n\n\n\n\n"
},

{
    "location": "system_of_equations/#EnergyExpressions.orbital_equation",
    "page": "System of equations",
    "title": "EnergyExpressions.orbital_equation",
    "category": "function",
    "text": "orbital_equation(E::EnergyExpression, orbital, integrals::Vector)\n\nGenerate the OrbitalEquation governing orbital by varying the EnergyExpression E, and storing common expressions in integrals.\n\n\n\n\n\n"
},

{
    "location": "system_of_equations/#EnergyExpressions.MCEquationSystem",
    "page": "System of equations",
    "title": "EnergyExpressions.MCEquationSystem",
    "category": "type",
    "text": "MCEquationSystem(equations, integrals)\n\nRepresents a coupled system of integro-differential equations, resulting from the variation of a multi-configurational EnergyExpression, with respect to all constituent orbitals. All integrals that are in common between the equations need only be computed once per iteration, for efficiency.\n\n\n\n\n\n"
},

{
    "location": "system_of_equations/#EnergyExpressions.pushifmissing!",
    "page": "System of equations",
    "title": "EnergyExpressions.pushifmissing!",
    "category": "function",
    "text": "pushifmissing!(vector, element)\n\nPush element to the end of vector, if not already present. Returns the index of element in vector.\n\n\n\n\n\n"
},

{
    "location": "system_of_equations/#Implementation-1",
    "page": "System of equations",
    "title": "Implementation",
    "category": "section",
    "text": "MCTerm\nOrbitalEquation\norbital_equation\nMCEquationSystem\npushifmissing!CurrentModule = nothing\nDocTestSetup = nothing"
},

]}
