# Spin-orbitals


$$\begin{equation}
\chi_a(\tau_i) \defd
\psi_a(\vec{r}_i)\sigma_a(s_i),
\end{equation}$$

where $\sigma(s)$ is either $\alpha$ (spin up) or $\beta$
(spin down).

# One-body matrix elements

$$\begin{equation}
I(a,b) \equiv \onebody{a}{b} \defd \matrixel{a}{\hamiltonian}{b},
\end{equation}$$

where

$$\begin{equation}
\hamiltonian \defd T + V
\end{equation}$$

is the (possibly time-dependent) one-body Hamiltonian.

# Two-body matrix elements

$$\begin{equation}
\twobodydx{ab}{cd} \defd
\twobody{ab}{cd} - \twobody{ab}{dc},
\end{equation}$$

where

$$\begin{equation}
\begin{aligned}
\twobody{ab}{cd} &\defd
\int\diff{\tau_1}\diff{\tau_2}
\conj{\chi_a}(\tau_1)
\conj{\chi_b}(\tau_2)
\frac{1}{r_{12}}
\chi_c(\tau_1)
\chi_d(\tau_2) \\
&=
\delta(\sigma_a,\sigma_c)
\delta(\sigma_b,\sigma_d)
\int\diff{\vec{r}_1}\diff{\vec{r}_2}
\conj{\chi_a}(\vec{r}_1)
\conj{\chi_b}(\vec{r}_2)
\frac{1}{r_{12}}
\chi_c(\vec{r}_1)
\chi_d(\vec{r}_2).
\end{aligned}
\end{equation}$$

The special case

$$\begin{equation}
F(a,b) \defd \twobody{ab}{ab}
\end{equation}$$

is called the *direct interaction* (gives rise to the screening
potential), and the other special case

$$\begin{equation}
G(a,b) \defd \twobody{ab}{ba}
\end{equation}$$

is called the *exchange interaction* (gives rise to the non-local
potential).

NB

$$\begin{equation}
\twobodydx{ii}{ii} = 0,
\end{equation}$$

such that

$$\begin{equation}
\frac{1}{2}\twobodydx{ij}{ij} \equiv
\sum_{i,j<i}\twobodydx{ij}{ij},
\end{equation}$$

i.e. we sum over $i,j$, divide by two to avoid double-counting
and avoid automatically the case $i=j$.

# Repulsion potential

We also have the repulsion potential formed from two orbitals $a,b$:

$$\begin{equation}
\twobody{a}{b} \defd
\int\diff{\tau_1}\conj{\chi_a}(\tau_1)\frac{1}{r_{12}}\chi_b(\tau_1)=
\delta(\sigma_a,\sigma_b)\int\diff{\vec{r}_1}\conj{\chi_a}(\vec{r}_1)\frac{1}{r_{12}}\chi_b(\vec{r}_1).
\end{equation}$$

