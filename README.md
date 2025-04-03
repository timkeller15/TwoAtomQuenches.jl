# TwoAtomQuenches.jl

This repository contains a Julia implementation of the code and data necessary for recreating the figures in the publication

*Probing the out-of-equilibrium dynamics of two interacting atoms*  
Tim Keller and Thomás Fogarty  
Phys. Rev. A **94**, 063620 (2016).  [doi: 10.1103/PhysRevA.94.063620](https://doi.org/10.1103/PhysRevA.94.063620)

## Physical Model and Quantities
We consider a system of two particles trapped in a one-dimensional harmonic potential and interacting through a point-like potential of strength $g$. At time $\tau=0$ the system is quenched by turning on $\delta$-shaped perturbations of height $\kappa$ at the harmonic trap center.

In dimensionless harmonic-oscillator units the Hamiltonian reads 

$$
\mathcal{H} =\sum_{j=1}^{2}  \left(  -\frac{1}{2}\frac{\partial^{2}}{\partial x _{j}^{2}}+\frac{1}{2}x_{j}^{2} + \kappa\delta (x_{j})\theta(\tau)\right) + g\delta (|x_{1}-x_{2}|) .
$$

The state of the system at some time $t>0$ after the quench is given by

$$
\Psi(x_{1},x_{2},t)=\sum_{n=0}^{\infty}a_{n}\phi_{n}(x_{1},x_{2})e^{-iE'_{n}t}
$$

where $\ket{\phi_n}$ and $E_n'$ are the eigenstates and eigenenergies of the perturbed system and the time-independent coefficients

$$
a_{n}=\int_{-\infty}^{+\infty}\int_{-\infty}^{+\infty}\phi_{n}^{*}(x_{1},x_{2})\psi_{0}(x_{1},x_{2})dx_{1}dx_{2}
$$

express the initial state $\ket{\psi_0}$ in the perturbed basis. 

We quantify the effect of the quench mainly by calculating the Loschmidt echo, which compares the perturbed and unperturbed evolution according to

$$
\mathcal{L}(t)=\left|\bra{\psi_0}e^{i\mathcal{H_f}t}e^{-i\mathcal{H_i}t}\ket{\psi_0}\right|^2=\left|\nu(t)\right|^2=\left|\sum_{n}\left|a_{n}\right|^{2}e^{i(E_0-E'_{n})t}\right|^{2} ,
$$

where $E_0$ is the energy of the unperturbed ground state. We also use the spectral function

$$
A(\omega)=2\pi \: \text{Re} \int_{-\infty}^{\infty}e^{i\omega t}\nu(t)dt = 2\pi\sum_{n}\left|a_{n}\right|^{2}\delta\left(\omega-\omega_{0}+\omega'_{n}\right) 
$$

to evaluate the response of the system. 

## Numerical Algorithm 
We calculate the energies and states of the perturbed basis by diagonalizing the Hamiltonian after the quench. The Hamiltonian is expressed on $N$ Cartesian position grid points $x_1,\ldots,x_N = -(N-1)/2,\ldots,(N-1)/2$ with Gauss quadrature basis functions

$$
f_{i}(x)=\frac{1}{N}\frac{\text{sin}[\pi(x-x_{i})]}{\text{sin}[\pi(x-x_{i})/N]}  .
$$

The kinetic energy term expressed in this basis reads

$$
T_{ij}=\begin{cases}\frac{\pi^2}{6}\left(1-\frac{1}{N^2}\right) &i=j \\
(-1)^{(i-j)}\frac{\pi^2}{N^2}\frac{\text{cos}[\pi(i-j)/N]}{\text{sin}^2[\pi(i-j)/N]} &i\neq j 
\end{cases}
$$

and the total Hamiltonian is described by

$$
\sum_{k,l=1}^{N}\left\lbrace\left(\frac{1}{h^{2}}T_{ik}+V(hx_{k})\delta_{ik}\right)\delta_{jl}+\left(\frac{1}{h^{2}}T_{jl}+V(hx_{l})\delta_{jl}\right)\delta_{ik}+\frac{g}{h}\delta_{il}\delta_{jl}\delta_{kl}\right\rbrace\Psi_{kl}=E\Psi_{ij} ,
$$

where the factor $h$ denotes the spacing between grid points after rescaling them to the region of interest. 

## Initial State

As the initial state we use the analytic solution for the ground state of the unperturbed system. 
By changing to center-of-mass and relative coordinates $X=\frac{1}{\sqrt{2}}(x_1+x_2)$ and $x=\frac{1}{\sqrt{2}}(x_1-x_2)$ the Hamiltonian reads

$$
\mathcal{H}=\left(-\frac{1}{2}\frac{\partial^2}{\partial X^2} + \frac{1}{2}X^2\right) + \left(-\frac{1}{2}\frac{\partial^2}{\partial x^2}+\frac{1}{2}x^2+\frac{g}{\sqrt{2}}\delta(x)\right)
$$

The odd eigenstates of the harmonic potential vanish at $x=0$ and therefore they are unaffected by the interaction. The even eigenstates of the relative motion are given by 

$$
\psi_{i}(x)=\mathcal{N}_{i}U\left(\frac{1}{4}-\frac{E_{i}}{2},\frac{1}{2},x^2\right)e^{-x^2/2}  ,  \qquad i=0,2,4,\ldots
$$

with the confluent hypergeometric function of the second kind $U(a,b,z)$ and a normalization constant $\mathcal{N}_i$. The corresponding eigenenergy is determined from the solutions of

$$
g =-2\sqrt{2}\frac{\Gamma\left(-\frac{E_i}{2}+\frac{3}{4}\right)}{\Gamma\left(-\frac{E_i}{2}+\frac{1}{4}\right)}   \qquad i=0,2,4,\ldots
$$

Overall, we find $\psi_0(x_1,x_2) \sim U\left(\frac{1}{4}-\frac{E_0}{2},\frac{1}{2},(x_1-x_2)^2/2\right)e^{-(x_1^2 + x_2^2)/2}$ for the unperturbed ground-state. 