# Problem 8.10, Sethna, on stochasticity in cell dynamics.

# Important note: The documentation of the implemented functions as well as the code used can be found in the notebook of this repository.

Stochastic cells (Biology, Computation)  Living cells are amazingly complex mixtures of a variety of complex molecules (RNA, DNA, proteins, lipids, . . . ) that are constantly undergoing
reactions with one another. This complex of reactions has been compared to computation; the cell gets input from external and internal sensors, and through an intricate series of reactions produces an appropriate response. Thus, for example, receptor cells in the retina ‘listen’ for light and respond by triggering a nerve impulse. The kinetics of chemical reactions are usually described using differential equations for the concentrations of the various chemicals, and rarely are statistical fluctuations considered important. In a cell, the numbers of molecules of a given type can be rather small; indeed, there is (often) only one copy of the relevant part of DNA for a given reaction. It is an important question whether and when we may describe the dynamics inside the cell using continuous concentration variables, even though the actual numbers of molecules are always integers.

<p align="center">
  <img src="https://raw.githubusercontent.com/EstebanM-98/Project_Stochastic_cells/refs/heads/main/Images/project.png" alt="Imagen" width="400">
</p>


Consider a dimerization reaction; a molecule M (called the ‘monomer’) joins up with another monomer and becomes a dimer D: 2M $ \leftrightarrow $ D. Proteins in cells often form dimers; sometimes (as here) both proteins are the same (homodimers) and sometimes they are different proteins (heterodimers). Suppose the forward reaction rate is $k_d$ and the backward reaction rate is $k_u$. Figure 8.11 shows this as a Petri net  with each reaction shown as a box, with incoming arrows showing species that are consumed by the reaction, and outgoing arrows showing species that are produced by the reaction; the number consumed or produced (the stoichiometry) is given by a label on each arrow. There are thus two reactions: the backward unbinding reaction rate per unit volume is $k_u$ [D] (each dimer disassociates with rate ku), and the forward binding reaction rate per unit volume is $k_b M^2$ (since each monomer must wait for a collision with another monomer before binding, the rate is proportional to the monomer concentration squared).

The brackets [.] denote concentrations. We assume that the volume per cell is such that one molecule per cell is 1 nM ($10^{−9}$ moles per liter). For convenience, we shall pick nanomoles as our unit of concentration, so [M] is also the number of monomers in the cell. Assume $k_b$ =1 $nM^{−1}s^{-1}$ and $k_u $= 2 $s^{-1}$, and that at t = 0 all N monomers are unbound.

**(a) Continuum dimerization. Write the differential equation for dM/dt treating M and D as continuous variables. (Hint: Remember that two M molecules are consumed in each reaction.) What are the equilibrium concentrations for [M] and [D] for N = 2 molecules in the cell, assuming these continuous equations and the values above for $k_b$ and $k_u$? For N = 90 and N = 10100 molecules?**


We set up differential equations

$$\frac{d[M]}{dt} = 2k_b[D] - 2k_b[M]^2, \quad \frac{d[D]}{dt} = k_b[M]^2 - k_u[D]$$

then, we use the conservation number of molecules $[M] + 2[D] = N$. Substituting out $[D]$ gives

$$ \frac{d[M]}{dt} = k_uN - k_u[M] - 2k_b[M]^2.$$

Setting $\frac{d[M]}{dt} = 0$ for equilibrium, we solve the quadratic equation, keeping only the positive solution,

$$[M]_0 = -\frac{k_u}{4k_b} + \sqrt{\left(\frac{k_u}{4k_b}\right)^2 + \frac{k_u}{2k_b}N} $$

For $k_u = 2$ $\mathrm{s}^{-1}$ and $k_b = 1$ $\mathrm{nM}^{-1}\mathrm{s}^{-1}$, we have $k_u/4k_b = \frac{1}{2}$ $\mathrm{nM}$. For $N = 2$ $\mathrm{nM}$ we find $[M]_0 = 1$ $\mathrm{nM}$. For $N = 90$ $\mathrm{nM}$ we get $[M]_0 = 9$ $\mathrm{nM}$. Finally for $N = 10100$ $\mathrm{nM}$ we have $[M]_0 = 100$ $\mathrm{nM}$.

**Numerically solve your differential equation for M(t) for N = 2 and N = 90, and verify that your solution settles down to the equilibrium values you found.**

The function solve_ivp of scipy was used for the solution of the differential equation. It is taken as an initial condition that the number of monomers at time $t = 0$ is equal to N.

<p align="center">
  <img src="https://raw.githubusercontent.com/EstebanM-98/Project_Stochastic_cells/refs/heads/main/Images/Nvst_continuum.png" alt="Imagen" width="600">
</p>


In the large time state the solution converges to the steady state.

# Monte Carlo Algorithm.

**For large numbers of molecules in the cell, we expect that the continuum equations may work well, but for just a few molecules there surely will be relatively large fluctuations. These fluctuations are called shot noise, named in early studies of electrical noise at low currents due to individual electrons in a resistor. We can implement a Monte Carlo algorithm to simulate this shot noise [36]. Suppose the reactions have rates $\Gamma_i$, with total rate $\Gamma_{tot}$ = $\sum_i \Gamma_i$. The idea is that the expected time to the next reaction is $1/\Gamma_{tot}$ , and the probability that the next reaction will be j is $\Gamma_j/\Gamma_{tot}$. To simulate until a final time $t_f$ , the algorithm runs as follows. To simulate until a final time $t_f$ , the algorithm runs as follows.**

**(1) Calculate a list of the rates of all reactions
in the system.**


**(2) Find the total rate Γtot.**

**(3) Pick a random time twait with probability distribution $\rho(t) = \Gamma_{tot} \times \text{exp}(− \Gamma_{tot} t)$.**

**(4) If the current time t plus twait is bigger than tf, no further reactions will take place; return.**

**(5) Otherwise,**

**– increment t by twait,**

**– pick a random number r uniformly distributed in the range [0, Γtot), $\sum_{i < j} \Gamma_i \leq r < \sum_{i < j+1} \Gamma_i$  (that is, $r$ lands in the $j$ th interval of the sum forming $\Gamma_{\text{tot}}$),**

**– execute that reaction, by incrementing each chemical involved by its stoichiometry.**

**(6) Repeat.**

**There is one important additional change the reaction rate for M total monomers binding is no longer $k_bM^2$ for discrete molecules; it is $k_bM(M − 1)$**


**(b) Stochastic dimerization. Implement this algorithm for the dimerization reaction of part (a).
Simulate for N = 2, N = 90, and N = 10100 and compare a few stochastic realizations with
the continuum solution.**

The implementation of the algorithm follows step by step the above described and a function called gillespie is created. For step 3 the function np.random.exponential($1/\Gamma_{tot}$) is used from the numpy library, which, according to the documentation, generates the desired random number distribution. For the comparison with the continuous case, the change of the reaction rate for M is taken into account and a function called dim_reaction2 is created.

<p align="center">
  <img src="https://raw.githubusercontent.com/EstebanM-98/Project_Stochastic_cells/refs/heads/main/Images/Nvst_stoc_comp.png" alt="Imagen" width="600">
</p>

In the case of a small number of monomers the fluctuations are evident with respect to the exact computational solution.

**How large a value of N do you need for the individual reactions to be well described by the continuum equations (say, fluctuations less than ±20% at late times)?**

We take two things into account. We will consider large times as those from which the monomers are a fraction of the initial value, since we know that these have an asymptotic behavior in large times. On the other hand, from this value the percentage error with respect to the exact solution is calculated for each time instant and an average is calculated. This has been implemented with the function calculate_relative_percentage_error.

<p align="center">
  <img src="https://raw.githubusercontent.com/EstebanM-98/Project_Stochastic_cells/refs/heads/main/Images/Mean_percentage_error_vsN_without_many_realizatios.png" alt="Imagen" width="600">
</p>


According to this method of analysis, determining the exact N value from which the percentage error is below 20 % is considerably random, since it may be the case, as is evident in the plot, that small values (e.g. less than 200) have a mean percentage error greater or less than 20 %. Therefore, we consider N from which all values have a mean percentage error of less than 20 %, i.e. approximately N = 600. It should be noted that this value can certainly fluctuate due to the statistical nature of the simulation.

**Measuring the concentrations in a single cell is often a challenge. Experiments often average over many cells. Such experiments will measure a smooth time evolution even though the individual cells are noisy. Let us investigate whether this ensemble average is well described by the continuum equations.**

**(c) Average stochastic dimerization. Find the average of many realizations of your stochastic dimerization in part (b), for N = 2 and N = 90.**

For each value of N, we used 10000 realizations of the Monte Carlo simulation and then took the average. We see, as expected, that for N large, the simulation approximates the solution of the differential equation.

<p align="center">
  <img src="https://raw.githubusercontent.com/EstebanM-98/Project_Stochastic_cells/refs/heads/main/Images/Mvst_many_realizations_comparison.png" alt="Imagen" width="700">
</p>

**How large a value of N do you need for the ensemble average of M(t) to be well described by the continuum equations (say, shifted by less than 5% at late times)**


<p align="center">
  <img src="https://raw.githubusercontent.com/EstebanM-98/Project_Stochastic_cells/refs/heads/main/Images/Mean_percentage_error_vsN_with_many_realizatios.png" alt="Imagen" width="600">
</p>

We did the same as in b) for computing the percentage error. We found that approximately for N = 80 the ensemble average of M(t) is well described by the continuum equations.
