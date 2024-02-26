The codes first computes the dispersion "e(k)" of Eq.3 of the SYK notes taking in input the "alpha" parameter and the sizes of the lattice. From e(k) the density of states (DoS) is computed over an energy grid. This imporoves numerical stability and allows to consider a very fine k-mesh in a computationally cheap way. The DoS has to be computed only once at the beginning, note that, especially for small alpha, the number of sites has to be quite large in order to have a smooth profile. Also, the finer the energy gird mesh the higher the number of k-point has to be. Once the DoS is computed it doesn't have to be updated duting the scan over temperature.

For each temperature Eq.17 and Eq.18 are solved self-consistently, namely from an initial guess for the Green's function "G(iw)" (the non interacting one) the self-energy "S(iw)" is computed. A new G(iw) is then produced and so on an so forth until the realtive distance between two subsequent G(iw) is below an user-provided threshold.

Once the G(iw) is converged the total energy is computed with Eq.37. The chemical potential "mu" parameter is provided by the user.

The code scans the low-temperature limit between the "TMIN" and "TMAX" parameters with a fine mesh of NT points, and if the NT_INF parameters is different than 0 it uses other NT_INF points to reach a temperature which equals the bandwidth "WBAND" of the model defined by alpha. This is because in the high-T regime the specific heat usually does not deviate from a power-law behaviour.

The total energy per particle "E/N" of NT+NT_INF temperature points is then interpolated with a cubic spline over "NT_intp" points between the same boundaries. The specific heat "Cv/N = 1/N dE/dT" is computed with a mid-point derivative over the new dataset.

The entropy per particle "S(T)/N" is computed as: 

S(Tinf)/N - S(T)/N = 1/N int ^{Tinf} _{T}  Cv/T dT

the choice of boundaries is due to the fact that the incognita of the model is the zero-T limit of S so an estimate for the Tinf limit is required. The latter is desumed from the fact that at T=inf (which we consider to be reached when T=WBAND) all the microstates are thermodinamically available which, for a 1D chain of length N are only empty/full. So the entropy at infinite temperature will be given by:

S(Tinf) = N log(2)

The quantity to be studied is then:

S(T)/N = log(2) -  1/N int ^{Tinf} _{T}  Cv/T dT

From the plot we can see that S reaches finite values for every alpha < 3/2 while it starts to deviate for higher alphas. The Fermi-liquid regime, where S goes to zero at T=0 is then reached for alpha > 3.

The code provides the data on the interpolated fine energy grid in the DeltaS_alpha[**].DAT file.
col 1: temperature in Kelvin
col 2: total energy per particle
col 3: specific heat
col 4: entropy difference 1/N int ^{Tinf} _{T}  Cv/T dT

In order to visualize the absolute entropy in the T-> limit one has to plot col1 against log(2)-col4. 

For each temperature a folder is created where the initial (non interacting) G(iw) and the inital S(iw) are stored. Upon convergence also the final G(iw) and S(iw) are stored. If the VERBOSE parameter is set to "T" then also the G(iw) and S(iw) for each self-consistency loop are printed.



