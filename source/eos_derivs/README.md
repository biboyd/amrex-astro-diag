# Convective gradients

This tool constructs gradients to estimate work done by 
convection moving electrons up a chemical gradient.  It is built
off of convective_grad

The gradients/values that are evaluated are:
- `("boxlib", "eta_e")` The electron degeneracy parameter $\eta$ which is defined as below. $\mu$ is the electron chemical potential, $m_e$ is the electron mass, $k_B$ is the boltzmann constant. units dimensionless erg/erg:
$$\eta= (\mu - m_e c^2)/k_B T$$

- `("boxlib", "d_chem_e")` The radial gradient of the electron chemical $\mu$, as derived from $\eta$ above. Units = erg/cm:
$$ \frac{d \mu}{d r}$$

- `("boxlib", "flux_e")` The Radial Flux density of electrons. With $\mathrm{N_A}$ is Avogadro number ie number of nucleons per gram, $Y_e$ is electron fraction number of e- per nucleon, $U_r$ is the radial velocity.  units # e- /cm^2:
$$f_e = \mathrm{N_A} \rho Y_e U_r$$

- `("boxlib", "eps_conv")` A somewhat motivated quantity of estimating the work done by convection in moving the electrons from the edge to the center (ie up the chemical gradient). Using the flux value described above and the electron chemical potential gradient. units erg / cm^3:

$$\epsilon_{\mathrm{conv}} = f_e * \frac{d \mu}{dr}$$

Derivatives are constructed radially from x,y,z like so for $dT/dr$:

$$\frac{dT}{dr}  = \left( \frac{x}{r} \frac{dT}{dx} + \frac{y}{r} \frac{dT}{dy} + \frac{z}{r} \frac{dT}{dz}\right)$$


To build, do:

```
make
```
It is also important that the network you build with matches
the one used for generating the plotfile.  This is set via
the `NETWORK_DIR` parameter in the `GNUmakefile`.

Runtime parameters are managed by AMReX's ParmParse. To run,
you specify the plotfile via `diag.plotfile`, either in an inputs
file or on the command line, e.g.:

```
./fconvwork.gnu.ex diag.plotfile=plt00000
```
