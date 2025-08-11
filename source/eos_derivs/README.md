# EOS Derivatives

This tool solves the EOS and calcs derivatives and other quantities
Takes in the density and pressure currently but can be changed to the 
density and temperature as well.

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
./feosderiv.gnu.ex diag.plotfile=plt00000
```
