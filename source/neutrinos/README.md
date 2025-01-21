# Neutrino Loss Rates

This tool calculates the neutrino losses, both thermal and 
from weak reactions. Currently only the A=23 Urca pair of
reactions are included


```
make 
```

It is important that the network you build with matches
the one used for generating the plotfile.  This is set via
the `NETWORK_DIR` parameter in the `GNUmakefile`.
Currently the only considered rates are the A=23 Urca rates included
in `URCA-simple` and `URCA-medium`

Runtime parameters are managed by AMReX's ParmParse. To run,
you specify the plotfile via `diag.plotfile`, either in an inputs
file or on the command line, e.g.:

```
./fneutrinos.gnu.ex diag.plotfile=plt00000
```
