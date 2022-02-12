# IGSE
Improved General Steinmetz Equation

This script implements the improved general Steinmetz equation core loss calculation method. For more information on
the techincal details, please see DOI: 10.1109/CIPE.2002.1196712.

Inputs:
    time : evenly sampled time vector in seconds
    flux : flux density vector in seconds
    a    : Steinmetz coefficient alpha
    b    : Steinmetz coefficient beta
    k    : Steinmetz coefficient k
    tol  : flux periocity check tolerance (optional, 0.1)
Output:
    loss : core loss in W/m^3


The optional argument "tol" is the tolerance used when checking flux periodicity. Ideally, the input should be a
single cycle of period signal and flux(1) == flux(end). However, if the input is from simulation, it is more often
than not that this is not true. This argument suggest how different is tolerable, with 0 being the most strict and 1
being the most loose.

The best execution speed is about 60x-120x faster than previous implementation. For best execution speed, provide a
over-sampled single cycle of time and flux signal with even time spacing. This is usually true for simulation results
coming from MATLAB/Simulink. For other simulators, consider resample the signal in MATLAB first.

Caveat: the implementation here does NOT perform interpolation at minor loop boundaries for the sake of speed. For
most cases where the signal is sufficiently sampled like in simulation, there isn't any issue. However, if the input
is not over-sampled, the error may be significant.
