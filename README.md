# Viscous Burger's Equation

Two solutions, written in MATLAB, for solving the viscous Burger's equation. They are both spectral methods: the first is a Fourier Galerkin method, and the second is Collocation on the Tchebyshev-Gauß-Lobatto points. 

##Nota Bene 

These codes can be easily adapted for just the heat equation or the inviscid Burger's equation. The former can be done unofficially by setting the diffusion coefficient to be very large (order unity) or officially by changing the the arrays associated with the non-linear term to be zeros within the code. The latter simply requires the diffusion coefficient to be set to 0. 