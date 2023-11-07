function dndt = species_DE(t,n,A,r,Press_amount,Press_species)
% this function defines the generalised Lotka-Volterra equations for a
% specific network. This includes a press perturbation applied to a species
%
% % t = time variable for ODE system
% % n = current abundances (vertical vector)
% % A = interaction matrix (square matrix with length = length(n))
% % r = growth rate vector (vertical vector with length = length(n))
% % Press_amount = strength of a press perturbation applied to a species
% % Press_species = species which the perturbation is applied to

% if no perturbation is defined, then apply a perturbation with zero
% strength to the first species (i.e. no effect)
if nargin == 4
    Press_amount = 0;
    Press_species = 1;
end

% update ODE system based on the generalised Lotka-Volterra equations
dndt = n.*r + (A*n).*n;

% Apply the perturbation proportional to species' abundance
dndt(Press_species) = dndt(Press_species) - Press_amount.*n(Press_species);
