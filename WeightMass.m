%%WeightMass
%Use the given equations to calculate the weight and mass of a sphere in
%fluid with the radius and density of the sphere and the density of the
%fluid.

%%Typically scripts begin with
%{
Rui Xu
005230642
December 14,2018
%}

function [W,m] = WeightMass(R,rho_metal,rho_fluid)
%W is the weight of a sphere in fluid
%m is the mass of a sphere in fluid
%R is the radius of the sphere
%rho_metal is the density of the sphere
%rho_fluid is the density of the fluid
g=9.8; %gravitational pull
W=4/3*pi*R^3*(rho_metal-rho_fluid)*g;
m=4/3*pi*R^3*rho_metal;
end

