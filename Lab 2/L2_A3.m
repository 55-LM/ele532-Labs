function [lambda] = L2_A3(R,C)
% L2_A3.m : Chapter 2, MATLAB Program 3
% Function M-file finds characteristic roots of op-amp circuit.
% INPUTS:   R = length-3 vector of resistances
%           C = length-2 vector of capacitances% OUTPUTS:  
% lambda = characteristic roots
% Determine coefficients for characteristic equation:
A = [1, (1/R(1)+1/R(2)+1/R(3))/C(2), 1/(R(1)*R(2)*C(1)*C(2))];
% Determine characteristic roots:
lambda = roots(A);