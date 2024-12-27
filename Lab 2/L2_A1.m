% L2_A1.m : Chapter 2, MATLAB Program 1
% Script M-file determines characteristic roots of op-amp circuit.

% Set component values:
R = [1e4, 1e4, 1e4]; C = [1e-6, 1e-6];
% Determine the coefficients for characteristic equation:
A = [1, (1/R(1)+1/R(2)+1/R(3))/C(2), 1/(R(1)*R(2)*C(1)*C(2))];
% Determine characteristic roots:
lambda = roots(A);

%The poly() function returns the coeffcieints of a polynomial that has its
%roots specified.
poly(lambda);


