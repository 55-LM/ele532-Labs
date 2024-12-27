% L2_A2.m : Chapter 2, MATLAB Program 2
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

%Plotting the impulse response of the system for t=0 to t=0.1 with 0.0005
%increments
t = 0:0.0005:0.1;
h = @(t) -0.0045*exp(lambda(1).*t) + 0.0045*exp(lambda(2).*t);

plot(t,h(t))
title('Impulse Response of the System');
xlabel('Time (s)');
ylabel('h(t)');
grid on;