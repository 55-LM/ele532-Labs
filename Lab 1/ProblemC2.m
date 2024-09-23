s = @(t) exp(-2)*exp(-2*t).*cos(4*pi*t).*u(t+1);
t = (-2:0.01:4);
plot(t, s(t));
xlabel('t');
ylabel('s(t) = g(t+1)');
grid;