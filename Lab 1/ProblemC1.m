f = @(t) exp(-2*t).*cos(4*pi*t);
g = @(t) f(t).*u(t);
t = (-2:0.01:2);
plot(t,g(t));
xlabel('t');
ylabel('g(t) = f(t)*u(t)');
grid;