f = @(t) exp(-t).*cos(2*pi*t);
t = (-2:2);
plot(t,f(t));
xlabel('t'); 
ylabel('f(t)')
grid;