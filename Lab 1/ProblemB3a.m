n1 = @(t) n(0.5*t);
plot(t, n1(t));
xlabel('t');
ylabel('n1(t) = n(0.5*t)');
grid;
axis([-1 5 -.1 1.1])