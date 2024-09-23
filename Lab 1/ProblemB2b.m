n = @(t) r(t) + r(-t+2);
plot(t,n(t));
xlabel('t');
ylabel('n(t) = r(t) + r(-t+2)');
grid;
axis([-1 2 -.1 1.1])