n2 = @(t) n1(t+0.5);
plot(t, n2(t));
xlabel('t');
ylabel('n2(t) = n1(t+0.5)');
grid;
axis([-1 4 -.1 1.1])