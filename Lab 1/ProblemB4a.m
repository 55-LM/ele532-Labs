n3 = @(t) n(t + 0.25);
plot(t, n3(t));
xlabel('t');
ylabel('n3(t) = n(t + 0.25)');
grid;
axis([-1 2.5 -.1 1.1])