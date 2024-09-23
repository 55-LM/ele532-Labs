n4 = @(t) n3(0.5*t);
plot(t, n4(t));
xlabel('t');
ylabel('n4(t) = n3(0.5*t)');
grid;
axis([-1 4 -.1 1.1])