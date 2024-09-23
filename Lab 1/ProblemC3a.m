t = (0:0.01:4);
alpha = [1, 3, 5, 7];
for a = alpha
    s_a = @(t) exp(-2)*exp(-a*t).*cos(4*pi*t).*u(t);
    plot(t, s_a(t));
    hold on;
end
xlabel('t');
ylabel('s_α(t)');
legend({'α = 1', 'α = 3', 'α = 5', 'α = 7'});
grid;
hold off;