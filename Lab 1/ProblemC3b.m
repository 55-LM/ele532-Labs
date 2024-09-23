t = (0:0.01:4)';
alpha = [1, 3, 5, 7];
T = t.*ones(1,4);
s_a = exp(-2)*exp(-T*diag(alpha)).*cos(4*pi*T).*u(T);
plot(t, s_a)
xlabel('t');
ylabel('s_α(t)');
legend({'α = 1', 'α = 3', 'α = 5', 'α = 7'});
grid;
