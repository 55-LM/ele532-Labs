r = @(t) t.*p(t);
plot(t,r(t));
xlabel('t');
ylabel('r(t) = t*p(t)');
grid;
axis([-1 2 -.1 1.1])