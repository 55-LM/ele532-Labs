u = @(t) 1.0*(t>=0);
h1 = @(t) exp(t/5).*u(t);
h2 = @(t) 4*exp(-t/5).*u(t);
h3 = @(t) 4*exp(-t).*u(t);
h4 = @(t) 4*(exp(-t/5) - exp(-t)).*u(t);

t = -1:0.001:5;

%Plotting h1(t)
figure;

subplot(4,1,1);
plot(t, h1(t));
title('h1(t) = exp(t/5) * u(t)');
xlabel('t');
ylabel('h1(t)');
grid on;

%Plotting h2(t)
subplot(4,1,2);
plot(t, h2(t));
title('h2(t) = 4 * exp(-t/5) * u(t)');
xlabel('t');
ylabel('h2(t)');
grid on;

%Plotting h3(t)
subplot(4,1,3);
plot(t, h3(t));
title('h3(t) = 4 * exp(-t) * u(t)');
xlabel('t');
ylabel('h3(t)');
grid on;

%Plotting h4(t)
subplot(4,1,4);
plot(t, h4(t));
title('h4(t) = 4 * (exp(-t/5) - exp(-t)) * u(t)');
xlabel('t');
ylabel('h4(t)');
grid on;