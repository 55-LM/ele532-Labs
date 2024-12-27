%Combined Script File (Highlight a section to run it)

%L2_A1
% L2_A1.m : Chapter 2, MATLAB Program 1
% Script M-file determines characteristic roots of op-amp circuit.

% Set component values:
R = [1e4, 1e4, 1e4]; C = [1e-6, 1e-6];
% Determine the coefficients for characteristic equation:
A = [1, (1/R(1)+1/R(2)+1/R(3))/C(2), 1/(R(1)*R(2)*C(1)*C(2))];
% Determine characteristic roots:
lambda = roots(A);

%The poly() function returns the coeffcieints of a polynomial that has its
%roots specified.
poly(lambda);
%% 

%L2_A2
% Determine the coefficients for characteristic equation:
A = [1, (1/R(1)+1/R(2)+1/R(3))/C(2), 1/(R(1)*R(2)*C(1)*C(2))];
% Determine characteristic roots:
lambda = roots(A);

%The poly() function returns the coeffcieints of a polynomial that has its
%roots specified.
poly(lambda);

%Plotting the impulse response of the system for t=0 to t=0.1 with 0.0005
%increments
t = 0:0.0005:0.1;
h = @(t) -0.0045*exp(lambda(1).*t) + 0.0045*exp(lambda(2).*t);

plot(t,h(t))
title('Impulse Response of the System');
xlabel('Time (s)');
ylabel('h(t)');
grid on;
%% 

%L2_A3
function [lambda] = L2_A3(R,C)
% L2_A3.m : Chapter 2, MATLAB Program 3
% Function M-file finds characteristic roots of op-amp circuit.
% INPUTS:   R = length-3 vector of resistances
%           C = length-2 vector of capacitances% OUTPUTS:  
% lambda = characteristic roots
% Determine coefficients for characteristic equation:
A = [1, (1/R(1)+1/R(2)+1/R(3))/C(2), 1/(R(1)*R(2)*C(1)*C(2))];
% Determine characteristic roots:
lambda = roots(A);
end
%% 

% L2_B1.m : Chapter 2, MATLAB Program 4
% Script M-file graphically demonstrates the convolution process.

figure(1)  % Create figure window and make visible on screen
u = @(t) 1.0*(t>=0);
x = @(t) 1.5*sin(pi*t).*(u(t)-u(t-1));
h = @(t) 1.5*(u(t)-u(t-1.5))-u(t-2)+u(t-2.5);

dtau = 0.005; tau = -1:dtau:4;
ti = 0; tvec = -.25:.1:3.75;

y = NaN*zeros(1,length(tvec));  % Pre-allocate memory

for t = tvec
    ti = ti+1;  % Time index
    xh = x(t-tau).*h(tau); 
    lxh = length(xh);
    y(ti) = sum(xh.*dtau);  % Trapezoidal approximation of convolution integral
    
    subplot(2,1,1), plot(tau,h(tau),'k-',tau,x(t-tau),'k--',t,0,'ok');
    axis([tau(1) tau(end) -2.0 2.5]);
    
    patch([tau(1:end-1);tau(1:end-1);tau(2:end);tau(2:end)], ...
        [zeros(1,lxh-1);xh(1:end-1);xh(2:end);zeros(1,lxh-1)], ...
        [.8 .8 .8],'edgecolor','none');
    
    xlabel('\tau'); 
    title('h(\tau) [solid], x(t-\tau) [dashed], h(\tau)x(t-\tau) [gray]');
    
    c = get(gca,'children'); 
    set(gca,'children',[c(2);c(3);c(4);c(1)]);
    
    subplot(2,1,2), plot(tvec,y,'k',tvec(ti),y(ti),'ok');
    xlabel('t'); 
    ylabel('y(t) = \int h(\tau)x(t-\tau) d\tau');
    axis([tau(1) tau(end) -1.0 2.0]); 
    grid;
    drawnow;
end
%% 

% L2_B2.m : Chapter 2, MATLAB Program 5
% Script M-file graphically demonstrates the convolution process.

figure(1)  % Create figure window and make visible on screen
u = @(t) 1.0*(t>=0);
x = @(t) u(t)-u(t-2); %2.4-28(a) x(t)=u(t)-u(t-2)
h = @(t) (t+1).*(u(t+1)-u(t)); %2.4-30 h(t)=(t+1)(u(t+1)-u(t))

dtau = 0.005; tau = -2:dtau:4;
ti = 0; tvec = -1:.1:3;

y = NaN*zeros(1,length(tvec));  % Pre-allocate memory

for t = tvec
    ti = ti+1;  % Time index
    xh = x(t-tau).*h(tau); 
    lxh = length(xh);
    y(ti) = sum(xh.*dtau);  % Trapezoidal approximation of convolution integral
    
    subplot(2,1,1), plot(tau,h(tau),'k-',tau,x(t-tau),'k--',t,0,'ok');
    axis([tau(1) tau(end) -2.0 2.5]);
    
    patch([tau(1:end-1);tau(1:end-1);tau(2:end);tau(2:end)], ...
        [zeros(1,lxh-1);xh(1:end-1);xh(2:end);zeros(1,lxh-1)], ...
        [.8 .8 .8],'edgecolor','none');
    
    xlabel('\tau'); 
    title('h(\tau) [solid], x(t-\tau) [dashed], h(\tau)x(t-\tau) [gray]');
    
    c = get(gca,'children'); 
    set(gca,'children',[c(2);c(3);c(4);c(1)]);
    
    subplot(2,1,2), plot(tvec,y,'k',tvec(ti),y(ti),'ok');
    xlabel('t'); 
    ylabel('y(t) = \int h(\tau)x(t-\tau) d\tau');
    axis([tau(1) tau(end) -1.0 2.0]); 
    grid;
    drawnow;
end
%% 

% L2_B3a.m : Chapter 2, MATLAB Program 6
% Script M-file graphically demonstrates the convolution process.

figure(1)  % Create figure window and make visible on screen
u = @(t) 1.0*(t>=0);
A=1; 
B=2; 
x = @(t) A*(u(t-4)-u(t-6)); %2.4-27(a) x1(t)=A(u(t-4)-u(t-6))
h = @(t) B*(u(t+5)-u(t+4)); %2.4-27(a) x2(t)=B(u(t+5)-u(t+4))

dtau = 0.005; tau = -6:dtau:6;
ti = 0; tvec = -5:.1:5;

y = NaN*zeros(1,length(tvec));  % Pre-allocate memory

for t = tvec
    ti = ti+1;  % Time index
    xh = x(t-tau).*h(tau); 
    lxh = length(xh);
    y(ti) = sum(xh.*dtau);  % Trapezoidal approximation of convolution integral
    
    subplot(2,1,1), plot(tau,h(tau),'k-',tau,x(t-tau),'k--',t,0,'ok');
    axis([tau(1) tau(end) -2.0 2.5]);
    
    patch([tau(1:end-1);tau(1:end-1);tau(2:end);tau(2:end)], ...
        [zeros(1,lxh-1);xh(1:end-1);xh(2:end);zeros(1,lxh-1)], ...
        [.8 .8 .8],'edgecolor','none');
    
    xlabel('\tau'); 
    title('h(\tau) [solid], x(t-\tau) [dashed], h(\tau)x(t-\tau) [gray]');
    
    c = get(gca,'children'); 
    set(gca,'children',[c(2);c(3);c(4);c(1)]);
    
    subplot(2,1,2), plot(tvec,y,'k',tvec(ti),y(ti),'ok');
    xlabel('t'); 
    ylabel('y(t) = \int h(\tau)x(t-\tau) d\tau');
    axis([tau(1) tau(end) -1.0 2.0]); 
    grid;
    drawnow;
end
%% 

% L2_B3b.m : Chapter 2, MATLAB Program 7
% Script M-file graphically demonstrates the convolution process.

figure(1)  % Create figure window and make visible on screen
u = @(t) 1.0*(t>=0);
A=0.5; 
B=1.0; 
x = @(t) A*(u(t-3)-u(t-5)); %2.4-27(b) x1(t)=A(u(t-3)-u(t-5))
h = @(t) B*(u(t+5)-u(t+3)); %2.4-27(b) x2(t)=B(u(t+5)-u(t+3))

dtau = 0.005; tau = -6:dtau:3;
ti = 0; tvec = -6:.1:3;

y = NaN*zeros(1,length(tvec));  % Pre-allocate memory

for t = tvec
    ti = ti+1;  % Time index
    xh = x(t-tau).*h(tau); 
    lxh = length(xh);
    y(ti) = sum(xh.*dtau);  % Trapezoidal approximation of convolution integral
    
    subplot(2,1,1), plot(tau,h(tau),'k-',tau,x(t-tau),'k--',t,0,'ok');
    axis([tau(1) tau(end) -2.0 2.5]);
    
    patch([tau(1:end-1);tau(1:end-1);tau(2:end);tau(2:end)], ...
        [zeros(1,lxh-1);xh(1:end-1);xh(2:end);zeros(1,lxh-1)], ...
        [.8 .8 .8],'edgecolor','none');
    
    xlabel('\tau'); 
    title('h(\tau) [solid], x(t-\tau) [dashed], h(\tau)x(t-\tau) [gray]');
    
    c = get(gca,'children'); 
    set(gca,'children',[c(2);c(3);c(4);c(1)]);
    
    subplot(2,1,2), plot(tvec,y,'k',tvec(ti),y(ti),'ok');
    xlabel('t'); 
    ylabel('y(t) = \int h(\tau)x(t-\tau) d\tau');
    axis([tau(1) tau(end) -1.0 2.0]); 
    grid;
    drawnow;
end
%% 

% L2_B3h.m : Chapter 2, MATLAB Program 8
% Script M-file graphically demonstrates the convolution process.

figure(1)  % Create figure window and make visible on screen
u = @(t) 1.0*(t>=0);
x = @(t) exp(t).*(u(t+2)-u(t)); %2.4-27(h) x(t)=(e^t)(u(t+2)-u(t))
h = @(t) exp(-2*t).*(u(t)-u(t-1)); %2.4-27(h) h(t)=(e^-2t)(u(t)-u(t-1))

dtau = 0.005; tau = -4:dtau:6;
ti = 0; tvec = -5:.1:4;

y = NaN*zeros(1,length(tvec));  % Pre-allocate memory

for t = tvec
    ti = ti+1;  % Time index
    xh = x(t-tau).*h(tau); 
    lxh = length(xh);
    y(ti) = sum(xh.*dtau);  % Trapezoidal approximation of convolution integral
    
    subplot(2,1,1), plot(tau,h(tau),'k-',tau,x(t-tau),'k--',t,0,'ok');
    axis([tau(1) tau(end) -2.0 2.5]);
    
    patch([tau(1:end-1);tau(1:end-1);tau(2:end);tau(2:end)], ...
        [zeros(1,lxh-1);xh(1:end-1);xh(2:end);zeros(1,lxh-1)], ...
        [.8 .8 .8],'edgecolor','none');
    
    xlabel('\tau'); 
    title('h(\tau) [solid], x(t-\tau) [dashed], h(\tau)x(t-\tau) [gray]');
    
    c = get(gca,'children'); 
    set(gca,'children',[c(2);c(3);c(4);c(1)]);
    
    subplot(2,1,2), plot(tvec,y,'k',tvec(ti),y(ti),'ok');
    xlabel('t'); 
    ylabel('y(t) = \int h(\tau)x(t-\tau) d\tau');
    axis([tau(1) tau(end) -1.0 2.0]); 
    grid;
    drawnow;
end
%% 

%L2_C1
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
%% 

%L2_C3a
u = @(t) 1.0*(t>=0);
%In order to truncate the function we multiply the function by a unit step
%subtacted by a 20 unit delayed unit step
h = @(t) exp(t/5).*(u(t)-u(t-20)); %h1(t)
x = @(t) (u(t)-u(t-3)).*sin(5*t);

dtau = 0.005; tau = 0:dtau:20; %Modified
ti = 0; tvec = 0:.1:20;%Modified

y = NaN*zeros(1,length(tvec));  % Pre-allocate memory

for t = tvec
    ti = ti+1;  % Time index
    xh = x(t-tau).*h(tau); 
    lxh = length(xh);
    y(ti) = sum(xh.*dtau);  % Trapezoidal approximation of convolution integral
    
    subplot(2,1,1), plot(tau,h(tau),'k-',tau,x(t-tau),'k--',t,0,'ok');
    axis([tau(1) tau(end) -2.0 2.5]);
    
    patch([tau(1:end-1);tau(1:end-1);tau(2:end);tau(2:end)], ...
        [zeros(1,lxh-1);xh(1:end-1);xh(2:end);zeros(1,lxh-1)], ...
        [.8 .8 .8],'edgecolor','none');
    
    xlabel('\tau'); 
    title('h(\tau) [solid], x(t-\tau) [dashed], h(\tau)x(t-\tau) [gray]');
    
    c = get(gca,'children'); 
    set(gca,'children',[c(2);c(3);c(4);c(1)]);
    
    subplot(2,1,2), plot(tvec,y,'k',tvec(ti),y(ti),'ok');
    xlabel('t'); 
    ylabel('y(t) = \int h(\tau)x(t-\tau) d\tau');
    axis([tau(1) tau(end) -1.0 2.0]); 
    grid;
    drawnow;
end
%% 

%L2_C3b
u = @(t) 1.0*(t>=0);
%In order to truncate the function we multiply the function by a unit step
%subtacted by a 20 unit delayed unit step
h = @(t) 4*exp(-t/5).*(u(t)-u(t-20)); %h2(t)
x = @(t) (u(t)-u(t-3)).*sin(5*t);

dtau = 0.005; tau = 0:dtau:20; %Modified
ti = 0; tvec = 0:.1:20;%Modified

y = NaN*zeros(1,length(tvec));  % Pre-allocate memory

for t = tvec
    ti = ti+1;  % Time index
    xh = x(t-tau).*h(tau); 
    lxh = length(xh);
    y(ti) = sum(xh.*dtau);  % Trapezoidal approximation of convolution integral
    
    subplot(2,1,1), plot(tau,h(tau),'k-',tau,x(t-tau),'k--',t,0,'ok');
    axis([tau(1) tau(end) -2.0 2.5]);
    
    patch([tau(1:end-1);tau(1:end-1);tau(2:end);tau(2:end)], ...
        [zeros(1,lxh-1);xh(1:end-1);xh(2:end);zeros(1,lxh-1)], ...
        [.8 .8 .8],'edgecolor','none');
    
    xlabel('\tau'); 
    title('h(\tau) [solid], x(t-\tau) [dashed], h(\tau)x(t-\tau) [gray]');
    
    c = get(gca,'children'); 
    set(gca,'children',[c(2);c(3);c(4);c(1)]);
    
    subplot(2,1,2), plot(tvec,y,'k',tvec(ti),y(ti),'ok');
    xlabel('t'); 
    ylabel('y(t) = \int h(\tau)x(t-\tau) d\tau');
    axis([tau(1) tau(end) -1.0 2.0]); 
    grid;
    drawnow;
end
%% 

%L2_C3c
u = @(t) 1.0*(t>=0);
%In order to truncate the function we multiply the function by a unit step
%subtacted by a 20 unit delayed unit step
h = @(t) 4*exp(-t).*(u(t)-u(t-20)); %h3(t)
x = @(t) (u(t)-u(t-3)).*sin(5*t);

dtau = 0.005; tau = 0:dtau:20; %Modified
ti = 0; tvec = 0:.1:20;%Modified

y = NaN*zeros(1,length(tvec));  % Pre-allocate memory

for t = tvec
    ti = ti+1;  % Time index
    xh = x(t-tau).*h(tau); 
    lxh = length(xh);
    y(ti) = sum(xh.*dtau);  % Trapezoidal approximation of convolution integral
    
    subplot(2,1,1), plot(tau,h(tau),'k-',tau,x(t-tau),'k--',t,0,'ok');
    axis([tau(1) tau(end) -2.0 2.5]);
    
    patch([tau(1:end-1);tau(1:end-1);tau(2:end);tau(2:end)], ...
        [zeros(1,lxh-1);xh(1:end-1);xh(2:end);zeros(1,lxh-1)], ...
        [.8 .8 .8],'edgecolor','none');
    
    xlabel('\tau'); 
    title('h(\tau) [solid], x(t-\tau) [dashed], h(\tau)x(t-\tau) [gray]');
    
    c = get(gca,'children'); 
    set(gca,'children',[c(2);c(3);c(4);c(1)]);
    
    subplot(2,1,2), plot(tvec,y,'k',tvec(ti),y(ti),'ok');
    xlabel('t'); 
    ylabel('y(t) = \int h(\tau)x(t-\tau) d\tau');
    axis([tau(1) tau(end) -1.0 2.0]); 
    grid;
    drawnow;
end
%% 

%L2_C3d
u = @(t) 1.0*(t>=0);
%In order to truncate the function we multiply the function by a unit step
%subtacted by a 20 unit delayed unit step
h = @(t) 4*(exp(-t/5)-exp(-t)).*(u(t)-u(t-20)); %h4(t)
x = @(t) (u(t)-u(t-3)).*sin(5*t);

dtau = 0.005; tau = 0:dtau:20; %Modified
ti = 0; tvec = 0:.1:20;%Modified

y = NaN*zeros(1,length(tvec));  % Pre-allocate memory

for t = tvec
    ti = ti+1;  % Time index
    xh = x(t-tau).*h(tau); 
    lxh = length(xh);
    y(ti) = sum(xh.*dtau);  % Trapezoidal approximation of convolution integral
    
    subplot(2,1,1), plot(tau,h(tau),'k-',tau,x(t-tau),'k--',t,0,'ok');
    axis([tau(1) tau(end) -2.0 2.5]);
    
    patch([tau(1:end-1);tau(1:end-1);tau(2:end);tau(2:end)], ...
        [zeros(1,lxh-1);xh(1:end-1);xh(2:end);zeros(1,lxh-1)], ...
        [.8 .8 .8],'edgecolor','none');
    
    xlabel('\tau'); 
    title('h(\tau) [solid], x(t-\tau) [dashed], h(\tau)x(t-\tau) [gray]');
    
    c = get(gca,'children'); 
    set(gca,'children',[c(2);c(3);c(4);c(1)]);
    
    subplot(2,1,2), plot(tvec,y,'k',tvec(ti),y(ti),'ok');
    xlabel('t'); 
    ylabel('y(t) = \int h(\tau)x(t-\tau) d\tau');
    axis([tau(1) tau(end) -1.0 2.0]); 
    grid;
    drawnow;
end

