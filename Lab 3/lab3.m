%%
%Combined Script File (Highlight a section to run it)   
%Problem A.3
%sinc = @(x) (x == 0) + (x ~= 0) .* (sin(x) ./ x); % Handle the case when x is 0
dn1 = @(n) 0.5*(abs(n)==3)+0.25*(abs(n)==1); %x1(t) has Dns of 1/2 for n=+-3 and 1/4 for n=+-1
dn2 = @(n) 0.5.*((sin(n*pi*0.5)./(n*pi*0.5))); %Sinc(x) = sin(x)/x Replace D0 manually
dn3 = @(n) 0.25.*((sin(n*pi*0.25)./(n*pi*0.25)));

%% 
%Problem A.4a 
%Plotting dn1
n = (-5:5); %Integer values for n from -5 to 5
figure;

%Plotting |Dn| (Magnitude of Dn)
subplot(1,2,1); %Row 1 Column 1 Plot
stem(n, abs(dn1(n))); %stem creates a straight line from x-axis, and abs(dn1(n)) is the magnitude of Dn so its the absolute value of Dn
xlabel('n');
ylabel('dn1');
grid on;

%Plotting ∠(Dn) (Angle of Dn)
subplot(1,2,2); %Row 1 Column 2 Plot
stem(n, angle(dn1(n))); %angle(dn1(n)) computes the phase angle of Dn
xlabel('n');
ylabel('∠dn1 (rad)');
grid on;

%% 
%Problem A.4a 
%Plotting dn2
n = (-5:5);
figure;
subplot(1,2,1);
stem(n, abs(dn2(n)));
xlabel('n');
ylabel('dn2');
grid on;

subplot(1,2,2);
stem(n, angle(dn2(n)));
xlabel('n');
ylabel('∠dn2 (rad)');
grid on;

%% 
%Problem A.4a 
%Plotting dn3
n = (-5:5);
figure;
subplot(1,2,1);
stem(n, abs(dn3(n)));
xlabel('n');
ylabel('dn3');
grid on;

subplot(1,2,2);
stem(n, angle(dn3(n)));
xlabel('n');
ylabel('∠dn3 (rad)');
grid on;

%% 
%Problem A.4b 
%Plotting dn1
n = (-20:20);
figure;

subplot(1,2,1);
stem(n, abs(dn1(n)));
xlabel('n');
ylabel('dn1');
grid on;

subplot(1,2,2);
stem(n, angle(dn1(n)));
xlabel('n');
ylabel('∠dn1 (rad)');
grid on;

%% 
%Problem A.4b 
%Plotting dn2
n = (-20:20);
figure;
subplot(1,2,1);
stem(n, abs(dn2(n)));
xlabel('n');
ylabel('dn2');
grid on;

subplot(1,2,2);
stem(n, angle(dn2(n)));
xlabel('n');
ylabel('∠dn2 (rad)');
grid on;

%% 
%Problem A.4b 
%Plotting dn3
n = (-20:20);
figure;
subplot(1,2,1);
stem(n, abs(dn3(n)));
xlabel('n');
ylabel('dn3');
grid on;

subplot(1,2,2);
stem(n, angle(dn3(n)));
xlabel('n');
ylabel('∠dn3 (rad)');
grid on;

%% 
%Problem A.4c 
%Plotting dn1
n = (-50:50);
figure;
subplot(1,2,1);
stem(n, abs(dn1(n)));
xlabel('n');
ylabel('dn1');
grid on;

subplot(1,2,2);
stem(n, angle(dn1(n)));
xlabel('n');
ylabel('∠dn1 (rad)');
grid on;

%% 
%Problem A.4c
%Plotting dn2
n = (-50:50);
figure;
subplot(1,2,1);
stem(n, abs(dn2(n)));
xlabel('n');
ylabel('dn2');
grid on;

subplot(1,2,2);
stem(n, angle(dn2(n)));
xlabel('n');
ylabel('∠dn2 (rad)');
grid on;

%% 
%Problem A.4c
%Plotting dn3
n = (-50:50);
figure;
subplot(1,2,1);
stem(n, abs(dn3(n)));
xlabel('n');
ylabel('dn3');
grid on;

subplot(1,2,2);
stem(n, angle(dn3(n)));
xlabel('n');
ylabel('∠dn3 (rad)');
grid on;

%% 
%Problem A.4d
%Plotting dn1
n = (-500:500);
figure;

subplot(1,2,1);
stem(n, abs(dn1(n)));
xlabel('n');
ylabel('dn1');
grid on;

subplot(1,2,2);
stem(n, angle(dn1(n)));
xlabel('n');
ylabel('∠dn1 (rad)');
grid on;

%% 
%Problem A.4d
%Plotting dn2
n = (-500:500);
figure;
subplot(1,2,1);
stem(n, abs(dn2(n)));
xlabel('n');
ylabel('dn2');
grid on;

subplot(1,2,2);
stem(n, angle(dn2(n)));
xlabel('n');
ylabel('∠dn2 (rad)');
grid on;

%% 
%Problem A.4d
%Plotting dn3
n = (-500:500);
figure;
subplot(1,2,1);
stem(n, abs(dn3(n)));
xlabel('n');
ylabel('dn3');
grid on;

subplot(1,2,2);
stem(n, angle(dn3(n)));
xlabel('n');
ylabel('∠dn3 (rad)');
grid on;

%% 
%Problem A.5
function A5(dn,a,To) % When synthesizing x(t) we need Dn, n-value and the fundamental frequency (wo) which can be calculated from the known To
t = (-300:300);
n = (-a:a); %vector of n-values
wo = 2*pi/To; %Fundamental Frequency
x = zeros(size(t)); %This is a vector that has the same size as t and all the values in it are currently set to zero but will be replaced with the signal values
for i = 1:length(n) %For loop used for recovering and storing value of x(t) into x for each value of n 
    x = x + dn(i).*exp(1j*n(i)*wo*t);
end
plot(t,x);
xlabel('t');
ylabel('x(t)');
grid on;
end

%% 
%Problem A.6 dn1, -5 <= n <= 5
A5(dn1,5,20);

%% 
%Problem A.6 dn1, -20 <= n <= 20
A5(dn1,20,20);

%% 
%Problem A.6 dn1, -50 <= n <= 50
A5(dn1,50,20);

%% 
%Problem A.6 dn1, -500 <= n <= 500
A5(dn1,500,20);

%% 
%Problem A.6 dn2, -5 <= n <= 5
A5(dn2,5,20);

%% 
%Problem A.6 dn2, -20 <= n <= 20
A5(dn2,20,20);

%% 
%Problem A.6 dn2, -50 <= n <= 50
A5(dn2,50,20);

%% 
%Problem A.6 dn2, -500 <= n <= 500
A5(dn2,500,20);

%% 
%Problem A.6 dn3, -5 <= n <= 5
A5(dn3,5,40);

%% 
%Problem A.6 dn3, -20 <= n <= 20
A5(dn3,20,40);

%% 
%Problem A.6 dn3, -50 <= n <= 50
A5(dn3,50,40);

%% 
%Problem A.6 dn3, -500 <= n <= 500
A5(dn3,500,40);




