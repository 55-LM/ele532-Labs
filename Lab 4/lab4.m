%%
%Problem A.1 (Computed manually)

% Define the time vector and the signal x(t)
%N = 100; PulseWidth = 10;
%t = [0:1:(N-1)];
%x = [ones(1,PulseWidth), zeros(1,N-PulseWidth)];

%Convolution of x(t) and x(t) produces z(t)
%z = conv(x, x, 'same');  

%Time vector for z(t)
%t_z = [0:1:(length(z)-1)];

%Plot x(t) 
%figure;
%subplot(211); 
%stairs(t, x); grid on; 
%title('x(t)');
%xlabel('t');
%ylabel('x(t)');

%Plot z(t)
%subplot(212); 
%plot(t_z, z); grid on;
%title('z(t) = x(t)*x(t)');
%xlabel('t');
%ylabel('z(t)');

%%   
%Combined Script File (Select a section to run it) 
%Problem A.2

% Define the time vector and the signal x(t)
N = 100; PulseWidth = 10;
t = [0:1:(N-1)]; %Time vector from 0 to N-1 corresponding to sample size of 100
x = [ones(1,PulseWidth), zeros(1,N-PulseWidth)];

%Computing fourier transform of x(t) (X(ω)) using fft function (fast fourier transform)
Xw = fft(x);
f = [-(N/2):1:((N/2)-1)]*(1/N); %Frequency vector from -0.5 ((-(N/2))/N)) to 0.49 ((1-N/2)/N))

%Computing Z(ω)=X(ω)*X(ω)
Zw = Xw.*Xw;

%% 
%Problem A.3

%Plot magnitude of Z(ω) (|Z(ω)|)
figure;
subplot(211); 
plot(f, fftshift(abs(Zw))); grid on;
title('Magnitude of Z(ω)');
xlabel('Frequency (rad/s)');
ylabel('|Z(ω)|');

%Plot phase of Z(ω) (∠Z(ω))
subplot(212); 
plot(f, fftshift(angle(Zw))); grid on;
title('Phase of Z(ω)');
xlabel('Frequency (rad/s)');
ylabel('Phase (rads)');

%% 
%Problem A.4
%Convolution in Time == Product in Frequency
%We calculated the product of X(w) and X(w) (FT of x(t)) in A.2 which yielded Z(w)
%z(t), the convolution of x(t) and x(t) should yield the same plot as in A.1
%This proves that computing the product of two signals in frequency-domain
%equals the convolution of two signals in the time-domain

%Plot z(t) in time-domain (compute z(t) by using inverse fourier transform on Z(ω))
z = ifft(Zw);
figure;
subplot(211); 
plot(t, z); grid on;
title('Time-Domain of z(t)');
xlabel('t');
ylabel('z(t)');

zf = fft(z);
%Plot z(t) in frequency-domain (compute z(t) in frequency-domain using fourier transform)
subplot(212); 
plot(f, fftshift(abs(zf))); grid on;
title('Frequency-Domain of z(t)');
xlabel('Frequency (rad/s)');
ylabel('Phase (degree)');


%% 
%Problem A.5
%The frequency plots of X(w) are scaled with respect to the pulse width

%Pulse width changed to 5 or 25
PulseWidth = 5;

%Defining x(t)
x = [ones(1,PulseWidth), zeros(1,N-PulseWidth)];

%Computing fourier transform of x(t) (X(ω)) using fft function (fast fourier transform)
Xw = fft(x);

%Plot magnitude of X(ω) (|X(ω)|)
figure;
subplot(211); 
plot(f, fftshift(abs(Xw))); grid on;
title('Magnitude of X(ω) with Pulse Width = 5');
xlabel('Frequency (rads/s)');
ylabel('|X(ω)|');

%Plot phase of X(ω) (∠X(ω))
subplot(212); 
plot(f, fftshift(angle(Xw))); grid on;
title('Phase of X(ω) with Pulse Width = 5');
xlabel('Frequency (rad/s)');
ylabel('Phase (rads)');

%Pulse width changed to 5 or 25
PulseWidth25 = 25;

%Defining x(t)
x1 = [ones(1,PulseWidth25), zeros(1,N-PulseWidth25)];

%Computing fourier transform of x(t) (X(ω)) using fft function (fast fourier transform)
Xw1 = fft(x1);

%Plot magnitude of X(ω) (|X(ω)|)
figure;
subplot(211); 
plot(f, fftshift(abs(Xw1))); grid on;
title('Magnitude of X(ω) with Pulse Width = 25');
xlabel('Frequency (rads/s)');
ylabel('|X(ω)|');

%Plot phase of X(ω) (∠X(ω))
subplot(212); 
plot(f, fftshift(angle(Xw1))); grid on;
title('Phase of X(ω) with Pulse Width = 25');
xlabel('Frequency (rad/s)');
ylabel('Phase (rads)');

%% 
%Problem A.6 
%The first two plots essentially have the angular frequency varied which
%shows frequency shifted plots of the same signal
%The third plot has a cos component which means that there are two angular
%frequencies, one positive and one negative

PulseWidth = 10;
x = [ones(1,PulseWidth), zeros(1,N-PulseWidth)];
t = [0:1:(N-1)];

%Defining w+(t)
w_plus = x.*exp(1j*(pi/3)*t); 

%Computing fourier transform of w+(t)
W_plus = fft(w_plus);

%Plot magnitude of W+(ω) (|W+(ω)|)
figure;
subplot(211); 
plot(f, fftshift(abs(W_plus))); grid on;
title('Magnitude of W_+(ω)');
xlabel('Frequency (rads/s)');
ylabel('|W_+(ω)|');

%Plot phase of W+(ω) (∠W+(ω))
subplot(212); 
plot(f, fftshift(angle(W_plus))); grid on;
title('Phase of W_+(ω)');
xlabel('Frequency (rad/s)');
ylabel('Phase (rads)');

%Defining w-(t)
w_minus = x.*exp(-1j*(pi/3)*t); 

%Computing fourier transform of w-(t)
W_minus = fft(w_minus);

%Plot magnitude of W-(ω) (|W-(ω)|)
figure;
subplot(211); 
plot(f, fftshift(abs(W_minus))); grid on;
title('Magnitude of W_-(ω)');
xlabel('Frequency (rads/s)');
ylabel('|W_-(ω)|');

%Plot phase of W-(ω) (∠W-(ω))
subplot(212); 
plot(f, fftshift(angle(W_minus))); grid on;
title('Phase of W_-(ω)');
xlabel('Frequency (rad/s)');
ylabel('Phase (rads)');

%Defining wc(t)
w_c = x.*cos((pi/3)*t); 

%Computing fourier transform of wc(t)
W_c = fft(w_c);

%Plot magnitude of Wc(ω) (|Wc(ω)|)
figure;
subplot(211); 
plot(f, fftshift(abs(W_c))); grid on;
title('Magnitude of W_c(ω)');
xlabel('Frequency (rads/s)');
ylabel('|W_c(ω)|');

%Plot phase of Wc(ω) (∠Wc(ω))
subplot(212); 
plot(f, fftshift(angle(W_c))); grid on;
title('Phase of W_c(ω)');
xlabel('Frequency (rad/s)');
ylabel('Phase (rads)');

%% 
%Problem B.1

%Plot xspeech 
figure;
MagSpect(xspeech);
title('Magnitude Spectrum of Speech Signal (xspeech)');

%Plot hLPF2000
figure;
MagSpect(hLPF2000);
title('Magnitude Spectrum of Lowpass Filter with Passband 2.0 kHz (hLPF2000)');

%Plot hLPF2500
figure;
MagSpect(hLPF2500);
title('Magnitude Spectrum of Lowpass Filter with Passband 2.5 kHz (hLPF2500)');

%Plot hChannel
figure;
MagSpect(hChannel);
title('Magnitude Spectrum of Bandpass Communications Channel (hChannel)');

%% 
%Problem B.1 Cont.
N = length(xspeech); % Number of samples

%Filter First

%MODULATION
%Modulate the signal, xspeech, by transmitting it through the channel, hChannel

%Frequency shift of hChannel (bandpass communication channel)
F_0 = 6000; 
%Modulation of signal producing a symmetrical signal
x_osc = osc(F_0, N); 

%Modulated signal is multplied by original signal so that it is able to pass through bandpass channel
x_encode = xspeech.*x_osc; 

%Plot modulated signal of xspeech
figure;
MagSpect(x_encode);
title('Magnitude Spectrum of Encoded Signal');

%Transmit the signal, xspeech, over the channel, hChannel
x_channel = conv(x_encode, hChannel, 'same'); 

%Plot transmitted signal of xspeech
figure;
MagSpect(x_channel);
title('Magnitude Spectrum of Transmitted Signal');

%DEMODULATION
%Demodulate the signal, xspeech, that was transmitted through the channel using the lowpass filters

%Signal is decoded using the same modulated signal producing symmetrical signal 
x_decode = x_channel.*x_osc; 

%Plot demodulated xspeech signal
figure;
MagSpect(x_decode);
title('Magnitude Spectrum of Decoded Signal');

%Use lowpass filter to recover the signal, xspeech (use hLPF2000 or hLPF2500)
x_recover = conv(x_decode, hLPF2500, 'same'); %The demodulate signal has a region where the signal is recovered, so convoluting it with a unit impulse where it is active in the same frequency range will recover the signal

%Plot the recovered signal of xspeech
figure;
MagSpect(x_recover);
title('Magnitude Spectrum of Recovered Signal');

