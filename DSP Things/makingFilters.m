% simple implementation of a band pass filter using poles in a z transform
% Q w and G can be our knobs to tune the filter response where Q adjusts Q
% value, w adjusts the center frequency and G is the gain of the filter
N = 512;
omega = linspace(0,2*pi,N);
Q = 0.8;
w = pi*6/8;
G = 1;
impulse = zeros(1,N);
impulse(1) = 1;
y = zeros(1,N*2);
for n=1:N
    m = n+N;
    y(m) = G*impulse(n)+2*Q*cos(w)*y(m-1)-Q^2*y(m-2);
end 
y = y(N+1:end);

figure
subplot(2,1,1)
stem(omega,abs(fft(y,N)));


% simple implementation of notch/bandrejection filter using zeros a z transform
Q = 1;
w = pi*3/8;
G = 1;
impulse = zeros(1,2*N);
impulse(N+1) = 1;
y = zeros(1,N*2);
for n=1:N
    m = n+N;
    y(m) = G*(impulse(m)-2*Q*cos(w)*impulse(m-1)+Q^2*impulse(m-2));
end 
y = y(N+1:end);
subplot(2,1,2)
stem(omega,abs(fft(y,N)));



% Biquad filter design:
% this uses the design method of billinear transforms to create digital
% filters that approximate the transfer functions of analog filters without
% aliasing due to digital sampling https://en.wikipedia.org/wiki/Bilinear_transform

% Biquad simply means that we're dealing with at most a time delay of 2
% samples in either feed forward or feedback configuration.
% https://en.wikipedia.org/wiki/Digital_biquad_filter

% designing using billinear tranforms involves the use of the formulas 
% s = (1-z^-1)/(1+z^-1) and Wc = 2/T*tan(Wd/2) where Wc is some continuous
% frequency and Wd is the corresponding digital frequency between -pi and pi.
% T is the sampling period

% This involves a lot of algebra in order to find the right feedback
% coefficients and we'd need to do it in a way where we can quickly change
% a few values in order to change the behavior of the filter. 
% specifically Gain, Q and cutoff/center frequency
% to make it easier on ourselves lets refer to this great resource
% https://webaudio.github.io/Audio-EQ-Cookbook/audio-eq-cookbook.html

impulse = zeros(1,N);
impulse(1) = 1;

y = LPF(impulse,pi/3,.8,1);

figure
stem(omega,abs(fft(y,N)));


% Low Pass Filter
function y = LPF(x,w,Q,G)
    N = length(x);
    y = zeros(1,N);
    x = [zeros(1,N) x ];
    alpha = sin(w)/(2*Q);
    a0 = 1 + alpha; b0 = (1-cos(w))/2;
    a1 = -2*cos(w); b1 = (1-cos(w));
    a2 = 1 - alpha; b2 = (1-cos(w))/2;
    for n=1:N
        m = n+N;
        y(m) = G * 1/a0 * (b0*x(m) + b1*x(m-1) + b2*x(m-2) - a1*y(m-1) - a2*y(m-2));
    end 
    x = x(N+1:end);
    y = y(N+1:end);
    return
end