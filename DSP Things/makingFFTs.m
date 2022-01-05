N = 512;
t = linspace(0,1,N);
omega = linspace(0,2*pi,N);
oscFreq = 4;
sig = square(oscFreq*t*2*pi);

figure
subplot(5,1,1)
stem(t,sig)

subplot(5,1,2)
stem(omega,abs(fft(sig,N)))

subplot(5,1,3)
stem(omega,abs(dft(sig,N)))

subplot(5,1,4)
stem(omega,abs(myFFT(sig,N)))
xfft = myFFT(sig,N);

subplot(5,1,5)
xresynth = iFFT(xfft,N);
stem(t,real(xresynth))

% this covers convolution using the FFT made here and compares it to
% matlab's conv function
figure
reverb = exp(-1/8*(1:N)).*rand(1,N);
sigVerb = myConv(sig,reverb);

subplot(2,1,1)
stem(linspace(0,2,length(sigVerb)),real(sigVerb));
subplot(2,1,2)
stem(linspace(0,2,2*N),[real(conv(sig,reverb)) 0])



% Simple implementation of the DFT
function F = dft(X,N)
    F = zeros(1,N);
    for k = 0:N-1
        for n = 0:N-1
            F(k+1) = F(k+1) + X(n+1)*exp(-1j*2*pi*k*n/N);
        end
    end    
    return
end

function F = myFFT(x,N)
    F = zeros(1,N);
    if(N<=1)
        F = x(1);
        return
    else
%       Recursively split into odd and even arrays and take the N/2 point
%       FFT of them
        fEven = myFFT(x(1:2:N-1),N/2);
        fOdd  = myFFT(x(2:2:N),N/2);
%       Odd array needs a corrective phase factor
        k = (0:N/2-1);
        phase = exp(-2*1j*pi*k/N);
%       add the even and odd index FFTs together.
        F(1:N/2) = fEven + phase.*fOdd;
%       the subtraction here is a result of the fact that our k index only
%       went half the length of N. If it were the length of N we could just
%       use even+odd*phase but the vector lengths in code wouldnt match and
%       result in more lines of code
        F(N/2+1:N) = fEven - phase.*fOdd;
        return
    end
end

% the inverse FFT utilizes the fact that two applications of a fourier
% transform results in a flipped signal in time and a scaling factor. so we
% just switch them back

function x = iFFT(F,N)
   x = 1/N*myFFT(F,N);
   x = flip(x);
   return
end

% Convolution is just mutliplication in frequency domain
function y = myConv(x1,x2)
% zeropad the two signals so that we get the 'tail' of the convolution in
% our output
    x1pad = [x1 zeros(1,length(x1))];
    x2pad = [x2 zeros(1,length(x2))];
    L1 = length(x1pad);
    L2 = length(x2pad);
    X1 = myFFT(x1pad,L1);
    X2 = myFFT(x2pad,L2);
    Y =  X1.*X2;
    y = iFFT(Y,length(Y));
end