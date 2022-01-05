% Parametric Signal Modeling lets us model a function with a smaller number
% of coefficients to describe the signal than whats necessary for the
% original signal. the lower the number of coefficients the less accurate
% the model is but the more space is saved for data transmission
N=512;
omega = linspace(0,2*pi,N);
sig = [0.5 1 2 1 0.5 0.25 0.125];
Rsig = autoCorrMatrix(sig,2);
[pCoeff, pG] = PSMCoe(sig,2);
pspect = PSMSpect(sig,3,N);

% here the first image shows the actual spectrum of a signal. the
% second image shows our model of that signal that's built with 4 coefficients
figure
subplot(2,1,1)
stem(omega,abs(fft(sig,N)))
subplot(2,1,2)
stem(omega,abs(pspect));


% These functions are for parametric signal modeling using inverse modeling
% Also known as inderect least squares modeling
% looking for online mathematical reference different than my notes.
% makes a P+1 by P+1 autocorrelation matrix 
function R = autoCorrMatrix(x,P)
% the entries of this matrix are just x convolved with it's flipped version
    N = length(x);
    if(P<N)
        v = conv(x,flip(x));
        R = zeros(P+1);
        for i=1:P+1
            R(i,:) = v(N+1-i:N+P+1-i);
        end
        return
    else
        R = 0;
        return
    end
end

% uses the autocorrelation matrix to solve for the model coefficients
% a1,a2,a3 and the gain G
function [a, G] = PSMCoe(x,P)
    R = autoCorrMatrix(x,P);
    M = R(2:P+1,2:P+1);
    z = -1*R(2:end,1);
    a = M\z;
    G = sqrt(R(1,2:end)*a + R(1,1));
end

% models the signal x with a P order parametric signal model and outputs
% it's spectrum from 0 to 2pi
% This model has the form 

% ______________________________1_________________________________________
% (1 + sum of (P coefficients) * (P complexexponentials))

% or

% ______________________________1_________________________________________
% (1 + sum_of_k_from_1_to_P_of{ (a[k]) * exp(-j*freqDomain*k) ) }
function s = PSMSpect(x,P,N)
    d = 0;
    [a,G] = PSMCoe(x,P);
    n = (0:N-1);
    for k=1:P
       d = d+a(k)*exp(1j*pi*2*n*k/N);
    end
    s = G./(1+d);
    return
end
