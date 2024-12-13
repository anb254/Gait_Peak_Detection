function [HRatio] = StrideHarmRatio(acc,evenoverodd)

%The function has two inputs: an acceleration vector (acc) and binary
%indicator of whether the desired HR is in the vertical/AP direction
%(evenoverodd=1) or the ML direction (evenoverodd=0). The function uses FFT to create a
%frequency spectrum. The Haromic Ratio is calculated from the odd and even
%harmonics.

f=fft(acc); %Take fast Fourier Transform of acceleration signal
norm=abs(f/length(acc)); %Consider only positive frequencies and normalize by length of signal

%f is a vector of harmonic coefficients (0/constant, 1, 2, 3...). Even
%harmonics start at 3rd index (2nd harmonic), and odd harmonics start at
%2nd index (1st harmonic). We are only intersted in first 20 harmonics (not including constant),
%which capture frequencies up to 10Hz, which is sufficient for gait
%(Bellanca, J Biomechanics 2013).

even_harm=sum(norm(3:2:21)); %starting at 2nd harmonic (3rd index)
odd_harm=sum(norm(2:2:21)); %starting at 1st harmonic (2nd index)

if evenoverodd == 1; 
    HRatio= even_harm/odd_harm;     % Vertical and Anterioposterior
elseif evenoverodd == 0;
    HRatio = odd_harm/even_harm;    % Medio-lateral is dominated monophasic signals, so HR=odd/even
end
