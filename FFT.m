function FFT(signal,fps)

Ts = 1/fps;                                             % Sampling Interval
Fs = 1/Ts;                                              % Sampling Frequency
Fn = Fs/2;                                              % Nyquist Frequency
L = size(signal,1);
NFFT = 2^nextpow2(L);                                   % For Efficiency
signal = signal - mean(signal);                        % Subtract Mean To Show Other Peaks
FTW = fft(signal, NFFT)/L;
Fv = linspace(0, 1, NFFT/2+1)*Fn;                       % Frequency Vector
Iv = 1:numel(Fv);                                       % Index Vector
figure
plot(Fv, abs(FTW(Iv,:))*2)                              % All Channels
grid
xlabel('Frequency (Hz)')
ylabel('Magnitude')