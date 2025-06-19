function [yfilt, yest, A, phi] = td_notch(y,t,freq,fs)
    %time-domain notch filter: determines phase of component at frequency
    %'freq' and subtracts a sine wave in time domain
    %
    %y: signal for filtering
    %t: array of time indices (e.g. 0-100 ms)
    %freq: frequency for notch filt
    %fs: sample rate
    %
    %yfilt: signal with subtracted sinusoid
    %yest: fit sine component
    %A: amplitude of fit sinusoid
    %phi: phase of fit sinusoid

    L = length(t);
    f = fs/L*(0:L-1);

    yfft = fft(y);
    phase = angle(yfft)* 180/pi+90;

    %ortho_proj = @(A,B) (sum(A.*B)/(norm(B)^2))*B;

    ind = find(f >=freq, 1);

    %yest = ortho_proj(y, sin(60*2*pi*t+phase(ind)*(pi/180)));
    yest = sin(60*2*pi*t+phase(ind)*(pi/180));
    yscale = abs(median(y./yest));
    yest = yest*yscale;

    %A = max(yest);
    A=yscale;

    phi = phase(ind);
    yfilt = y - yest;

end
