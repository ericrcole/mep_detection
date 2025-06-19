%design time-domain notch filter

tt = linspace(0,0.1,2200);
fs = 22000;

phi = 125; freq = 60; A = 5;
ytest = A*sin(freq*2*pi*tt+phi*(pi/180)) + randn(size(tt));

figure
plot(tt,ytest)

%% 
L = 2200;
f = fs/L*(0:L-1);

phase = angle(fft(ytest))* 180/pi+90;
pow = fft(ytest).*conj(fft(ytest))/L/L/2;

figure
subplot(3,1,1)
plot(f, pow)
hold on
xline(60, 'k--')
ylabel('Power')
xlim([0,100])
subplot(3,1,2)
plot(f,phase)
ylabel('Phase')
hold on
xline(60, 'k--')
xlabel('Frequency (Hz)')
xlim([0,100])

ortho_proj = @(A,B) (sum(A.*B)/(norm(B)^2))*B;

ind = find(f >=60, 1);
yest = ortho_proj(ytest, sin(60*2*pi*tt+phase(ind)*(pi/180)));

figure
plot(tt,ytest,'r')
hold on
plot(tt,yest,'k')
xlabel('Time (s)')

fprintf('True parameters: amp = %.2f, phase = %.2f\n\n', A, phi)
fprintf('Power at f=60: %.2f\n', max(yest))
fprintf('Phase at f=60: %.2f\n', phase(ind))
%%

tic
[yfilt, yest, A, phi] = td_notch2(ytest,tt,60,fs);
toc

figure
plot(ytest, 'r');
hold on
plot(yfilt, 'k')

%%

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

    ortho_proj = @(A,B) (sum(A.*B)/(norm(B)^2))*B;

    ind = find(f >=freq, 1);
    yest = ortho_proj(y, sin(60*2*pi*t+phase(ind)*(pi/180)));

    A = max(yest);
    phi = phase(ind);
    yfilt = y - yest;

end

function [yfilt, yest, A, mean_phase] = td_notch_mc(y,t,freq,fs)
    %y: signal for filtering
    %t: array of time indices (e.g. 0-100 ms)
    %freq: frequency for notch filt
    %fs: sample rate
    %
    %yfilt: signal with subtracted sinusoid
    %yest: fit sine component
    %A: amplitude of fit sinusoid
    %phi: phase of fit sinusoid 
    %
    %multichannel variant: estimates average phi from multi channels
    %before fitting and subtracting on individual channels
    
    if size(y,2) > size(y,1)
        y = transpose(y);
    end

    L = length(t);
    f = fs/L*(0:L-1);
    
    yfft = fft(y);
    
    phase = angle(yfft)* 180/pi+90;
    

    ortho_proj = @(A,B) (sum(A.*B)/(norm(B)^2))*B;

    ind = find(f >=freq, 1);
    mean_phase = mean(phase(ind,:)* 180/pi+90);

    y = transpose(y);
    yest = zeros(size(y));
    A = zeros(size(y,1));

    for kk = 1:size(y,1)
        yest(kk,:) = ortho_proj(y, sin(60*2*pi*t+mean_phase*(pi/180)));
        A(kk) = max(yest);
    end
    yfilt = y - yest;

end

function [yfilt, yest, A, phi] = td_notch2(y,t,freq,fs)
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

    ortho_proj = @(A,B) (sum(A.*B)/(norm(B)^2))*B;

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
