# Seiche-Detection
# A simple way of detecting coastal seiches will be presented in this document implementing Empirical Mode Decomposition, Hilbert Tranform, and comparing results with previous works in this area of research.
# Matlab Code
clc,clear,close
addpath('C:\Users\Astro E. Munoz\Documents\INGE5996\INGE5996 HW7');
expdata = load('S7.txt');
fs = 200;
dt = 1/fs;
%accelerations in cm/s^2
sensor1 = expdata(:,1);
sensor2 = expdata(:,2);
sensor3 = expdata(:,3);
np = length(sensor1);
t = 0:dt:(np-1)*dt;
bl = zeros(size(t));
%Maximum accelerations in three sensors
[maxampa1,maxposa1] = findmaxpos(t,sensor1);
[maxampa2,maxposa2] = findmaxpos(t,sensor2);
[maxampa3,maxposa3] = findmaxpos(t,sensor3);
figure;
subplot(3,1,1);plot(t,sensor1,'b-',maxposa1,maxampa1,'ro',t,bl,'r-');ylabel('accel.[cm/s^2]');title('Signals with Noise');
subplot(3,1,2);plot(t,sensor2,'b-',maxposa2,maxampa2,'ro',t,bl,'r-');ylabel('accel.[cm/s^2]');
subplot(3,1,3);plot(t,sensor3,'b-',maxposa3,maxampa3,'ro',t,bl,'r-');ylabel('accel.[cm/s^2]');xlabel('time[s]');

%Fourier Transform
[freqs1,mag1,phase1] = OSFS(t,sensor1.',0,0.001,'none');
[freqs2,mag2,phase2] = OSFS(t,sensor2.',0,0.001,'none');
[freqs3,mag3,phase3] = OSFS(t,sensor3.',0,0.001,'none');
figure;
subplot(211); plot(freqs1,mag1); title('FFT sensor 1 results'); ylabel('Magnitude');xlim([0,100]);
subplot(212); plot(freqs1,phase1); ylabel('Phase'); xlabel('Freq [hz]');xlim([0,100]);
figure;
subplot(211); plot(freqs2,mag2); title('FFT sensor 2 results'); ylabel('Magnitude');xlim([0,100]);
subplot(212); plot(freqs2,phase2); ylabel('Phase'); xlabel('Freq [hz]');xlim([0,100]);
figure;
subplot(211); plot(freqs3,mag3); title('FFT sensor 3 results'); ylabel('Magnitude');xlim([0,100]);
subplot(212); plot(freqs3,phase3); ylabel('Phase'); xlabel('Freq [hz]');xlim([0,100]);


%Filtered Fourier Coefficients
nF1 = fft(sensor1);
nF2 = fft(sensor2);
nF3 = fft(sensor3);
nm = 0:np-1;
fr = nm*fs/np;
figure;
subplot(3,1,1);plot(fr,abs(nF1));title('Fourier Coefficients');
subplot(3,1,2);plot(fr,abs(nF2));
subplot(3,1,3);plot(fr,abs(nF3));xlabel('freq.[Hz]');
%frequency used to filter 6Hz
nF1(122:np-122) = 0;
nF2(122:np-122) = 0;
nF3(122:np-122) = 0;
figure;
subplot(3,1,1);plot(fr,abs(nF1));title('Filtered Fourier Coefficients');
subplot(3,1,2);plot(fr,abs(nF2));
subplot(3,1,3);plot(fr,abs(nF3));xlabel('freq.[Hz]');
%Inverse Fourier Transform for filtering
nsensor1 = ifft(nF1);
nsensor2 = ifft(nF2);
nsensor3 = ifft(nF3);
[maxampan1,maxposan1] = findmaxpos(t,nsensor1);
[maxampan2,maxposan2] = findmaxpos(t,nsensor2);
[maxampan3,maxposan3] = findmaxpos(t,nsensor3);
figure;
subplot(3,1,1);plot(t,nsensor1,'b-',maxposan1,maxampan1,'ro',t,bl,'r-');ylabel('accel.[cm/s^2]');title('Filtered Signals in the Frequency Domain');
subplot(3,1,2);plot(t,nsensor2,'b-',maxposan2,maxampan2,'ro',t,bl,'r-');ylabel('accel.[cm/s^2]');
subplot(3,1,3);plot(t,nsensor3,'b-',maxposan3,maxampan3,'ro',t,bl,'r-');ylabel('accel.[cm/s^2]');xlabel('time[s]');


%Spectrogram and STFT
[S1,F1,T1,P1] = spectrogram(nsensor1,1000,990,2001,fs);
[S2,F2,T2,P2] = spectrogram(nsensor2,1000,990,2001,fs);
[S3,F3,T3,P3] = spectrogram(nsensor3,1000,990,2001,fs);
%Sensor 1
figure;
subplot(1,2,1);imagesc(T1,F1,10*log10(P1));title('spectrogram sensor #1');
xlabel('time[s]');ylabel('freq[Hz]');set(gca,'YDir','normal');ylim([0,10]);
subplot(1,2,2);imagesc(T1,F1,abs(S1));title('STFT sensor #1');
xlabel('time[s]');ylabel('freq[Hz]');set(gca,'YDir','normal');ylim([0,10]);
%Sensor 2
figure;
subplot(1,2,1);imagesc(T2,F2,10*log10(P2));title('spectrogram sensor #2');
xlabel('time[s]');ylabel('freq[Hz]');set(gca,'YDir','normal');ylim([0,10]);
subplot(1,2,2);imagesc(T2,F2,abs(S2));title('STFT sensor #2');
xlabel('time[s]');ylabel('freq[Hz]');set(gca,'YDir','normal');ylim([0,10]);
%Sensor 3
figure;
subplot(1,2,1);imagesc(T3,F3,10*log10(P3));title('spectrogram sensor #3');
xlabel('time[s]');ylabel('freq[Hz]');set(gca,'YDir','normal');ylim([0,10]);
subplot(1,2,2);imagesc(T3,F3,abs(S3));title('STFT sensor #3');
xlabel('time[s]');ylabel('freq[Hz]');set(gca,'YDir','normal');ylim([0,10]);

%HHT Implementation on Filtered Signal
nimfs    = 6; % number of imfs desired
maxit    = 100; % max number of iterations
stopcrit = 0.5*[0.05,0.5,0.05]; % rec- [0.05,0.5,0.05] [theta1 theta2 alpha]
displ    = 1;
[imf1,ORT1,NB_ITERATIONS1] = emd(nsensor1,'T',t,'STOP',stopcrit,'MAXITERATIONS',maxit,'MAXMODES',nimfs,'display',displ);
[imf2,ORT2,NB_ITERATIONS2] = emd(nsensor2,'T',t,'STOP',stopcrit,'MAXITERATIONS',maxit,'MAXMODES',nimfs,'display',displ);
[imf3,ORT3,NB_ITERATIONS3] = emd(nsensor3,'T',t,'STOP',stopcrit,'MAXITERATIONS',maxit,'MAXMODES',nimfs,'display',displ);
%Sensor #1
figure;
subplot(611);plot(t,imf1(1,:),'b'); title('IMF 1 Sensor # 1'); 
subplot(612);plot(t,imf1(2,:),'b'); title('IMF 2 Sensor # 1'); 
subplot(613);plot(t,imf1(3,:),'b'); title('IMF 3 Sensor # 1'); 
subplot(614);plot(t,imf1(4,:),'b'); title('IMF 4 Sensor # 1'); 
subplot(615);plot(t,imf1(5,:),'b'); title('IMF 5 Sensor # 1'); 
subplot(616);plot(t,imf1(6,:),'b'); title('IMF 6 Sensor # 1'); 

figure;
subplot(211);plot(t,imf1(7,:)); title('residue');%Residue of the Signal
subplot(212);plot(t,imf1(1,:)+imf1(2,:)+imf1(3,:)+imf1(4,:)+imf1(5,:)+imf1(6,:)+imf1(7,:)); title('reconstructed');%Reconstructed Signal
h1 = hilbert(imf1(1,:)); h2 = hilbert(imf1(2,:));h3 = hilbert(imf1(3,:));
if1 = angle(h1(2:end).*conj(h1(1:end-1)))/(2*pi)/dt;
if2 = angle(h2(2:end).*conj(h2(1:end-1)))/(2*pi)/dt;
if3 = angle(h3(2:end).*conj(h3(1:end-1)))/(2*pi)/dt;

%Instantaneous Frequencies Best Line Fits
nt = t(400:3600);
n = length(nt);
[yquadfit,a,b,c] = quadfit(nt,if1(400:3600),n);
[yfit1,m1,ifb1] = linearfit(nt,if2(400:3600));
[yfit2,m2,ifb2] = linearfit(nt,if3(400:3600));


nif1 = interp1(nt,yquadfit,t,'spline');
nif2 = interp1(nt,yfit1,t,'spline');
nif3 = interp1(nt,yfit2,t,'spline');

[nyquadfit,na1,nb1,nc1] = quadfit(t,nif1,np);

if1quadfiteq = na1 + nb1*t + nc1*t.^2;

figure;
plot(t(2:end),if1,t(2:end),if2,t(2:end),if3,t,nif1,'r-',t,nif2,'k-',t,nif3,'b-'); xlabel('time [s]'); ylabel('freq. [Hz]')

%Damped Signal Analysis
%Finding the instant amplitude and instant phase:
IAmp1 = abs(h1); % Instant amplitude
IPhase1 = angle(h1); % Instant phase
IAmp2 = abs(h2); % Instant amplitude
IPhase2 = angle(h2); % Instant phase
IAmp3 = abs(h3); % Instant amplitude
IPhase3 = angle(h3); % Instant phase

figure;
subplot(2,1,1);plot(t,IAmp1);title('IMF1 HT Analysis');ylabel('Instant Amplitude');
subplot(2,1,2);plot(t,IPhase1);ylabel('Instant Phase');xlabel('time[s]');
figure;
subplot(2,1,1);plot(t,IAmp2);title('IMF2 HT Analysis');ylabel('Instant Amplitude');
subplot(2,1,2);plot(t,IPhase2);ylabel('Instant Phase');xlabel('time[s]');
figure;
subplot(2,1,1);plot(t,IAmp3);title('IMF3 HT Analysis');ylabel('Instant Amplitude');
subplot(2,1,2);plot(t,IPhase3);ylabel('Instant Phase');xlabel('time[s]');

%Using the instant amplitude as the enveloping signal
figure;
subplot(3,1,1);plot(t,imf1(1,:),t,IAmp1,'r-');title('IMFs with Instantaneous Amplitudes');
subplot(3,1,2);plot(t,imf1(2,:),t,IAmp2,'r-');
subplot(3,1,3);plot(t,imf1(3,:),t,IAmp3,'r-');xlabel('time[s]');


ampn1 = length(t(540:3436));
[ampycubefit1,a1,b1,c1,d1] = cubefit(t(540:3436),IAmp1(540:3436),ampn1);
nIAmp1 = interp1(t(540:3436),ampycubefit1,t,'spline');


ampn2 = length(t(400:3600));
[ampycubefit2,a2,b2,c2,d2] = cubefit(t(400:3600),IAmp2(400:3600),ampn2);
nIAmp2 = interp1(t(400:3600),ampycubefit2,t,'spline');

ampn3 = length(t(600:4000));
[ampycubefit3,a3,b3,c3,d3] = cubefit(t(600:4000),IAmp3(600:4000),ampn3);
nIAmp3 = interp1(t(600:4000),ampycubefit3,t,'spline');

%Finding the amplitude & damping ratio
%IMF 1
ifint = cumtrapz(t,if1quadfiteq);
s1 = sin(2*pi*ifint).*nIAmp1;
figure;
plot(t,imf1(1,:),t,s1,'r');title('IMF1 with enveloping signal equation');xlabel('time[s]');

%IMF 2
IA2np = length(nIAmp2);
invIAmp2 = zeros(size(nIAmp2));
invIAmp2(1,np:-1:1) = nIAmp2(1,1:IA2np);

linamp2 = log(invIAmp2);
[ampyfit2,ampm2,ampb2] = linearfit(t,linamp2);
A2 = exp(ampb2);
damprat2 = sqrt(1/(2*pi*(mean(nif2)^2)/(ampm2^2)+1));
wn2 = -1*(ampm2/damprat2);
invmod2 = A2*exp(-1*damprat2*wn2*t);
mod2 = zeros(size(invmod2));
mod2(1,1:length(invmod2)) = invmod2(1,length(invmod2):-1:1);
s2 = mod2.*sin(2*pi*mean(nif2)*t);

figure;
plot(t,imf1(2,:),t,s2,'r');title('IMF2 with enveloping signal equation');xlabel('time[s]');

%IMF 3
linamp3 = log(nIAmp3);
[ampyfit3,ampm3,ampb3] = linearfit(t,linamp3);
A3 = exp(ampb3);
damprat3 = sqrt(1/(2*pi*(mean(nif3)^2)/(ampm3^2)+1));
wn3 = -1*(ampm3/damprat3);
mod3 = A3*exp(-1*damprat3*wn3*t);
s3 = mod3.*sin(2*pi*mean(nif3)*t);
figure;
plot(t,imf1(3,:),t,s3,'r');title('IMF3 with enveloping signal equation');xlabel('time[s]');

figure;
subplot(2,1,1);plot(t,nIAmp2,t,imf1(2,:),'r');title('Interpolated Amplitudes over IMF1 & IMF2')
subplot(2,1,2);plot(t,nIAmp3,t,imf1(3,:),'r');xlabel('time[s]');




%Sensor #2
figure;
subplot(611);plot(t,imf2(1,:),'b'); title('IMF 1 Sensor # 2'); 
subplot(612);plot(t,imf2(2,:),'b'); title('IMF 2 Sensor # 2'); 
subplot(613);plot(t,imf2(3,:),'b'); title('IMF 3 Sensor # 2'); 
subplot(614);plot(t,imf2(4,:),'b'); title('IMF 4 Sensor # 2'); 
subplot(615);plot(t,imf2(5,:),'b'); title('IMF 5 Sensor # 2'); 
subplot(616);plot(t,imf2(6,:),'b'); title('IMF 6 Sensor # 2'); 
figure;
subplot(211);plot(t,imf2(7,:)); title('residue');%Residue of the Signal
subplot(212);plot(t,imf2(1,:)+imf2(2,:)+imf2(3,:)+imf2(4,:)+imf2(5,:)+imf2(6,:)+imf2(7,:)); title('reconstructed');%Reconstructed Signal
%Sensor #3
figure;
subplot(611);plot(t,imf3(1,:),'b'); title('IMF 1 Sensor # 3'); 
subplot(612);plot(t,imf3(2,:),'b'); title('IMF 2 Sensor # 3'); 
subplot(613);plot(t,imf3(3,:),'b'); title('IMF 3 Sensor # 3'); 
subplot(614);plot(t,imf3(4,:),'b'); title('IMF 4 Sensor # 3'); 
subplot(615);plot(t,imf3(5,:),'b'); title('IMF 5 Sensor # 3'); 
subplot(616);plot(t,imf3(6,:),'b'); title('IMF 6 Sensor # 3'); 
figure;
subplot(211);plot(t,imf3(7,:)); title('residue');%Residue of the Signal
subplot(212);plot(t,imf3(1,:)+imf3(2,:)+imf3(3,:)+imf3(4,:)+imf3(5,:)+imf3(6,:)+imf3(7,:)); title('reconstructed');%Reconstructed Signal

%velocities in cm/s
vel1 = cumtrapz(sensor1)*dt;
vel2 = cumtrapz(sensor2)*dt;
vel3 = cumtrapz(sensor3)*dt;
%Maximum velocities in three sensors
[maxampv1,maxposv1] = findmaxpos(t,vel1);
[maxampv2,maxposv2] = findmaxpos(t,vel2);
[maxampv3,maxposv3] = findmaxpos(t,vel3);
figure;
subplot(3,1,1);plot(t,vel1,'b-',maxposv1,maxampv1,'ro',t,bl,'r-');ylabel('velocity[cm/s]');
subplot(3,1,2);plot(t,vel2,'b-',maxposv2,maxampv2,'ro',t,bl,'r-');ylabel('velocity[cm/s]');
subplot(3,1,3);plot(t,vel3,'b-',maxposv3,maxampv3,'ro',t,bl,'r-');ylabel('velocity[cm/s]');xlabel('time[s]');

%displacement in cm
disp1 = cumtrapz(vel1)*dt;
disp2 = cumtrapz(vel2)*dt;
disp3 = cumtrapz(vel3)*dt;
figure;
subplot(3,1,1);plot(t,disp1,'b-',t,bl,'r-');ylabel('disp.[cm]');
subplot(3,1,2);plot(t,disp1,'b-',t,bl,'r-');ylabel('disp.[cm]');
subplot(3,1,3);plot(t,disp1,'b-',t,bl,'r-');ylabel('disp.[cm]');xlabel('time[s]');

%Sum of all signals without noise contamination
S = s1 + s2 + s3;
figure;
subplot(2,1,1);plot(t,nsensor1);title('Original Signal without Noise')
subplot(2,1,2);plot(t,S);title('Created Signal without Noise');xlabel('time[s]');


%Noise Analysis
noiseF1 = fft(sensor1);
noiseF2 = fft(sensor2);
noiseF3 = fft(sensor3);
%Chirp frequencies 3Hz-6Hz
avef2 = mean(nif2);
avef3 = mean(nif3);
%Filter out the signals frequencies 
noiseF1(1:121) = 0;
noiseF1(np-121:np) = 0;
noiseF2(1:121) = 0;
noiseF2(np-121:np) = 0;
noiseF3(1:121) = 0;
noiseF3(np-121:np) = 0;

figure;
subplot(3,1,1);plot(fr,abs(noiseF1));title('Signal Filtered Fourier Coefficients');
subplot(3,1,2);plot(fr,abs(noiseF2));
subplot(3,1,3);plot(fr,abs(noiseF3));xlabel('freq.[Hz]');
%Inverse Fourier Transform for filtering the signal frequencies
noise1 = ifft(noiseF1);
noise2 = ifft(noiseF2);
noise3 = ifft(noiseF3);

figure;
subplot(3,1,1);plot(t,noise1);title('Noise in Instruments');ylim([-20,20]);
subplot(3,1,2);plot(t,noise2);ylim([-20,20]);
subplot(3,1,3);plot(t,noise3);xlabel('time[s]');ylim([-20,20]);

%SNR for 3 instruments
SNR1 = var(S)/var(noise1);
SNR2 = var(S)/var(noise2);
SNR3 = var(S)/var(noise3);

stdnoise1 = sqrt(var(noise1));
stdnoise2 = sqrt(var(noise2));
stdnoise3 = sqrt(var(noise3));

cnoise1 = stdnoise1*randn(size(t));
cnoise2 = stdnoise2*randn(size(t));
cnoise3 = stdnoise3*randn(size(t));

%Proposed Signals in Top of Original Signals with Noise

Inst1 = S + cnoise1;
Inst2 = S + cnoise2;
Inst3 = S + cnoise3;

figure;
plot(t,sensor1,t,Inst1,'r');title('Instrument # 1 Proposed and Original Signal');xlabel('time[s]');
figure;
plot(t,sensor2,t,Inst2,'r');title('Instrument # 2 Proposed and Original Signal');xlabel('time[s]');
figure;
plot(t,sensor3,t,Inst3,'r');title('Instrument # 3 Proposed and Original Signal');xlabel('time[s]');

%High Order Statistical Moments for Accelerations
%Original Signal
accmean1 = mean(sensor1);%mean
accmean2 = mean(sensor2);%mean
accmean3 = mean(sensor3);%mean
allmean = [accmean1,accmean2,accmean3];
accvar1 = var(sensor1);%variance
accvar2 = var(sensor2);%variance
accvar3 = var(sensor3);%variance
allvar = [accvar1,accvar2,accvar3];
accstd1 = sqrt(accvar1);%standard deviation
accstd2 = sqrt(accvar2);%standard deviation
accstd3 = sqrt(accvar3);%standard deviation
accskew1 = ((1/np)*sum((sensor1 - accmean1).^3))/accstd1^3;%skewness 
accskew2 = ((1/np)*sum((sensor2 - accmean2).^3))/accstd2^3;%skewness
accskew3 = ((1/np)*sum((sensor3 - accmean3).^3))/accstd3^3;%skewness
allskew = [accskew1,accskew2,accskew3];
acckurto1 = ((1/np)*sum((sensor1 - accmean1).^4))/accstd1^4;%kurtosis 
acckurto2 = ((1/np)*sum((sensor2 - accmean2).^4))/accstd2^4;%kurtosis 
acckurto3 = ((1/np)*sum((sensor3 - accmean3).^4))/accstd3^4;%kurtosis
allkurto = [acckurto1,acckurto2,acckurto3];
accRMS1 = sqrt(sum(sensor1.^2)/np);%RMS 
accRMS2 = sqrt(sum(sensor2.^2)/np);%RMS 
accRMS3 = sqrt(sum(sensor3.^2)/np);%RMS 
allRMS = [accRMS1,accRMS2,accRMS3];
accCF1 = max(abs(sensor1))/accRMS1;%Crest Factor
accCF2 = max(abs(sensor2))/accRMS2;%Crest Factor
accCF3 = max(abs(sensor3))/accRMS3;%Crest Factor
allCF = [accCF1,accCF2,accCF3];
accKF1 = max(abs(sensor1))*accRMS1;%K-Factor
accKF2 = max(abs(sensor2))*accRMS2;%K-Factor
accKF3 = max(abs(sensor3))*accRMS3;%K-Factor
allKF = [accKF1,accKF2,accKF3];

%Proposed Signal
accmeanp1 = mean(Inst1);%mean
accmeanp2 = mean(Inst2);%mean
accmeanp3 = mean(Inst3);%mean
allmeanp = [accmeanp1,accmeanp2,accmeanp3];
accvarp1 = var(Inst1);%variance
accvarp2 = var(Inst2);%variance
accvarp3 = var(Inst3);%variance
allvarp = [accvarp1,accvarp2,accvarp3];
accstdp1 = sqrt(accvarp1);%standard deviation
accstdp2 = sqrt(accvarp2);%standard deviation
accstdp3 = sqrt(accvarp3);%standard deviation
accskewp1 = ((1/np)*sum((Inst1 - accmeanp1).^3))/accstdp1^3;%skewness 
accskewp2 = ((1/np)*sum((Inst2 - accmeanp2).^3))/accstdp2^3;%skewness
accskewp3 = ((1/np)*sum((Inst3 - accmeanp3).^3))/accstdp3^3;%skewness
allskewp = [accskewp1,accskewp2,accskewp3];
acckurtop1 = ((1/np)*sum((Inst1 - accmeanp1).^4))/accstdp1^4;%kurtosis 
acckurtop2 = ((1/np)*sum((Inst2 - accmeanp2).^4))/accstdp2^4;%kurtosis 
acckurtop3 = ((1/np)*sum((Inst3 - accmeanp3).^4))/accstdp3^4;%kurtosis
allkurtop = [acckurtop1,acckurtop2,acckurtop3];
accRMSp1 = sqrt(sum(Inst1.^2)/np);%RMS 
accRMSp2 = sqrt(sum(Inst2.^2)/np);%RMS 
accRMSp3 = sqrt(sum(Inst3.^2)/np);%RMS 
allRMSp = [accRMSp1,accRMSp2,accRMSp3];
accCFp1 = max(abs(Inst1))/accRMSp1;%Crest Factor
accCFp2 = max(abs(Inst2))/accRMSp2;%Crest Factor
accCFp3 = max(abs(Inst3))/accRMSp3;%Crest Factor
allCFp = [accCFp1,accCFp2,accCFp3];
accKFp1 = max(abs(Inst1))*accRMSp1;%K-Factor
accKFp2 = max(abs(Inst2))*accRMSp2;%K-Factor
accKFp3 = max(abs(Inst3))*accRMSp3;%K-Factor
allKFp = [accKFp1,accKFp2,accKFp3];


figure;
subplot(2,1,1);bar(allmean);title('mean for the three instruments Original & Proposed');xlabel('instrument #');
subplot(2,1,2);bar(allmeanp);xlabel('instrument #');
figure;
subplot(2,1,1);bar(allvar);title('variance for the three instruments Original & Proposed');xlabel('instrument #');
subplot(2,1,2);bar(allvarp);xlabel('instrument #');
figure;
subplot(2,1,1);bar(allskew);title('skewness for the three instruments Original & Proposed');xlabel('instrument #');
subplot(2,1,2);bar(allskewp);xlabel('instrument #');
figure;
subplot(2,1,1);bar(allkurto);title('kurtosis for the three instruments Original & Proposed');xlabel('instrument #');
subplot(2,1,2);bar(allkurtop);xlabel('instrument #');
figure;
subplot(2,1,1);bar(allRMS);title('RMS for the three instruments Original & Proposed');xlabel('instrument #');
subplot(2,1,2);bar(allRMSp);xlabel('instrument #');
figure;
subplot(2,1,1);bar(allCF);title('crest factor for the three instruments Original & Proposed');xlabel('instrument #');
subplot(2,1,2);bar(allCFp);xlabel('instrument #');
figure;
subplot(2,1,1);bar(allKF);title('K-factor for the three instruments Original & Proposed');xlabel('instrument #');
subplot(2,1,2);bar(allKFp);xlabel('instrument #');



%Cross-Correlation Between the Proposed and Original Signals
%Instrument # 1
[cc1, lags1] = xcorr(Inst1,sensor1); 
tlags1 = dt*lags1;
cc1p = cc1(tlags1>=0);
figure;
subplot(211);plot(tlags1,cc1,'k');
title('Crosscorrelation between Original and Proposed Signal # 1')
subplot(212);plot(tlags1,cc1,':k'); xlabel('time [x]');
hold on;
plot(t,cc1p,'-b','linewidth',2);
hold off
%Instrument # 2
[cc2, lags2] = xcorr(Inst2,sensor2); 
tlags2 = dt*lags2;
cc2p = cc1(tlags2>=0);
figure;
subplot(211);plot(tlags2,cc2,'k');
title('Crosscorrelation between Original and Proposed Signal # 2')
subplot(212);plot(tlags2,cc2,':k'); xlabel('time [x]');
hold on;
plot(t,cc2p,'-b','linewidth',2);
hold off
%Instrument # 3
[cc3, lags3] = xcorr(Inst3,sensor3); 
tlags3 = dt*lags3;
cc3p = cc3(tlags3>=0);
figure;
subplot(211);plot(tlags3,cc3,'k');
title('Crosscorrelation between Original and Proposed Signal # 3')
subplot(212);plot(tlags3,cc3,':k'); xlabel('time [x]');
hold on;
plot(t,cc3p,'-b','linewidth',2);
hold off
