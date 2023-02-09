%%  CT time history   CT_time.m
% defines a time history using the Von Karman spectra

close all
clear all
vrms = 25.908;            %  rms turbulence   m/s
U = 250;                 %  airspeed m/s
L = 762;                 %  characteristic scale wavelength  m
a = 1.339;

%%%%%%%%%%%%%%%%
fnyq = 10;   %   Nyquist Frequency  Hz
npts = 10000;  %  number of time points
nd2 = npts/2;
nd2p1 = nd2+1;  % number of freq points including zero Hz (DC)
dt = 1/2/fnyq;    %   sampling time
df = 1/npts/dt;   %   frequency increment
tt = dt*(1:npts);



%%%  PSDgust(w) = sigma_g2 2L/V  (1 + (8/3) (1.339 sigma L)^2  / (1 + (1.339sigma L)^2 )^(11/6)
PSDgust = zeros(size(0:df:fnyq));
PSDgustom = zeros(size(0:df:fnyq));
f = zeros(size(0:df:fnyq));
om = zeros(size(0:df:fnyq));
icount = 0;
vrms2 = vrms^2;
figure(1)
for f = 0:df:fnyq;
   icount = icount + 1;
   ff(icount) = f;
   numer = 1 + 8/3*(2*pi*a*f*L/U)^2;
   denom = (1 + (2*pi*a*f*L/U)^2)^(11/6);
   PSDgust(icount) = vrms2*2*L/U * numer / denom;
end
loglog(ff,PSDgust)
xlabel('freq Hz')
ylabel('PSD of gust velocity  (m/s)^2 / Hz')
figure
loglog(ff,PSDgust.*ff)
xlabel('freq Hz')
ylabel('PSD of gust velocity  (m/s)^2')


icount = 0;
figure
for f = 0:df:fnyq;
   icount = icount + 1;
   om(icount) = f/U;
   numer = 1 + 8/3*(a*L*f/U)^2;
   denom = (1 + (a*L*f/U)^2)^(11/6);
   PSDgustom(icount) = vrms2*L/pi * numer / denom;
end
loglog(om,PSDgustom)
xlabel('scaled freq (rad/m)')
ylabel('PSD of gust velocity')

%%%%%%%%%%%%%%%%%%%%  
%%%%%%   Generate equivalent time history with same amplitude
%%%%%%%%%%%%%%%%%%%%

amp = sqrt(PSDgust.*ff);      %  define amplitude   %%  complex values a + jb
a = zeros(size(amp));
b = zeros(size(amp));
a(1) = 0;    %%%  set dc value to zero
b(1) = 0;
for ii = 2:nd2
    a(ii) = amp(ii)*(-1 + 2*rand);
    ssign = 2 * (rand > 0.5) - 1;
    b(ii) = ssign * sqrt(amp(ii)^2 - a(ii)^2);
    ab(ii) = sqrt(a(ii)^2 + b(ii)^2);
    phas(ii) = atan2(b(ii),a(ii));
end
ssign = 2 * (rand > 0.5) - 1;
a(nd2p1) = ssign * amp(nd2p1);
    
    b(nd2p1) = 0;
    ab(nd2p1) = sqrt(a(nd2p1)^2 + b(nd2p1)^2);
    phas(nd2p1) = atan2(b(nd2p1),a(nd2p1));


figure
subplot(211)
amp(1) = 0;
plot(ff,ab,'r',ff,amp)
ylabel('Amplitude  m/s')
xlabel('Freq (Hz)')
subplot(212)
plot(ff,phas*180/pi)
ylabel('Phase')
xlabel('Degrees')
   
%% Flip and add the frequency values

aa = [a fliplr(a(2:nd2))];
bb = [b -fliplr(b(2:nd2))];

y = ifft(aa+1i*bb);

figure
plot(tt,y)
ylabel('gust vel m/s')
xlabel('time')

%%%  check on the time history
YF = fft(y);
YF = YF(1:nd2p1);
% figure
% subplot(211)
% plot(ff,abs(YF))
% subplot(212)
% plot(ff,angle(YF)*180/pi)

%%%%%  Loop to reduce the crest factor


for jj = 1:10
crest(jj) = (max(y) - min(y)) / rms(y);
clip = 0.7*max(abs(y));    
    
for ii=1:npts
    if y(ii) > clip
        y(ii) = clip;
    elseif y(ii) < -clip
        y(ii) = -clip;
    end
end

%%%   sliced time history
% figure
% plot(tt,y)
% ylabel('gust vel m/s')
% xlabel('time')

YF = fft(y);
YF = YF(1:nd2p1);
% figure
% subplot(211)
% plot(ff,abs(YF))
% subplot(212)
% plot(ff,angle(YF)*180/pi)

%%%% scale so that the amplitude is the same as the original - multiply by
%%%% the same amount for the real and imag parts so keep the phase
YF = YF./abs(YF).*amp;

% figure  
% subplot(211)
% plot(ff,abs(YF))
% subplot(212)
% plot(ff,angle(YF)*180/pi)
aa = real(YF);
bb = imag(YF);
aa = [aa fliplr(aa(2:nd2))];
bb = [bb -fliplr(bb(2:nd2))];

y = ifft(aa+1i*bb);
figure
plot(tt,y)
ylabel('gust vel m/s')
xlabel('time')
end

figure
plot(crest,'xr')
xlabel('iteration')
ylabel('crest factor')