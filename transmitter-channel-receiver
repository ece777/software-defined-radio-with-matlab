clear all
clf
 
% specification of impairments
 
cng=input('channel noise gain: try 0, 0.6 or 2 :: '); 
cdi=input('channel multipath: 0 for none, 1 for mild or 2 for harsh ::  ');
fo =input('tranmsitter mixer freq offset in %: try 0 or 0.01 ::  ');  
po =input('tranmsitter mixer phase offset in rad: try 0, 0.7 or 0.9 ::  ');
toper=input('baud timing offset as % of symb period: try 0, 20 or 30 ::  ');
so=input('symbol period offset: try 0 or 1 ::  ');
 
%TRANSMITTER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% encode text string as T-spaced PAM (+/-1, +/-3) sequence
str2='01234 I wish I were an Oscar Mayer wiener 56789 ';
str3 = 'A0Oh well whatever Nevermindl ';
str4 = 'Awell ';
n=text2bin(str3);          % change text into 7 bit binary using text2bin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% encode 
coded_txt = blockcode52_encode(n);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% biner to 4-PAM
j=1; mpam=zeros(1,ceil(length(coded_txt)/2)); 
for i=1:2:length(coded_txt)-1                 % and then into 4-PAM
  if coded_txt(i:i+1)==[0,0], mpam(j)=-3; end 
  if coded_txt(i:i+1)==[0,1], mpam(j)=-1; end
  if coded_txt(i:i+1)==[1,0], mpam(j)=1; end
  if coded_txt(i:i+1)==[1,1], mpam(j)=3; end
  j=j+1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=length(mpam);                                % 4-level signal of length N
% zero pad T-spaced symbol sequence to create upsampled T/M-spaced
%      sequence of scaled T-spaced pulses (with T = 1 time unit)
M=200-so; mup=zeros(1,N*M); mup(1:M:end)=mpam; % oversampling factor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SRRC pulse filter with T/M-spaced impulse response
L=0.5; p=(14.0)*srrc(L,0.3,M,0.4);             % blip pulse of width M
x=filter(p,1,mup);                  % convolve pulse shape with data
figure(1), plotspec(x,1/M)          % baseband signal spectrum
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% am modulation
t=1/M:1/M:length(x)/M;              % T/M-spaced time vector
fc=20;                              % carrier frequency
c=cos(2*pi*(fc*(1+0.01*fo))*t+po);  % carrier with offsets relative to rec osc
r=c.*x;                             % modulate message with carrier
figure(2),plot(r(1:10*M));
figure(3), plotspec(r,1/M);
figure(4), plotspec(p,1/M);
% OK!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CHANNEL
 
%IMPAIRMENT : FADING
ds=pow(r);                      % desired average power of signal
lr=length(r);                   % length of transmitted signal vector
fp=[ones(1,floor(0.2*lr)),0.5*ones(1,lr-floor(0.2*lr))]; % flat fading profile
r=r.*fp;                        % apply profile to transmitted signal vector
 
if cdi < 0.5,                  % channel ISI
  mc=[1 0 0];                  % distortion-free channel
elseif cdi<1.5,
  mc=[1 zeros(1,M) 0.28 zeros(1,2.3*M) 0.11]; % mild multipath channel  
else
  mc=[1 zeros(1,M) 0.28 zeros(1,1.8*M) 0.44]; % harsh multipath channel  
end
mc=mc/(sqrt(mc*mc'));          % normalize channel power
dv=filter(mc,1,r);             % filter transmitted signal through channel
nv=dv+cng*(randn(size(dv)));   % add Gaussian channel noise
to=floor(0.01*toper*M);        % fractional period delay in  sampler
rnv=nv(1+to:end);              % delay in on-symbol designation
rt=(1+to)/M:1/M:length(nv)/M;  % modified time vector with delayed message start 
rM=M+so;                       % receiver sampler timing offset (delay)
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RECEIVER
 
% automatic gain control (AGC)
g=zeros(1,lr); g(1)=1;          % initialize gain
nr=zeros(1,lr);
mu=0.0003;                      % stepsize
for i=1:lr-1
 nr(i)=g(i)*r(i);               % AGC output
 g(i+1)=g(i)-mu*(nr(i)^2-ds);   % adapt gain
end
rnv=nr;                           % received signal is still called r
 
% pllconverge.m simulate costas loop
% input rsc from pulrecsig.m
rpll=rnv;                                    % rsc is from pulrecsig.m 
fl=100; ff=[0 .01 .02 1]; fa=[1 1 0 0];   
h=remez(fl,ff,fa);                        % LPF design
mu=.003;                                  % algorithm stepsize
fc=20;                                  % assumed freq. at receiver
theta=zeros(1,length(t)); theta(1)=0;     % initialize estimate vector
zs=zeros(1,fl+1); zc=zeros(1,fl+1);       % initialize buffers for LPFs
for k=1:length(t)-1                       % z's contain past fl+1 inputs
  zs=[zs(2:fl+1), 2*rpll(k)*sin(2*pi*fc*t(k)+theta(k))];
  zc=[zc(2:fl+1), 2*rpll(k)*cos(2*pi*fc*t(k)+theta(k))];
  lpfs=fliplr(h)*zs'; lpfc=fliplr(h)*zc'; % new output of filters
  theta(k+1)=theta(k)-mu*lpfs*lpfc;       % algorithm update
end
 
figure(5),plot(t,theta),
title('Phase Tracking via the Costas Loop')
xlabel('time'); ylabel('phase offset')
 
theta;
phoff = theta;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%am demodulation of received signal sequence r
ram = rnv;
c2=cos(2*pi*fc*t + phoff);                   % synchronized cosine for mixing
x2=ram.*c2;                            % demod received signal
% LPF design
fl=100;                               % LPF length
fbe=[0 0.1 0.2 1]; damps=[1 1 0 0 ]; % design of LPF parameters
b=remez(fl,fbe,damps);               % create LPF impulse response
x3=2*filter(b,1,x2);                 % LPF and scale downconverted signal
figure(6),plotspec(x3,1/M)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % matchfilt.m: test of SNR maximization
%x3 = r; fl = 100;
recfilt=(15)*srrc(L,0.3,M,0.0);                    % receive filter H sub R
%recfilt=recfilt/sqrt(sum(recfilt.^2));  % normalize the pulse shape  
v=1/180*filter(fliplr(recfilt),1,x3);          % matched filter with data
figure(7), ul=floor((length(v)-124)/(4*rM)); 
plot(reshape(v(125:ul*4*rM+124),4*rM,ul)) % plot the eye diagram
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clock recovery algorithm
xcl = v; n = N; l = L; 
tnow=l*M+1; tau=0; xs=zeros(1,n);          % initialize variables
tausave=zeros(1,n); tausave(1)=tau; i=0;
mu=0.003;                                  % algorithm stepsize
delta=0.1;                                 % time for derivative
while tnow<length(xcl)-2*l*M                 % run iteration
  i=i+1;
  xs(i)=interpsinc(xcl,tnow+tau,l);          % interpolated value at tnow+tau
  x_deltap=interpsinc(xcl,tnow+tau+delta,l); % get value to the right
  x_deltam=interpsinc(xcl,tnow+tau-delta,l); % get value to the left
  dx=x_deltap-x_deltam;                    % calculate numerical derivative  
  qx=quantalph(xs(i),[-3,-1,1,3]);         % quantize xs to nearest 4-PAM symbol
  tau=tau+mu*dx*(qx-xs(i));                % alg update: DD 
  tnow=tnow+M; tausave(i)=tau;             % save for plotting
end
figure(8), subplot(2,1,1), plot(xs(1:i-2),'b.')    % plot constellation diagram
title('constellation diagram');
ylabel('estimated symbol values')
subplot(2,1,2), plot(tausave(1:i-2))               % plot trajectory of tau
ylabel('timing offset estimates'), xlabel('iterations')
tausave;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% downsample to symbol rate 
z=v(0.5*fl+rM:rM+so:end);  % set delay to first symbol-sample and increment by M
figure(9), plot([1:length(z)],z,'.')      % soft decisions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LSequalizer.m find a LS equalizer f for the channel b
n=3;      re = n;                        % length of equalizer - 1
delta=3;                          % use delay <=n*length(b) 
p=length(re)-delta;
RE=toeplitz(re(n+1:p),re(n+1:-1:1)); % build matrix R 
SE=re(n+1-delta:p-delta)';          % and vector SE
f=inv(RE'*RE)*RE'*SE                  % calculate equalizer f
Jmin=SE'*SE-SE'*RE*inv(RE'*RE)*RE'*SE     % Jmin for this f and delta
ye=filter(f,1,re);                  % equalizer is a filter
dec=sign(ye);                      % quantize and find errors
err=0.5*sum(abs(dec(delta+1:end)-re(1:end-delta)))
z = ye;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% decision device and symbol matching performance assessment
mprime=quantalph(z,[-3,-1,1,3])';         % quantize to +/-1 and +/-3 alphabet
cluster_variance=(mprime-z)*(mprime-z)'/length(mprime),     % cluster variance
 
j=1; y=zeros(1,2*length(mprime));  
rq=quantalph(mprime,[-3,-1,1,3]);             % quantize back to 4-PAM
for i=1:length(mprime)                        % translate back into 0-1 binary
  if rq(i)==3, y(j:j+1)=[1,1]; end 
  if rq(i)==1, y(j:j+1)=[1,0]; end 
  if rq(i)==-1, y(j:j+1)=[0,1]; end 
  if rq(i)==-3, y(j:j+1)=[0,0]; end 
  j=j+2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% decode 
y = blockcode52_decode(y);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% binary to text
ytext=bin2text(y)
numerr=length(find(ytext~=text));
 
lmp=length(mprime);   % number of recovered symbol estimates
percentage_symbol_errors=100*sum(abs(sign(mprime-mpam(1:lmp))))/lmp % symb err
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
