close all
clear
% Partitioned block frequency domain adaptive filtering NLMS and 
fid=fopen('aecFar.pcm', 'rb'); % Load far end
rrin=fread(fid,inf,'int16');
fclose(fid); 
x =rrin;
fid=fopen('aecNear.pcm', 'rb'); % Load near end
ssin=fread(fid,inf,'int16');
fclose(fid);
d=ssin;

% x = audioread('C:\Users\garry\Desktop\ref.wav');
% d = audioread('C:\Users\garry\Desktop\ref_query_22_R.wav');
% x = x * 2^15;
% d = d * 2^15;

[Mx,Nx] = size(x);
[Md,Nd] = size(d);

%  Variable initialization
ntr = length(x);            %  temporary number of iterations 
M = 10;
N = 320;                   %  block length
TN = 2*N;

ntrB = floor(length(x)/N);  %  temporary number of block iterations
y = zeros(size(x));         %  initialize output signal vector
e = y;                      %  initialize error signal vector
FFTX = zeros(TN,M);         %  initialize temporary FFT input signal buffer
E = zeros(TN,1);            %  initialize temporary error signal buffer
nnM = 1:M-1;                %  index variable used for FFT input signal buffer
nnMp1 = nnM + 1;            %  index variable used for FFT input signal buffer
nnN = 1:N;                  %  index variable used for input signal buffer
nnNpN = N+1:TN;             %  index variable used for input signal buffer
nnf = ntrB*N+1:ntr;         %  index variable used for final block
Nf = length(nnf);           %  length of final block
FFTW = zeros(TN,M); %  initialize and assign frequency domain Coefficients
normFFTX = ones(TN,1);         %  initialize and assign FFT input signal powers
FFTX(:,nnM) = 0;  %  initialize and assign FFT input signal States
mu = 1/(2*N);            %  assign step size
bet = 0.999;%(1-1/(3*M*N))^N;          %  assign averaging factor
ombet = 1 - bet;            %  compute (1 - h.AvgFactor)
lam = 1;            %  assign Leakage
Offset = 0.001;          %  assign Offset
ZNM = zeros(N,M);           %  assign (N x M) zero matrix-dimensional zero vector
OMt = ones(1,M);            %  assign M-dimensional zero vector
fs =16000;             %the default frequecy

%  Initialize temporary input signal buffer

if isreal(x),
    X = real(ifft(FFTX(:,1)));
else
    X = ifft(FFTX(:,1));
end;
energy =0;
%  Main loop 

for n=1:ntrB,
    nn = ((n-1)*N+1):(n*N); %  index for current signal blocks
    FFTX(:,nnMp1) = FFTX(:,nnM);  %  shift temporary FFT input signal buffers over
    X(nnN) = X(nnNpN);      %  shift temporary input signal buffer up
    X(nnNpN) = x(nn);       %  assign current input signal block
    FFTX(:,1) = fft(X);     %  compute FFT of input signal vector
    Y = ifft(sum(FFTW.*FFTX,2));  %  compute current output signal vector
    y(nn) = Y(nnNpN);       %  assign current output signal block
    e(nn) = d(nn) - y(nn);  %  assign current error signal block
    E(nnNpN) = mu*e(nn);    %  assign current error signal vector
    normFFTX = bet*normFFTX + ombet*real(FFTX(:,1).*conj(FFTX(:,1))); 
    energy(n+1) = 0.9*energy(n) + 0.1* max(real(FFTX(:,1).*conj(FFTX(:,1))));
                    %  update FFT input signal powers
    FFTEinvnormFFTX = (fft(E)./(normFFTX + 1e4*N))*OMt;         
                            %  compute power-normalized FFT of error signal matrix 
                            
    G = ifft(FFTEinvnormFFTX.*conj(FFTX));
 
    
                            %  compute gradient matrix
    G(nnNpN,:) = ZNM;       %  impose gradient constraint
    FFTW = FFTW + fft(G);  %  update frequency domain coefficient matrix
    %FFTW(:,3:end)=0;
end
% fid=fopen('aecout.pcm', 'wb'); % Load far end
% fwrite(fid,e,'int16');
% fclose(fid);
audiowrite('a.wav',e,fs);

