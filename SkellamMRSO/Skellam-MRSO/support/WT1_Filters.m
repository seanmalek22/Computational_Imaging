function [Hd,Hs] = WT1_Filters(Htype,M)

% generate filter coefficients for the wavelet transform

[L0,H0,L1,H1] = wfilters(Htype);  % Wavelet filters
% Lo_D, the decomposition low-pass filter
% Hi_D, the decomposition high-pass filter
% Lo_R, the reconstruction low-pass filter
% Hi_R, the reconstruction high-pass filter

% hd = sum(2.^[0:M-1])*(length(L0))+1; % size of the largest filter
hd = (2^M-1) * (length(L0))+1; % size of the largest filter

Hd = zeros(hd,M+1);

% initialize

LL = [L0 0]';
ll = [L0 0]';
hh = [H0 0]';
Hd((hd-size(hh,1))/2+1:end-(hd-size(hh,1))/2,2)=hh;
k  = 3;
for m=1:M-1
    ll = conv2(up(ll),LL); 
    hh = conv2(up(hh),LL); 
    Hd((hd-size(hh,1))/2+1:end-(hd-size(hh,1))/2,k)=hh;
    k = k+1;
end
Hd(:,1)=ll;

hd = sum(2.^[0:M-1])*(length(L1))+1; % size of the largest filter
Hs = zeros(hd,M+1);
% initialize
LL = [0 L1]';
ll = [0 L1]';
hh = [0 H1]';
Hs((hd-size(hh,1))/2+1:end-(hd-size(hh,1))/2,2)=hh;
k  = 3;
for m=1:M-1
    ll = conv2(up(ll),LL); 
    hh = conv2(up(hh),LL); 
    Hs((hd-size(hh,1))/2+1:end-(hd-size(hh,1))/2,k)=hh;
    k = k+1;
end
Hs(:,1)=ll;


% compensate for overcomplete wavelet redundancies
Hs(:,1)=Hs(:,1)/(2^(M));
for m=1:M
    Hs(:,m+1)=Hs(:,m+1)/(2^(m));
end

% reorder to go from coarse to fine
Hd(:,2:end)=Hd(:,end:-1:2);
Hs(:,2:end)=Hs(:,end:-1:2);
    
function H = up(h)
H=zeros(size(h,1)*2-1,1);
H(1:2:end)=h;

