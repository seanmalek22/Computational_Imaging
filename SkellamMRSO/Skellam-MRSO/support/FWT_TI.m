function [w,v]=FWT_TI(y,WType,WLevel,varargin)
% Poisson Rate Estimation Demostration Toolbox 
% version: 3.0, August 2008
% authors: Keigo Hirakawa & Patrick J. Wolfe, Harvard University
%
% Translation-Invariant Wavelet Transform
%    W=FWT_TI(X,...) performs forward wavelet transform. X can be 1D, 2D (image), 
%    or 3D (color image). Regardless of dimensionality, subbands are indexed by W(:,:,:,m). 
%    Y=IWT_TI(W,...) is an inverse of FWT_TI. Pad(X,sz) is used to pad the border, while 
%    the output is cropped with Pad(Y,-sz).
% 
%    [W,V]=FWT_TI(X,WType,WLevel,Modes)
%      'X'      - time/space domainsignal
%      'WType'  - type of wavelets (see wfilters)
%      'WLevel' - level of Haar wavelet decomposition
%      'Modes'  - optional tags (see below)
%      'W'      - translation-invariant wavelet coefficients (default)
%               - translation-invariant scaling coefficients (Mode='scaling')
%      'V'      - variance of valid part of the wavelet coefficients (default)
%               - mean of valid part of the scaling coefficients (Mode='scaling')
% 
%    [F]=IWT_TI(W,WType,WLevel,Modes)
%      'W'      - translation-invariant wavelet coefficients (default)
%      'WType'  - type of wavelets (see wfilters)
%      'WLevel' - level of Haar wavelet decomposition
%      'Modes'  - optional tags (see below)
%      'F'      - time/space domainsignal
%
%   Admissible 'Modes' are
%      'scaling' - W is a scaling coefficient rather than wavelets (Haar only)
%      'norm'    - 0,-1,+1 comprise transform (Haar only, still orthogonal but no longer orthonormal)
%      'valid'   - W contains valid parts of the transform 
%   'scaling' and 'valid' not available in IWT_TI.
%

if WLevel==0
    w=y;
    v=mean(y(:).^2);
    return
end

if size(y,2)==1 % 1D
    [Hd,Hs]=WT1_Filters(WType,WLevel); % generate filters
    Hd = permute(Hd,[1 3 2]); % subband index in 3rd dimension
    y=Pad(y,[size(Hd,1) 0]); % pad
else %2D/3D
    [Hd,Hs]=WT2_Filters(WType,WLevel); % generate filters
    y=Pad(y,size(Hd,1)); % pad
end

% optional modes tags
valid=0;
sc = 0;
for i=4:nargin
    switch varargin{i-4+1}
    case 'valid'
        valid=1;
    case 'scaling'
        sc = 1;
        Hd=abs(Hd);
    case 'square' % hidden option, for general wavelets.
        sc = 1;
        Hd=Hd.^2;
    case 'norm'
        for m=1:size(Hd,3);
            Hd(:,:,m)=Hd(:,:,m)/max(max(abs(Hd(:,:,m))));
        end
    end
end 

w=zeros([size(y,1),size(y,2),size(y,3),size(Hd,3)]);
for m=1:size(Hd,3);
    for c=1:size(y,3)
        w(:,:,c,m)=conv2(y(:,:,c),Hd(:,:,m),'same');
    end
end

s = size(Hd,1);


if sc == 0
    if size(y,2)==1
        v = mean(mean(Pad(w.^2,[-s,0]),2),1);
    else
        v = mean(mean(Pad(w.^2,-s)));
    end    
else
    if size(y,2)==1
        v = mean(mean(Pad(w,[-s,0]),2),1);
    else
        v = mean(mean(Pad(w,-s)));
    end
end

if valid==1
    if size(y,2)==1
        w=Pad(w,[-s,0]);
    else
        w=Pad(w,-s);
    end
end 

