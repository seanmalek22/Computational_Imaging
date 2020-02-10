function F=IWT_TI(w,WType,WLevel,varargin)
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
    F=w; 
    return;
end


if size(w,2)==1
    [Hd,Hs]=WT1_Filters(WType,WLevel);
    Hd = permute(Hd,[1 3 2]);
    Hs = permute(Hs,[1 3 2]);
else
    [Hd,Hs]=WT2_Filters(WType,WLevel);
end

valid=0;
for i=4:nargin
    switch varargin{i-4+1}
    case 'norm'
        for m=1:size(Hd,3);
            Hs(:,:,m)=Hs(:,:,m)*max(max(abs(Hd(:,:,m))));
        end
    end
end 

F=zeros([size(w,1),size(w,2),size(w,3)]);
for m=1:size(Hs,3)
for c=1:size(w,3)
    F(:,:,c)=F(:,:,c)+conv2(w(:,:,c,m),Hs(:,:,m),'same');
end
end

if size(w,2)==1
    F=Pad(F,[-size(Hd,1),0]);
else
    F=Pad(F,-size(Hd,1));
end

