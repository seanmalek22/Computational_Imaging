function X_hat = PD_MRSO_Bi_wav(Y,T)
% Bivariate Minimum Risk Shrinkage Operator
% Version 1.0, 2014
%
% INPUT PARAMETERS:
%    Y:         wavelet coefficient of image
%    T:         scaling coefficient of image
%
% RETURN PRARMETERS:
%    F_hat:     denoiseed wavelet coefficient
%
% Author:     Wu Cheng, Keigo Hirakawa
%             {chengw1, khirakawa1}@udayton.edu
%             Intelligent Signal Processing Laboratory
%             University of Dayton

% initialize
X_hat = Y;              % estimated wavelet coefficients
for c=1:size(Y,3)       % for each color
    for m=2:size(Y,4)   % for each subband
        
        y = Y(:,:,c,m);
        t = T(:,:,c,m);
        
        % skellam => poisson
        g0 = (t+y)/2; % g+
        g1 = (t-y)/2; % g-
        mG = max([g0(:); g1(:)]);
        
        % 2D histogram
        p = accumarray([g0(:),g1(:)]+1, 1, [mG+1 mG+1]);
        p = (p+p')/2; % impose symmetry: p(y)=p(-y)
        p = p/sum(p(:));
        
        % PY0: p(Y = y, T = T)      <=>     p(g0,g1)
        % PY1: p(Y = y+1, T = T+1)  <=>     p(g0+1,g1)
        % PY2: p(Y = y-1, T = T+1)  <=>     p(g0,g1+1)
        
        PY0 = p(g0+g1*size(p,1)+1);
        PY1 = p(min(g0+1,mG)+g1*size(p,1)+1);
        PY2 = p(g0+min(g1+1,mG)*size(p,1)+1);
        
        x_hat = (y.*(PY1+PY2)+(t+2).*(PY1-PY2))./(2*PY0);
        x_hat(sign(x_hat)~=sign(y))=0;
        x_hat = sign(y).*min(abs(x_hat),abs(y));
        
        X_hat(:,:,c,m) = x_hat;
    end
end

end