function [Y] = Pad(X,n)
% Poisson Rate Estimation Demostration Toolbox 
% version: 3.0, August 2008
% authors: Keigo Hirakawa & Patrick J. Wolfe, Harvard University
%
% Reflexive Padding/Unpadding
%    Pad(X,sz) will pad the borders of X with reflexive padding.
%    The resulting signal will increase in size by sz in all directions.
% 
%    Y = Pad(X,sz)
%      'Y'     - padded signal
%      'X'     - input signal
%      'sz'    - increase by sz in each direction
%
%    If sz is negative, X will be cropped by sz in each direction.
%

% set parameters
if length(n)==1
    if ndims(X)==1
        n=[n 0];  % 1D expansion
    else
        n=[n,n];  % 2D expansion
    end
end

% padding
if n(1)>0 % if expanding...
    
    m=[size(X,1),size(X,2)]-1;

    % if the signal is too small for reflection, we use recursion
    if n(1)>m(1) % signal too small, begin recursion
        Y = Pad(X,[m(1),n(2)]);
        Y = Pad(Y,[n(1)-m(1),n(2)]);
    elseif n(2)>m(2) % signal too small, begin recursion
        Y = Pad(X,[n(1),m(2)]);
        Y = Pad(Y,[n(1),n(2)-m(2)]);        
    else % signal big enough. Do reflexive padding
        Y = [X(n(1)+1:-1:2,n(2)+1:-1:2,:,:)       X(n(1)+1:-1:2,:,:,:)        X(n(1)+1:-1:2,end-1:-1:end-n(2),:,:)
            X(:,n(2)+1:-1:2,:,:)                  X                           X(:,end-1:-1:end-n(2),:,:)
            X(end-1:-1:end-n(1),n(2)+1:-1:2,:,:)  X(end-1:-1:end-n(1),:,:,:)  X(end-1:-1:end-n(1),end-1:-1:end-n(2),:,:)];
    end
else % if cropping...
    n=-n;
    Y = X(n(1)+1:end-n(1),n(2)+1:end-n(2),:,:); %take the center area
end