function [fH,varargout] = ske_mrso(G,wavLvl,varargin)
%SKE_MRSO: Minimum Risk Shrinkage Operator for Wavelet
%
% INPUT PARAMETERS:
%    G:         input image
%    wavLvl:    level of wavlet transform
%
% RETURN PRARMETERS:
%    F_hat:     image restoration
%
% Example:
%
%     FH = PD_MRSO(G, WAVLVL), return HMRSO denoising result with default
%     settings
% 
%     [FH,FB,FU] = PD_MRSO(G, WAVLVL), return results of HMRSO, BMRSO and
%     UMRSO with default settings
% 
%     PD_MRSO(G, WAVLVL, THRE, ALPHAS), allow customization settings,
%     for each wavlet subband, selection between HMRSO / UMRSO is
%     determined by threshold THRE. The smoothness of HMRSO is specified by
%     the vector ALPHAS.
%
%   Authors: 
%       Wu Cheng (ciivvi@gmail.com)
%       Keigo Hirakawa (khirakawa1@udayton.edu)
% 
%   Distribution Statement:
%       This software is for demonstration purpose only. Use of software 
%       for commercial purposes without a prior agreement with the authors 
%       is strictly prohibited. We do not guarantee the code's accuracy. 
%
%   Please acknowledge this work by citing the following paper:
%       Cheng, Wu; Hirakawa, Keigo (2015): Minimum Risk Wavelet 
%             Shrinkage Operator For Poisson Image Denoising. In: IEEE 
%             Transactions on Image Processing, 2015.

%   Copyright: 
%       University of Dayton, Intelligent Signal Systems Laboratory
%       campus.udayton.edu/~issl      

% default settings
thre = 80;
if mean2(G) < 6.4
    alphas = 0.9 * ones(3*wavLvl + 1);
else
    alphas = ones(3,1) * linspace(0.2,0.6,wavLvl);
    alphas = [1; alphas(:)];
end

% user settings
if nargin > 2
    thre = varargin{1};
end
if nargin > 3
    alphas = varargin{2};
    alphas = alphas(:);
end

% Hybrid Method of Optimal Skellam Shrink
fprintf('\nSkellam MRSO\n'); tic;

% wavelet transform
fprintf('   Wavelet transform\t\t');
[Y,~] = FWT_TI(G,'haar',wavLvl,'norm');
[T,~] = FWT_TI(G,'haar',wavLvl,'scaling','norm');
fprintf('... Done\n');

fprintf('   Denoising\t\t\t');
U = PD_MRSO_Uni_wav(Y,T,alphas);
B = PD_MRSO_Bi_wav(Y,T);

isU = mean(mean(T,1),2) > thre;
isU(1) = 0; % ignore DC channel
isU = find(isU == 1);

H = B;
H(:,:,:,isU) = U(:,:,:,isU);
fprintf('... Done\n');

fprintf('   Inverse wavelet transform\t');
fH = IWT_TI(H,'haar',wavLvl,'norm');
varargout{1} = IWT_TI(B,'haar',wavLvl,'norm');
varargout{2} = IWT_TI(U,'haar',wavLvl,'norm');
fprintf('... Done\n');
fprintf('Done! MRSO finished in %.0f sec.\n\n', toc);