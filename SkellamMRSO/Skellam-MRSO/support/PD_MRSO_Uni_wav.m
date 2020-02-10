function X_hat = PD_MRSO_Uni_wav(Y,T,alphas)
% Univariate Minimum Risk Shrinkage Operator
% Version 1.0, 2014
%
% INPUT PARAMETERS:
%    Y:         wavelet coefficient of image
%    T:         scaling coefficient of image
%    alphas:    Irwin mixture parameter
%
% RETURN PRARMETERS:
%    F_hat:     denoiseed wavelet coefficient
%
% Author:     Wu Cheng, Keigo Hirakawa
%             {chengw1, khirakawa1}@udayton.edu
%             Intelligent Signal Processing Laboratory
%             University of Dayton

% wavelet transform
X_hat = Y;

% initialization
ET_Y_cell = cell(size(Y,3),size(Y,4)); % E[T|Y]
mY_cell = cell(size(Y,3),size(Y,4));
PY_cell = cell(size(Y,3),size(Y,4));

%% Estimate E[T|Y=y] and P[Y=y]
% Merge horizontal and vertical subbands
for c = 1:size(Y,3)
    for m = 2:size(Y,4)
        
        y = Y(:,:,c,m);
        t = T(:,:,c,m);
        y = y(:);
        t = t(:);
        mY = max(abs(y));
        sz = 2*mY+1;
        
        % add [y,t] pair (-mY,0) and (+mY,0) to force range of T_Y be symmetric
        y = [y;-mY;mY];
        t = [t;0;0];
        
        y = y + mY +1;
        PY = max(accumarray(y,1,[sz,1]),1);
        ET_Y = accumarray(y,t,[sz,1]) ./ PY;
        
        % force symmetric
        ET_Y = make_symmetric(ET_Y);
        PY = make_symmetric(PY);
        
        % smoothness
        win_size = 7;
        lpf = [1:win_size,win_size-1:-1:1]';
        lpf = lpf/sum(lpf);
        
        ET_Y= conv(ET_Y,lpf,'same');
        PY= conv(PY,[1:4 3:-1:1],'same');
        
        PY = PY/sum(PY);
        PY_cell{c,m} = PY;
        ET_Y_cell{c,m} = ET_Y;
        mY_cell{c,m} = mY;
        
    end
end

%% Denoising
for c=1:size(Y,3) % for each color
    for m=2:size(Y,4) % for each subband
        
        y = Y(:,:,c,m);
        t = T(:,:,c,m);
        y = y(:);
        t = t(:);
        ET_Y = ET_Y_cell{c,m};
        mY = mY_cell{c,m};
        PY = PY_cell{c,m};
        
        % Irwin fix
        PY_Irwin= pY_integral(y,t);
        a = alphas(m);
        PY = a*PY + (1-a)*PY_Irwin;
        
        % P[Y=y]
        idx = y + mY + 1;
        PY0 = PY(idx);
        
        % P[Y=y-1]
        idx = y + mY;
        idx = max(idx,1);
        PYm = PY(idx);
        
        % E[T|Y=y-1]
        ET_Ym = ET_Y(idx);
        
        % P[Y=y+1]
        idx = y + mY + 2;
        idx = min(idx,2*mY+1);
        PYp = PY(idx);
        
        % E[T|Y=y1]
        ET_Yp = ET_Y(idx);
        
        x_hat = ( (y+1+ET_Yp).*PYp + (y-1-ET_Ym).*PYm ) ./ (2*PY0 + eps);
        X_hat(:,:,c,m) = reshape(x_hat,size(Y(:,:,c,m)));
    end
end
end

%% --------------------------------------------
%         make pdf symmetric
%  --------------------------------------------
function y = make_symmetric(x)
% force x to be symmetric
% x is assumed to be properly padded
x = x(:);
x1 = flipud(x);

y = x+x1;
idx = x&x1;
y(idx) = (x(idx) + x1(idx)) / 2;
end

%% --------------------------------------------
%         Skellam Distribution
%  --------------------------------------------
function p = Distribution_Skellam_x0(y,s)
% u =-s + y/2.*(log(max(s+x,eps))-log(max(s-x,eps)))+log(besseli(y,sqrt(s.^2-x.^2),1));
% p = exp(-s)*besseli(y,s);
p = besseli(y,s,1);
end

%% --------------------------------------------
%         Skellam Distribution
%  --------------------------------------------
function py = pY_integral(y,s)
% p[Y=y|X=0] = integral { p[Y=y|X=0,S=s] * p[S=s] } ds
% where distribution of y is Skellam: p[Y|X,S] = Skellam(y; (s+x)/2, (s-x)/2))
y = y(:);
s = s(:);
% add [y,t] pair (-mY,0) and (+mY,0) to force range of T_Y be symmetric
mY = max(abs(y));
maxS = max(s);
minY = min(y);
maxY = max(y);

% val range: 0~mS, ind range: 1:mS+1
ps = accumarray(s+1,1,[maxS+1,1]);
ps = ps/sum(ps);

sval = find(ps~=0);
ps = ps(sval);
sval = sval - 1;

y_unique = unique(y);

% Implementation 1: for loop
py_unique = zeros(size(y_unique));
for i = 1:numel(y_unique)
    py_i = Distribution_Skellam_x0(y_unique(i),sval) .* ps;
    py_unique(i) = mean(py_i);
end
py_unique = mean(py_unique,2);

% % Implementation 2: matrix multiplication
% py_unique = Distribution_Skellam_x0(repmat(y_unique,[1,numel(ps)]), repmat(sval.',[numel(y_unique),1]));
% py_unique = py_unique .* repmat(ps.',[numel(y_unique),1]);
% py_unique = mean(py_unique,2);

py = zeros(maxY-minY+1,1);
py(y_unique-minY+1) = py_unique;
py = py/sum(py);

% make py domain be symmetric
padsize = abs(maxY+minY);
if -minY < maxY
    py = [zeros(padsize,1); py];
else
    py = [py; zeros(padsize,1)];
end

% make py range be symmetric
py = make_symmetric(py);
end
