% synthetic image denoising demo
addpath(genpath('support'));
wavLvl = 4;
im = double(imread('test_image/cafe.png'));

% generate Poisson noise
targetMean = 1.6;
scale = 1/mean(im(:))*targetMean;
imNsy = poissrnd(im*scale); 

% denoising
[fhat1,fhat2,fhat3] = ske_mrso(imNsy,wavLvl);
fhat1 = fhat1/scale;
fhat2 = fhat2/scale;
fhat3 = fhat3/scale;

savDir = 'demo1_output/';
if ~exist(savDir,'dir'); mkdir(savDir); end;
imwrite(im/255,[savDir,'input.png'],'png');
imwrite(fhat1/255,[savDir,'HMRSO.png'],'png');
imwrite(fhat2/255,[savDir,'BMRSO.png'],'png');
imwrite(fhat3/255,[savDir,'UMRSO.png'],'png');