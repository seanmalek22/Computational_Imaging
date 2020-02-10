% sensor image denoising demo
addpath(genpath('support'));
wavLvl = 5;
load('test_image/testImageRaw','imRaw','black');

% denoising
[fhat1,fhat2,fhat3] = ske_mrso(imRaw,wavLvl,40);
im = imRaw-black;
fhat1 = fhat1-black;
fhat2 = fhat2-black;
fhat3 = fhat3-black;

scale = median(fhat1(:))*4; % scale for display
savDir = 'demo2_output/';
if ~exist(savDir,'dir'); mkdir(savDir); end;
imwrite(im/scale,[savDir,'input.png'],'png');
imwrite(fhat1/scale,[savDir,'HMRSO.png'],'png');
imwrite(fhat2/scale,[savDir,'BMRSO.png'],'png');
imwrite(fhat3/scale,[savDir,'UMRSO.png'],'png');