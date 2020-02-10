clear all;
close all;
clc;

x = double(imread('Chicago.png'));
x(675,:) = [];
%figure(1);
%imshow(x/255);

x = mean(x,3);

x = x*10; % fix me!
% figure(2);
% imshow(x/255);

m1 = randn(size(x));
n1 = randn(size(x));
sigma_m = input("What is the value for sigma_m? ");
sigma_n = input("What is the value for sigma_n? ");

y1 = x + m1*sigma_m;
z1 = x + n1*sigma_n;

center_y = size(y1)/2;

smaller_x = x(center_y(1) + [-40:40] + 25, center_y(2) + [-40:40] + 200);
figure(5);
imshow(smaller_x/max(smaller_x(:)));

smaller_y = y1(center_y(1) + [-40:40] + 25, center_y(2) + [-40:40] + 200);
figure(6);
imshow(smaller_y/max(smaller_x(:)));

smaller_z = z1(center_y(1) + [-40:40] + 25, center_y(2) + [-40:40] + 200);
figure(7);
imshow(smaller_z/max(smaller_x(:)));

patch_size = 4;
%sig = 700;
sig = 700;
w = zeros(size(smaller_y));

for i = patch_size+1:(size(smaller_y,1)-patch_size)
    for j = patch_size+1:(size(smaller_y,1) - patch_size)
        numerator = 0;
        denominator = 0;
        y_patch = smaller_y(i + [-patch_size:patch_size], j + [-patch_size:patch_size]);
        for m = patch_size+1:(size(y1,1)-patch_size)
            for n = patch_size+1:(size(y1,1) - patch_size)
                comparing_patch = y1(m + [-patch_size:patch_size], n + [-patch_size:patch_size]);
                patch_diff = sum((y_patch(:) - comparing_patch(:)).^2);
                patch_phi = exp(-patch_diff/(2*(sig^2)));
                numerator = numerator + patch_phi*z1(m,n);
                denominator = denominator + patch_phi;
            end
        end
        w(i,j) = numerator/denominator;
    end
end

figure(10);
imshow(w/max(smaller_x(:)));

%figure(3);
%imshow(y_patch/max(y_patch(:)));
        
% y1 = reshape(poissrnd(x(:)),size(x));
% z1 = reshape(poissrnd(x(:)),size(x));

%sigma_n2_check1 = mean((z1(:) - x(:)).^2);

%avg_sigma_n2_check = (sigma_n2_check1 + sigma_n2_check2 + sigma_n2_check3)/3;% + sigma_n2_check4 + sigma_n2_check5 + sigma_n2_check6 + sigma_n2_check7 + sigma_n2_check8 + sigma_n2_check9 + sigma_n2_check10)/10;

%sigma_m2_check = mean((y(:) - x(:)).^2);

%sigma_m_check = sqrt(sigma_m2_check);
%sigma_n_check = sqrt(sigma_n2_check);

% figure(3);
% imshow(y/255);
% figure(4);
% imshow(z/255);

%% Anscombe transform

%xa = 2.*sqrt(x + (3/8));
%ya = 2.*sqrt(y + (3/8));
%za = 2.*sqrt(z + (3/8));

%% Phi Calculations

% User input
  % phi1 = input("What does phi represent (ex: y, y.^2, sqrt(y), etc)? ");
  phi1 = smaller_y;

% BM3D denoising method
% third parameter is sigma_m for regular BM3D, but 1 for Anscombe transform
%[PSNR1,BM3D_phi1] = BM3D(x, y1, 1, 'np', 0);
% BM3D_phi1 = BM3D_phi1 * 255;
%% inverse Anscombe
%BM3D_phi = (BM3D_phi/2).^2 - (3/8);

%% Other denoising methods

% DNN denoising method
%net = denoisingNetwork('DnCNN');
%DNN_phi1 = denoiseImage(y1/255,net);
%DNN_phi1 = DNN_phi1 * 255;


% Skellam Shrink denoising method
% waveLvl = 2; %waveLvl = 5 produces the best MSE?
% skellam_phi = ske_mrso(y, waveLvl);
%% Mean square error

MSE_result1 = MSE(smaller_x,phi1);

%avg_MSE = (MSE_result1 + MSE_result2)/2;% + MSE_result3 + MSE_result4 + MSE_result5 + MSE_result6 + MSE_result7 + MSE_result8 + MSE_result9 + MSE_result10)/10;
%% Changed mean square error

[MYMSE_result1, sigma_n_squared1, sigma_m_squared1] = MYMSE(smaller_y, smaller_z,phi1);

%avg_sigma_n_squared = (sigma_n_squared1 + sigma_n_squared2 + sigma_n_squared3)/3;% + sigma_n_squared4 + sigma_n_squared5 + sigma_n_squared6 + sigma_n_squared7 + sigma_n_squared8 + sigma_n_squared9 + sigma_n_squared10)/10;
%avg_MYMSE = (MYMSE_result1 + MYMSE_result2 + MYMSE_result3)/3;% + MYMSE_result4 + MYMSE_result5 + MYMSE_result6 + MYMSE_result7 + MYMSE_result8 + MYMSE_result9 + MYMSE_result10)/10;

%% Functions

function [result] = MSE(x,phi)
    result = mean((x(:) - phi(:)).^2);
end

function [result2, sigma_n_squared, sigma_m_squared] = MYMSE(y,z,phi)
    %sigma_z_squared = mean((z(:) - mean(z(:))).^2);
    %sigma_y_squared = mean((y(:) - mean(y(:))).^2);
    sigma_n_squared = (mean(z(:).^2) - mean(y(:).^2) + mean((y(:) - z(:)).^2))/2;
    sigma_m_squared = (mean(y(:).^2) - mean(z(:).^2) + mean((y(:) - z(:)).^2))/2;
    result2 = mean((z(:) - phi(:)).^2) - sigma_n_squared;
end