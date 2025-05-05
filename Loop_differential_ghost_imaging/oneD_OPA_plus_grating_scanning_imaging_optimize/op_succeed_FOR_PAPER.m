clc;
clear;

obj = imread('littleGI64.tif');  
n = size(obj);                 
if size(obj, 3) == 3
    obj = rgb2gray(obj);        
end
obj = im2double(obj);           
m = n(1);                       
sampling_rates = [0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.75, 1.0];  
num_rates = length(sampling_rates); 

% Initialization
lam = 1550 * 10^(-9);            
A = 1;                          
Nx = 128;                     
Ny = 128;                     
neff = 1.56;                    
Zx = 2000 * 10^-9;              
Zy = 2000 * 10^-9;             
a1 = 1000 * 10^-9;           
a2 = 1000 * 10^-9;               
L1 = 0.1413;                   
kk_max = 5000;                
kkO = 10;                   
LDGI_results = cell(1, num_rates);

for rate_idx = 1:num_rates
    sampling_rate = sampling_rates(rate_idx); 
    kk = round(kk_max * sampling_rate);        
    c = zeros(n(1), n(1));  
    c1 = zeros(n(1), n(1)); 
    sjsball = zeros(n(1), n(1)); 
    B = [];  

    fprintf('Processing sampling rate: %.0f%%\n', sampling_rate * 100);

   
    for i = 1:round(kk/kkO)
        lamx = 100 * 10^(-9) * (rand(1,1) - 0.5); 
        for O = 1:kkO
            sjjz = opa_grating_scan(lamx, lam, A, Nx, Ny, neff, Zx, Zy, a1, a2, L1, m);
            obj_after_speckle = sjjz .* obj;      
            total_intensity = sum(obj_after_speckle(:)); 
            c = c + total_intensity * sjjz; 
            B(O) = total_intensity;           
            sjsball = sjsball + sjjz;        
        end
        % LDGI
        c1 = c1 + (c / kkO) - (sjsball / kkO .* sum(B) / kkO);
        sjsball = zeros(n(1), n(1));
        c = zeros(n(1), n(1));
    end
    LDGI_results{rate_idx} = c1;
end
ref_image = obj;

psnr_values_LDGI = zeros(1, num_rates);
ssim_values_LDGI = zeros(1, num_rates);
rmse_values_LDGI = zeros(1, num_rates);
for rate_idx = 1:num_rates
    current_image = abs(LDGI_results{rate_idx});
    current_image_norm = (current_image - min(current_image(:))) / (max(current_image(:)) - min(current_image(:)));
    ref_image_norm = (ref_image - min(ref_image(:))) / (max(ref_image(:)) - min(ref_image(:)));
    psnr_values_LDGI(rate_idx) = psnr(current_image_norm, ref_image_norm);
    ssim_values_LDGI(rate_idx) = ssim(current_image_norm, ref_image_norm);
    rmse_values_LDGI(rate_idx) = sqrt(mean((current_image_norm(:) - ref_image_norm(:)).^2));
end

figure(1);
for rate_idx = 1:num_rates
    subplot(1, num_rates, rate_idx);
    imagesc(LDGI_results{rate_idx});   
    colormap('gray');
    title(sprintf('LDGI (%.0f%%)\nPSNR: %.2f\nSSIM: %.3f\nRMSE: %.4f', ...
          sampling_rates(rate_idx) * 100, psnr_values_LDGI(rate_idx), ...
          ssim_values_LDGI(rate_idx), rmse_values_LDGI(rate_idx)));
    xlabel('X Position');
    ylabel('Y Position');
    colorbar;
end

figure(2);
subplot(3,1,1);
plot(sampling_rates * 100, psnr_values_LDGI, '-o');
title('PSNR vs Sampling Rate');
xlabel('Sampling Rate (%)');
ylabel('PSNR (dB)');

subplot(3,1,2);
plot(sampling_rates * 100, ssim_values_LDGI, '-o');
title('SSIM vs Sampling Rate');
xlabel('Sampling Rate (%)');
ylabel('SSIM');

subplot(3,1,3);
plot(sampling_rates * 100, rmse_values_LDGI, '-o');
title('RMSE vs Sampling Rate');
xlabel('Sampling Rate (%)');
ylabel('RMSE');

function sjjz = opa_grating_scan(lamx, lam, A, Nx, Ny, neff, Zx, Zy, a1, a2, L1, m)
    k = 2 * pi / lam;
    dx1 = L1 / m;
    thix = -L1/2 : dx1 : (L1/2 - dx1);  
    thiy = thix';                      
    Eall = 0;
    for h = 0:Ny-1
        E = exp(1i * (k * h * Zy .* sin(thiy) - (rand(1, 1) - 0.5) * 2 * pi));
        Eall = Eall + E;
    end
    Eall = Eall .* sin(Nx / 2 * (k * Zx .* sin(thix) - k * neff * lamx)) ./ ...
                   sin(1 / 2 * (k * Zx .* sin(thix) - k * neff * lamx));
    Eall = abs(Eall);
    T = A * Nx * Ny .* ((sinc((pi * a1) / lam .* sin(thix))) .* ...
                        (sinc((pi * a2) / lam .* sin(thiy))) .* Eall).^2;
    sjjz = abs(T);
end