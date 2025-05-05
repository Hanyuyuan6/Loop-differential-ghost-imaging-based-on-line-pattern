clc;
clear;

obj = imread("littleGI64.tif");   
n = size(obj);                    
if size(obj, 3) == 3
    obj = rgb2gray(obj);  
end
obj = im2double(obj);        
m = n(1);                     

% Initialization
max_k = 5000;                     
sampling_rates = [0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.75, 1.0];  
num_rates = length(sampling_rates);    
c_corrs = cell(1, num_rates);     
c_diffs = cell(1, num_rates);    


PSNR_corr = zeros(1, num_rates);   % PSNR for correlation imaging
SSIM_corr = zeros(1, num_rates);   % SSIM for correlation imaging
RMSE_corr = zeros(1, num_rates);   % RMSE for correlation imaging

PSNR_diffs = zeros(1, num_rates);  % PSNR for differential imaging
SSIM_diffs = zeros(1, num_rates);  % SSIM for differential imaging
RMSE_diffs = zeros(1, num_rates);  % RMSE for differential imaging

% Ground truth 
ground_truth = obj;

lamx_perturbation = 100 * 10^(-9) * (rand(1, round(max_k / 10)) - 0.5);
for rate_idx = 1:num_rates
    k = round(max_k * sampling_rates(rate_idx)); 
    c_corr = zeros(m, m);           
    sjsball = zeros(m, m);          
    BB = zeros(1, k);           
    B = zeros(1, 10);                

    fprintf('Processing sampling rate: %.0f%%\n', sampling_rates(rate_idx) * 100);

    for i = 1:round(k/10)
        for O = 1:10
            sjjz = opa_grating_scan(lamx_perturbation(i), 1550, 1, 128, 128, 1.56, 2000, 2000, 1000, 1000, 0.1413, m);
            in1 = sjjz .* obj;          
            in2 = sum(in1(:));           
            c_corr = c_corr + in2 * sjjz;
            B(O) = in2;           
            sjsball = sjsball + sjjz;
        end
        BB(10*(i-1)+1:10*i) = B; 
    end
    c_corr = abs(c_corr);
    c_diff = abs((c_corr / k) - (sjsball / k .* sum(BB) / k));
    c_corr_norm = (c_corr - min(c_corr(:))) / (max(c_corr(:)) - min(c_corr(:)));
    c_diff_norm = (c_diff - min(c_diff(:))) / (max(c_diff(:)) - min(c_diff(:)));
    ground_truth_norm = (ground_truth - min(ground_truth(:))) / (max(ground_truth(:)) - min(ground_truth(:)));

    c_corrs{rate_idx} = c_corr_norm;
    c_diffs{rate_idx} = c_diff_norm;

    PSNR_corr(rate_idx) = psnr(c_corr_norm, ground_truth_norm);
    SSIM_corr(rate_idx) = ssim(c_corr_norm, ground_truth_norm);
    RMSE_corr(rate_idx) = sqrt(mean((c_corr_norm(:) - ground_truth_norm(:)).^2));

    PSNR_diffs(rate_idx) = psnr(c_diff_norm, ground_truth_norm);
    SSIM_diffs(rate_idx) = ssim(c_diff_norm, ground_truth_norm);
    RMSE_diffs(rate_idx) = sqrt(mean((c_diff_norm(:) - ground_truth_norm(:)).^2));
end


figure;
for rate_idx = 1:num_rates
    subplot(2, num_rates, rate_idx);
    imagesc(c_corrs{rate_idx});       
    colormap('gray');
    title(sprintf('Corr Imaging (%.0f%%)', sampling_rates(rate_idx) * 100));
    xlabel('X Position');
    ylabel('Y Position');
    colorbar;

    subplot(2, num_rates, num_rates + rate_idx);
    imagesc(c_diffs{rate_idx});      
    colormap('gray');
    title(sprintf('Diff Imaging (%.0f%%)', sampling_rates(rate_idx) * 100));
    xlabel('X Position');
    ylabel('Y Position');
    colorbar;
end

figure;
subplot(3, 1, 1);
plot(sampling_rates * 100, PSNR_corr, '-o', 'LineWidth', 2);
hold on;
plot(sampling_rates * 100, PSNR_diffs, '-o', 'LineWidth', 2);
xlabel('Sampling Rate (%)');
ylabel('PSNR');
legend('Correlation Imaging', 'Differential Imaging');
title('PSNR vs Sampling Rate');
grid on;

subplot(3, 1, 2);
plot(sampling_rates * 100, SSIM_corr, '-o', 'LineWidth', 2);
hold on;
plot(sampling_rates * 100, SSIM_diffs, '-o', 'LineWidth', 2);
xlabel('Sampling Rate (%)');
ylabel('SSIM');
legend('Correlation Imaging', 'Differential Imaging');
title('SSIM vs Sampling Rate');
grid on;

subplot(3, 1, 3);
plot(sampling_rates * 100, RMSE_corr, '-o', 'LineWidth', 2);
hold on;
plot(sampling_rates * 100, RMSE_diffs, '-o', 'LineWidth', 2);
xlabel('Sampling Rate (%)');
ylabel('RMSE');
legend('Correlation Imaging', 'Differential Imaging');
title('RMSE vs Sampling Rate');
grid on;

% generate OPA speckle patterns
function sjjz = opa_grating_scan(lamx, lam, A, Nx, Ny, neff, Zx, Zy, a1, a2, L1, m)

lam = lam * 10^-9;              
k = 2 * pi / lam;
Zx = Zx * 10^-9;         
Zy = Zy * 10^-9;             
a1 = a1 * 10^-9;           
a2 = a2 * 10^-9;          

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