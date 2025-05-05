clear;            
clc;              

lam = 1550e-9;    
a = 1000e-9;      
d = 2000e-9;      
N = 128;          
k = 2*pi/lam;     % k = 2π / λ
p = 100;            
phi = -4.05*pi/180 : 0.01*pi/180 : 4.05*pi/180;  % scanning angle - phi

for i = 1 : 500

neff = 1.56;   
nct = 1;       
lamx = 100*10^(-9)*(rand(1,1)-0.5) ;  
p = sin(N/2*(k*d.*sin(phi) - k*neff*lamx)) ./ sin(1/2*(k*d.*sin(phi) - k*neff*lamx));
p = abs(p);
P = N^2 .* (sinc((pi.*a)./lam.*sin(phi)).^2) .* (p.^2);

angle = rad2deg(phi);
figure(1);
plot(angle, P);      
xlabel('Position');  
ylabel('Intensity'); 

end
% figure(2);
% plot(angle,p);      
% xlabel('Position'); 
% ylabel('Intensity');
% aaaa = P';