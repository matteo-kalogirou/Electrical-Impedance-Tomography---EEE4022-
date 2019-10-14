
%
%
%   Hom - Inh signal analysis
%
%



run import_hom_exp_data;

p = '/Users/matteokalogirou/Documents/MATLAB/Electrical Impedance Tomography - EEE4022/Results';

p1 = [p '/1kHz'];
p2 = [p '/300Hz'];

%% --- 1Khz
% Load inhomogenoud data
inh_exp_left = importEITinput([p1 '/inh_test_cylinder_left/']);
inh_exp_right = importEITinput([p1 '/inh_test_cylinder_right/']);
inh_exp_top = importEITinput([p1 '/inh_test_cylinder_top/']);
inh_exp_bottom = importEITinput([p1 '/inh_test_cylinder_bottom/']);


%% --- 300Hz

% Load inhomogenous data
inh_exp_c1 = importEITinput([p2 '/inh_test_cylinder/']);
inh_exp_c2 = importEITinput([p2 '/inh_test_cylinder2/']);
inh_exp_2c1 = importEITinput([p2 '/inh_test_2cylinders1/']);
inh_exp_2c2 = importEITinput([p2 '/inh_test_2cylinders2/']);
inh_exp_gb = importEITinput([p2 '/inh_test_glassbottle/']);

%%

figure()
plot(1:1:256, hom_exp_data_1k, 1:1:256, inh_exp_top); title('1k Hom');

%%
figure()
subplot(5,1,1); plot(hom_exp_data_1k); title('1k Hom');
subplot(5,1,2); plot(hom_exp_data_1k-inh_exp_left); title('1k left');
subplot(5,1,3); plot(hom_exp_data_1k-inh_exp_right); title('1k right');
subplot(5,1,4); plot(hom_exp_data_1k-inh_exp_top); title('1k top');
subplot(5,1,5); plot(hom_exp_data_1k-inh_exp_bottom); title('1k bottom');
%%

figure()
subplot(6,1,1); plot(hom_exp_data_300); title('Hom 300Hz'); axis([0 256 0 1]);
subplot(6,1,2); plot(hom_exp_data_300-inh_exp_c1); title('hom - inh 1'); axis([0 256 0 1]);
subplot(6,1,3); plot(hom_exp_data_300-inh_exp_c2); title('hom - inh 2'); axis([0 256 0 1]);
subplot(6,1,4); plot(hom_exp_data_300-inh_exp_2c1); title('hom - inh 3'); axis([0 256 0 1]);
subplot(6,1,5); plot(hom_exp_data_300-inh_exp_2c2); title('hom - inh 4'); axis([0 256 0 1]);
subplot(6,1,6); plot(hom_exp_data_300-inh_exp_gb); title('hom - inh 4'); axis([0 256 0 1]);


%% New Results


p = '/Users/matteokalogirou/Documents/MATLAB/Electrical Impedance Tomography - EEE4022/NewResults/';
el = 1:1:256;

hom_data = importEITinput([p 'hom_1k_1vpp_1/']);
hom_data = [hom_data importEITinput([p 'hom_1k_1vpp_2/'])] ;
hom_data = [hom_data importEITinput([p 'hom_1k_1vpp_3/'])] ;
hom_data = mean(hom_data, 2);


b6 = importEITinput([p 'banana_6oclk_1k_1vpp/']);
b12 = importEITinput([p 'banana_12oclk_1k_1vpp/']);

figure()
subplot(1,2,1); 
plot(el, hom_data, el, b6); axis tight
xlabel('Measurement number');
ylabel('\Delta V [V]');
set(gca, 'FontSize', 14);
legend('Hom', 'Inh', 'location', 'East');
axis([1 256 -0.1 0.4])
subplot(1,2,2); plot(el, hom_data-b6); axis tight
xlabel('Measurement number');
ylabel('\Delta V [V]');
set(gca, 'FontSize', 14);
axis([1 256 -0.1 0.4])
%%
figure()
% subplot(1,2,1);
plot(el, hom_data, el, b12); axis tight
xlabel('Measurement number');
ylabel('\Delta V [V]');
set(gca, 'FontSize', 14);
legend('Homogeneous', 'Inhomogeneous', 'location', 'south');
axis([1 256 -0.1 0.4])
% subplot(1,2,2); 
figure()
plot(el, hom_data-b12); axis tight
xlabel('Measurement number');
ylabel('\Delta V [V]');
set(gca, 'FontSize', 14);
axis([1 256 -0.1 0.4])





