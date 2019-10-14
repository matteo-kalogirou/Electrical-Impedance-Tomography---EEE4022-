%
% Create a 32 electrode model of the experiental tank. System has 16
% stimulations and current injecting elctrodes are made on odd numbered
% electrodes. Measurements are made on even electrodes. An adjacent type
% drive is used, but only between adjacent current electrodes.
%
%


clear hom_img hom_data imdl fmdl stim meas_sel

n_elec = 32;
I = 0.001;


%---MEASUREMENT PATTERN -> 16x32 array with 1 and -1
meas_sel = zeros(16,32);
for i = 1:16
    meas_sel(i,i*2) = 1;
    if(i == 16)
       meas_sel(i,2) = -1;
    else 
       meas_sel(i,2*i+2) = -1;
    end
end


%%
%---STIMUATION PATTERN
for i=1:16
    s = zeros(32,1);   
    if(i*2-1 < 31)
        s(i*2-1) = I;                        
        s(i*2+1) = -I;
    else
        s(i*2-1) = I;                        
        s(1) = -I;
    end

    stim(i).stimulation = 'A';
    stim(i).stim_pattern = s;
    stim(i).meas_pattern = meas_sel;
end


%% --- FWD MODEL

fmdl = mk_circ_tank(8, [], 32);
fmdl.stimulation = stim;
fmdl.solve = @fwd_solve_1st_order;
fmdl.system_mat = @eidors_default;
fmdl.jacobian = @eidors_default;

if(~valid_fwd_model(fmdl))
    display('INVALID FWD MODEL');
end

%% --Homogeneous Image

hom_img = mk_image(fmdl, 1);    %constant conductivity ~ 55 S/m for saline
hom_img.fwd_solve.get_all_meas = 1;

hom_data = fwd_solve(hom_img);

% Simulate Inhomogeneous Object
t_circ = makeCircle(); %create random circle object
circle = @(x,y,z) (x-t_circ.x).^2 + (y-t_circ.y).^2 < t_circ.r.^2;
inh_img = hom_img;
inh_img.elem_data = hom_img.elem_data + elem_select(inh_img.fwd_model,circle);

inh_data = fwd_solve(inh_img);

%% ------------------------Plotting 
figure();
plot(1:1:256,hom_data.meas, '-r', 1:256,inh_data.meas, '-b'); axis tight;
ylabel('Voltage');
xlabel('Measurement number');
legend('Homogeneous', 'Inhomogeneous', 'location', 'West');
set(gca, 'FontSize', 14);

figure();
plot(hom_data.meas - inh_data.meas); axis tight;
ylabel('Differential Measurement [V]');
xlabel('Measurement number');
set(gca, 'FontSize', 14);

%%
% figure();
% [hax, hL] = plotyy(1:256, [hom_data.meas inh_data.meas] , ... 
%            1:256, hom_data.meas - inh_data.meas);
% hL(1).Color = 'magenta';
% hL(2).Color = 'blue';
% axis tight;
% ylabel('Voltage');
% xlabel('Measurement number');
% legend('Homogeneous', 'Inhomogeneous', 'location', 'West');
% set(gca, 'FontSize', 14);


%% --- INV MODEL

imdl = eidors_obj('inv_model', '32e EIT Inverse model');
imdl.fwd_model = fmdl;
imdl.reconst_type = 'difference';
imdl.jacobian_bkgnd.value = 1;
imdl.solve = @inv_solve_diff_GN_one_step;
imdl.hyperparameter.value = 3e-3;
imdl.RtR_prior = @prior_tikhonov;

%% --- IMAGE RECONSTRUCTION

rimg = inv_solve(imdl, hom_data, inh_data);

subplot(121); show_fem(inh_img); title('Inhomogenous Image'); axis off;
subplot(122); show_fem(rimg); title('Reconstructed Image'); axis off;

%% ========================================================================
%              --- RECONSTRUCT OLD EXPERIMENTAL DATA ---
%==========================================================================

% Load DATA
run import_hom_exp_data

p = '/Users/matteokalogirou/Documents/MATLAB/Electrical Impedance Tomography - EEE4022/Results';
f = '/1kHz/';
p = [p f];

%1Khz
inh_data_left = importEITinput([p '/inh_test_cylinder_left/']);
inh_data_right = importEITinput([p '/inh_test_cylinder_right/']);
inh_data_top = importEITinput([p '/inh_test_cylinder_top/']);
inh_data_bottom = importEITinput([p '/inh_test_cylinder_bottom/']);

%300Hz
p = '/Users/matteokalogirou/Documents/MATLAB/Electrical Impedance Tomography - EEE4022/Results';
f = '/300Hz/';
p = [p f];
inh_data_cylinder = importEITinput([p '/inh_test_cylinder/']);
inh_data_cylinder2 = importEITinput([p '/inh_test_cylinder2/']);
inh_data_2cylinder1 = importEITinput([p '/inh_test_2cylinders1/']);
inh_data_2cylinder2 = importEITinput([p '/inh_test_2cylinders2/']);
inh_data_glassbottle = importEITinput([p '/inh_test_glassbottle/']);

% SOLVE
img_left        = inv_solve(imdl, hom_exp_data_1k, inh_data_left);
img_right       = inv_solve(imdl, hom_exp_data_1k, inh_data_right);
img_top         = inv_solve(imdl, hom_exp_data_1k, inh_data_top);
img_bottom      = inv_solve(imdl, hom_exp_data_1k, inh_data_bottom);
img_cylinder    = inv_solve(imdl, hom_exp_data_300, inh_data_cylinder);
img_cylinder2   = inv_solve(imdl, hom_exp_data_300, inh_data_cylinder2);
img_2cylinder1  = inv_solve(imdl, hom_exp_data_300, inh_data_2cylinder1);
img_2cylinder2  = inv_solve(imdl, hom_exp_data_300, inh_data_2cylinder2);
img_glassbottle = inv_solve(imdl, hom_exp_data_300, inh_data_glassbottle);

figure(); title('1 KHz');
h11 = subplot(1,4,1); show_fem(img_top); title('Top');
h12 = subplot(1,4,2); show_fem(img_left); title('Left');
h13 = subplot(1,4,3); show_fem(img_right); title('Right');
h14 = subplot(1,4,4); show_fem(img_bottom); title('Bottom');


figure(7); title('300Hz');
subplot(3,3,2); show_fem(img_cylinder); title('Cylinder');
subplot(3,3,4); show_fem(img_cylinder2); title('Cylinder');
subplot(3,3,8); show_fem(img_2cylinder1); title('2 Cylinders');
subplot(3,3,6); show_fem(img_2cylinder2); title('2 Cylinders');
subplot(3,3,8); show_fem(img_glassbottle); title('Glass Bottle');



%% ========================================================================
%              --- RECONSTRUCT NEW EXPERIMENTAL DATA ---
%==========================================================================

Hom_path = './NewResults/hom_1k_1vpp_';
b6_path = './NewResults/banana_6oclk_1k_1vpp/';
b12_path = './NewResults/banana_12oclk_1k_1vpp/';

% Load Hom data
rmsAmp = [];
for j = 1:3
    for i=0:15
        comma2point_overwrite([Hom_path num2str(j) '/frame' num2str(i) '.txt']);
        Frames(:,:,i+1) =  dlmread([Hom_path num2str(j) '/frame' num2str(i) '.txt'], ';', 0,0);
        RmsFrames(:,:,i+1) = rms(Frames(:,:,i+1), 1);
    end
    x = reshape(RmsFrames, [1,256]);
    rmsAmp = [rmsAmp ; x];
end

Hom_data = mean(rmsAmp);

% Load Inh data
b6 = importEITinput(b6_path);
b12 = importEITinput(b12_path);

% Reconstruct Image 3e-6 good
imdl.hyperparameter.value = 4e-6; 

img_b6 = inv_solve(imdl, Hom_data', b6);
img_b12 = inv_solve(imdl, Hom_data', b12);
   
figure(2);
subplot(1,2,1); show_fem(img_b6); axis off; title('Banana at 6');
subplot(1,2,2); show_fem(img_b12); axis off; title('Banana at 12');

%%
imdl.hyperparameter.value = 4e-6; 
figure(); show_fem(img_b12, [0,1,0]); axis off;
    
%% ========================================================================
%              --- RECONSTRUCT NEW EXPERIMENTAL DATA ---
%==========================================================================

Hom_path = './NewResults/no_mux_hom_';
b6_path = './NewResults/no_mux_ban_6oclk/';
b12_path = './NewResults/no_mux_ban_12oclk/';

% Load Hom data
rmsAmp = [];
for j = 1:2
    for i=0:15
        comma2point_overwrite([Hom_path num2str(j) '/frame' num2str(i) '.txt']);
        Frames(:,:,i+1) =  dlmread([Hom_path num2str(j) '/frame' num2str(i) '.txt'], ';', 0,0);
        RmsFrames(:,:,i+1) = rms(Frames(:,:,i+1), 1);
    end
    x = reshape(RmsFrames, [1,256]);
    rmsAmp = [rmsAmp ; x];
end

Hom_data = mean(rmsAmp);
%SCALEING
Hom_data = Hom_data;

% Load Inh data
b6 = importEITinput(b6_path);
b12 = importEITinput(b12_path);

%  Reconstruct Image

imdl.hyperparameter.value = 5e-7; 

img_b6 = inv_solve(imdl, Hom_data', b6');
img_b12 = inv_solve(imdl, Hom_data', b12');
   
figure(3);
subplot(1,2,1); show_fem(img_b6); axis off; title('Banana at 6');
subplot(1,2,2); show_fem(img_b12); axis off; title('Banana at 12');
    




