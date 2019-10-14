
%
%   Plot the data from the no Mux system
%
%   Inh = banana at 6 and 12 oclock
%
%
clc;
clear;

results_path = './NewResults/no_mux_hom_';  % double delimitted by ';'
for k=1:2
    % Convert all commas in files to decimal points
    
    for i=0:15
        comma2point_overwrite([results_path num2str(k) '/frame' num2str(i) '.txt']);        
    end

    % Read in the files
    for i = 1:16
       Frames(:,:,i) =  dlmread([results_path num2str(k) '/frame' num2str(i-1) '.txt'], ';', 0,0);
    end

    % Calculate the RMS value of each channel
    for i=1:16
       RmsFrames(:,:,i) = rms(Frames(:,:,i), 1);
    end      
end

hom_exp_data = reshape(RmsFrames, [256,1]);

ban_6oclk_data = importEITinput('./NewResults/no_mux_ban_6oclk/');

ban_12oclk_data = importEITinput('./NewResults/no_mux_ban_12oclk/');

ban_6oclk_2mA_data = importEITinput('./NewResults/no_mux_ban_6oclk_2mA/');

%% Plot Data

% 6 o'clock Data
figure()
subplot(211);plot(1:1:256,hom_exp_data, 1:1:256, ban_6oclk_data);
title('Homogenous and Inhomogeneous For Banana at 6 o`clock')
xlabel('Sample Number'); ylabel('RMS differential voltage');
legend('Hom', 'Inh');
set(gca, 'FontSize', 12);
axis([-1 257 0 0.4]);
subplot(212);plot(hom_exp_data - ban_6oclk_data);
title('Difference plot');
xlabel('Sample Number'); ylabel('RMS differential voltage');
axis([-1 257 -0.1 0.1]);
set(gca, 'FontSize', 12);

% 12 o'clock data
figure()
subplot(211);plot(1:1:256,hom_exp_data, 1:1:256, ban_12oclk_data);
title('Homogenous and Inhomogeneous For Banana at 12 o`clock')
xlabel('Sample Number'); ylabel('RMS differential voltage');
legend('Hom', 'Inh');
set(gca, 'FontSize', 12);
axis([-1 257 0 0.4]);
subplot(212);plot(hom_exp_data - ban_12oclk_data);
title('Difference plot');
xlabel('Sample Number'); ylabel('RMS differential voltage');
axis([-1 257 -0.1 0.1]);
set(gca, 'FontSize', 12);

%% % 2mA
figure()
subplot(211);plot(1:1:256,hom_exp_data, 1:1:256, ban_6oclk_2mA_data);
title('Homogenous and Inhomogeneous For Banana at 6 o`clock')
xlabel('Sample Number'); ylabel('RMS differential voltage');
legend('Hom', 'Inh');
set(gca, 'FontSize', 12);
axis([-1 257 0 0.5]);
subplot(212);plot(hom_exp_data - ban_6oclk_2mA_data);
title('Difference plot');
xlabel('Sample Number'); ylabel('RMS differential voltage');
axis([-1 257 -0.1 0.1]);
set(gca, 'FontSize', 12);

%% Reconstruct an image

imdl = mk_common_model('b2c', 16);
himg = mk_image(imdl,1);

options = {'meas_current', 'no_rotate_meas'};
[stim, meas_sel] = mk_stim_patterns(16, 1, '{ad}', '{ad}', options, 0.001);

himg.fwd_model.stimulation = stim;
himg.fwd_model.meas_select = meas_sel;

imdl.jacobian_bkgnd.value = 12;
imdl.solve = @inv_solve_diff_GN_one_step;
imdl.hyperparameter.value = 3e-3;
imdl.RtR_prior = @prior_tikhonov;

%%
figure();
img_ban_6 = inv_solve(imdl, hom_exp_data, ban_6oclk_data);
subplot(121); show_fem(img_ban_6);
img_ban_12 = inv_solve(imdl, hom_exp_data, ban_12oclk_data);
subplot(122); show_fem(img_ban_12);

%%

% Create Inverse Model
inv2d= eidors_obj('inv_model', 'EIT inverse');
inv2d.reconst_type= 'difference';
inv2d.jacobian_bkgnd.value= 1;

% This is not an inverse crime; inv_mdl != fwd_mdl
imb=  mk_common_model('b2c',16);
options = {'meas_current', 'no_rotate_meas'};
[stim, meas_sel] = mk_stim_patterns(16, 1, '{ad}', '{ad}', options, 0.001);
imb.fwd_model.stimulation = stim;
imb.fwd_model.meas_select = meas_sel;

inv2d.fwd_model= imb.fwd_model;

% Guass-Newton solvers
inv2d.solve=       @inv_solve_diff_GN_one_step;

% Tikhonov prior
inv2d.hyperparameter.value = 3e-3;
inv2d.RtR_prior=   @prior_tikhonov;
imgr(1)= inv_solve( inv2d, hom_exp_data, ban_6oclk_data);
imgb(1)= inv_solve( inv2d, hom_exp_data, ban_12oclk_data);

% NOSER prior
inv2d.hyperparameter.value = .1;
inv2d.RtR_prior=   @prior_noser;
imgr(2)= inv_solve( inv2d, hom_exp_data, ban_6oclk_data);
imgb(2)= inv_solve( inv2d, hom_exp_data, ban_12oclk_data);

% Laplace image prior
inv2d.hyperparameter.value = .1;
inv2d.RtR_prior=   @prior_laplace;
imgr(3)= inv_solve( inv2d, hom_exp_data, ban_6oclk_data);
imgb(3)= inv_solve( inv2d, hom_exp_data, ban_12oclk_data);
% 
% % Automatic hyperparameter selection
inv2d.hyperparameter = rmfield(inv2d.hyperparameter,'value');
inv2d.hyperparameter.func = @choose_noise_figure;
inv2d.hyperparameter.noise_figure= 0.5;
inv2d.hyperparameter.tgt_elems= 1:4;
inv2d.RtR_prior=   @prior_gaussian_HPF;
inv2d.solve=       @inv_solve_diff_GN_one_step;
imgr(4)= inv_solve( inv2d, hom_exp_data, ban_6oclk_data);
imgb(4)= inv_solve( inv2d, hom_exp_data, ban_12oclk_data);

inv2d.hyperparameter = rmfield(inv2d.hyperparameter,'func');

% Total variation using PDIPM
inv2d.hyperparameter.value = 1e-5;
inv2d.solve=       @inv_solve_TV_pdipm;
inv2d.R_prior=     @prior_TV;
inv2d.parameters.max_iterations= 10;
inv2d.parameters.term_tolerance= 1e-3;

%Vector of structs, all structs must have exact same (a) fields (b) ordering
imgr5= inv_solve( inv2d, hom_exp_data, ban_6oclk_data);
imgr5=rmfield(imgr5,'type'); imgr5.type='image';
imgr(5)=imgr5;
 
imgb5= inv_solve( inv2d, hom_exp_data, ban_12oclk_data);
imgb5=rmfield(imgb5,'type'); imgb5.type='image';
imgb(5)=imgb5;



subplot(251); show_fem(imgr(1)); set(gca, 'visible', 'off');
subplot(252); show_fem(imgr(2)); set(gca, 'visible', 'off');
subplot(253); show_fem(imgr(3)); set(gca, 'visible', 'off');
subplot(254); show_fem(imgr(4)); set(gca, 'visible', 'off');
subplot(255); show_fem(imgr(5)); set(gca, 'visible', 'off');

subplot(256); show_fem(imgb(1)); set(gca, 'visible', 'off');
subplot(257); show_fem(imgb(2)); set(gca, 'visible', 'off');
subplot(258); show_fem(imgb(3)); set(gca, 'visible', 'off');
subplot(259); show_fem(imgb(4)); set(gca, 'visible', 'off');
subplot(2,5,10); show_fem(imgb(5)); set(gca, 'visible', 'off');


