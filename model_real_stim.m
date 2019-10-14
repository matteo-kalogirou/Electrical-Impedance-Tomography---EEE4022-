
%
% Modelling the 32 electrode tank with an adjacent current drive. Current
% electrodes are odd numebred and voltage measuring electrodes are even
% numbered. Differential measurements are made on every adjacent pair of
% voltage measuring electrodes (256 measurements total). The input signal
% is a 1kHz sinusoidal current signal with amplitude 1mA.
%

% DOES NOT WORK AND I DO NOT KNOW WHY

clc;
clear;

nElec = 32;


% Tank forward model
fmdl = mk_circ_tank(16, [], nElec);


%%---Stimulation signal
Fs = 10000;                  % [samples/second]
dt = 1/Fs;                   % [seconds/sample]
num_samples = 10;           % [samples]
StopTime = dt*num_samples;   % seconds
t = (0:dt:StopTime-dt)';     % seconds
%%Sine wave:
Fc = 1000;                   % [Hz]
I = 0.001*sin(2*pi*Fc*t);          % Stimulation Signal


%%---Measurement Pattern
meas_pattern = zeros(32,32);
for i=2:2:30
    meas_pattern(i, i) = 1;
    meas_pattern(i, i+2) = -1;
end
meas_pattern(32, 32) = 1;
meas_pattern(32, 2) = -1;


%%--Create stim object
% 16 stimulations
j =1;

%%  
    SP = [];
    for k = 1:2:32
        sp = zeros(32,num_samples);
        sp(k, :) = I;
        if(k+2 < 32)
            sp(k+2,:) = -I;
        else
           sp(1,:) = -I;
        end
        SP = [SP, sp];
    end    
    
%%
    for i = 1:num_samples*16
        stim(i).stim_pattern = SP(:, i);
        stim(i).meas_pattern = meas_pattern;
        stim(i).stimulation = 'Amps';
    end

%%

fmdl.stimulation = stim;

%%

fmdl.solve = @fwd_solve_1st_order;
fmdl.system_mat = @system_mat_1st_order;
mdl_2d = eidors_obj('fwd_model', fmdl);


mat = ones( size(mdl_2d.elems, 1));
hom_img = eidors_obj('image', 'hom image', ...
                        'elem_data', mat, ...
                        'fwd_model', mdl_2d);

hom_data = fwd_solve(hom_img);          % takes very long to calculate
%%
show_fem(hom_img);
show_slices(hom_img);

%%

sl_fn = inline('(x.^2+(y+0.4).^2<0.2^2)', 'x', 'y', 'z');
inh_img = hom_img;

inh_img.elem_data = 1+ elem_select(inh_img.fwd_model, sl_fn);
inh_data = fwd_solve(inh_img);
%%
show_fem(inh_img);
 show_slices(inh_img);

%% SOLVE FWD PROBLEM


v_h = fwd_solve(hom_img);

v_inh= fwd_solve(inh_img);


%% DIFF Reconstruction

clear inv2d params mdl_2d_2;

% params = mk_circ_tank(16, [], nElec);
% 
% params.stimualtion = stim;
% params.solve = @fwd_solve_1st_order;
% params.system_mat = @system_mat_1st_order;
% params.jacobian = 'jacobian_adjoint';
% mdl_2d_2 = eidors_obj('fwd_model', params);
% 
% %--Inverse model
% inv2d.name = 'EIT inverse';
% inv2d.solve = @inv_solve_diff_GN_one_step;
% inv2d.hyperparameter.value= 3e-3;
% inv2d.RtR_prior = @prior_tikhonov;
% inv2d.reconst_type= 'difference';
% in2d.jacobian_bkgnd.value = 1;
% 
% %--- My way = fmdl; Sammy says mdl_2d_2
% inv2d.fwd_model = mdl_2d_2;

invmdl = eidors_obj('inv_model', 'EIT inverse model');
invmdl.reconst_type = 'difference';
invmdl.jacobian_bkgnd.value=1;
inv_model.solve = @inv_solve_diff_GN_one_step;
invmdl.hyperparameter.value = 3e-3;
invmdl.solve = @inv_solve_diff_GN_one_step;
invmdl.RtR_prior = @prior_tikhonov;
invmdl.fwd_model = fmdl;
invmdl.fwd_model.jacobian = @jacobian_adjoint;

%% reconstruct

img = inv_solve(invmdl, v_inh, v_h);
%%
show_slices(img);



