%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Matteo Kalogirou
%   2019
%   EEE4022 - EIT
%
%   This script is used to generate the simualted EIT data to be used in
%   training a neural network for the reconstruction process.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
close all;
clear;

run /Users/matteokalogirou/Documents/MATLAB/eidors-v3.9.1-ng/eidors/startup.m

%% --- Modelling
% This model represents the experimental setup, which is a single plane 
% 16/16-electrode, adjacent current injection (1mA) and measurement system.

num_trials = 1;



%% ========================================================================
%%========================= CUSTOM STIMULAITON ============================
%%=========================================================================
clear hom_img hom_data inv_model mdl

n_elec = 32;
I = 0.001;
n_meas = 256;
load('sp_mp.mat');      % Custom stimulation and measurement pattern for tank

stim = stim_meas_list(sp_mp, n_elec, I, 1);

for i =1:size(stim, 2)
   stim(i).meas_pattern = -stim(i).meas_pattern;
end

mdl = mk_common_model('b2c', n_elec);
z_model = 1;
hom_img = mk_image(mdl, z_model);
%show_fem(hom_img, [0, 1, 1]);

hom_img.fwd_model.stimulation = stim;
hom_img.fwd_model.meas_select = true(n_meas,1);
hom_img.fwd_solve.get_all_meas = 1;

% FWD SOLVE
hom_data = fwd_solve(hom_img);

%Create inverse model same as forward
inv_model = eidors_obj('inv_model', 'EIT_inverse');
f_mdl = mk_common_model('b2c', 16);
inv_model.fwd_model = f_mdl.fwd_model;

%Solver parameters
inv_model.reconst_type = 'difference'; %difference imaging
inv_model.jacobian_bkgnd.value = 1;
inv_model.solve = @inv_solve_diff_GN_one_step;
inv_model.hyperparameter.value = 3e-3;
inv_model.RtR_prior = @prior_tikhonov;

%% Test
    test_img = hom_img;
    t_circ = makeCircle();
    circle = @(x,y,z) (x-t_circ.x).^2 + (y-t_circ.y).^2 < t_circ.r.^2;

    test_img.elem_data = hom_img.elem_data + elem_select(test_img.fwd_model,circle);
    test_data = fwd_solve(test_img);

    subplot(121);show_fem(test_img); title('Simulated Inhomogeneity');    
    subplot(122);show_fem(inv_solve(inv_model, hom_data, test_data), [0,0,0]); title('Reconstructed Image')



%% ========================================================================
%%====================== 16 Electrode Adjacent ============================
%%=========================================================================
clear hom_img hom_data inv_model mdl

n_elec = 16;
n_rings = 1;
I = 0.001;
n_meas = 256;

mdl = mk_common_model('b2c', 16);
z_model = 1;
hom_img = mk_image(mdl,z_model);

options = {'no_meas_current', 'no_rotate_meas'};
[stim, meas_sel] = mk_stim_patterns(n_elec, n_rings, '{ad}', '{ad}', options, I);

hom_img.fwd_model.stimulation = stim;
hom_img.fwd_model.meas_select = meas_sel;
hom_img.fwd_solve.get_all_meas = 1;

%FWD SOLVE
hom_data = fwd_solve(hom_img);

%Create inverse model same as forward
inv_model = eidors_obj('inv_model', 'EIT_inverse');
inv_model.fwd_model = mdl.fwd_model;

%Solver parameters
inv_model.reconst_type = 'difference'; %difference imaging
inv_model.jacobian_bkgnd.value = 1;
inv_model.solve = @inv_solve_diff_GN_one_step;
inv_model.hyperparameter.value = 3e-3;
inv_model.RtR_prior = @prior_tikhonov;


%% Test
    test_img = hom_img;
    t_circ = makeCircle();
    circle = @(x,y,z) (x-t_circ.x).^2 + (y-t_circ.y).^2 < t_circ.r.^2;

    test_img.elem_data = hom_img.elem_data + elem_select(test_img.fwd_model,circle);
    test_data = fwd_solve(test_img);
    
    set(gca, 'visible', 'off');
    subplot(121);show_fem(test_img); title('Simulated Inhomogeneity');
    set(gca, 'visible', 'off');
    subplot(122);show_fem(inv_solve(inv_model, hom_data, test_data), [0,0,0]); title('Reconstructed Image')

    
    
%% =============== View the Stimulation pattern and save to GIF============

h = figure;
axis tight manual
set(gca,'visible','off');
filename = 'stim_pattern.gif';
vh = fwd_solve(hom_img);
img_v = rmfield(hom_img, 'elem_data');
for n = 1:16
    
    img_v.node_data = vh.volt(:,n);
    show_fem(img_v, [0, 1, 0]);    
    drawnow 
    %UNCOMMENT to create and save GIF
%       % Capture the plot as an image 
%       frame = getframe(h); 
%       im = frame2im(frame); 
%       [imind,cm] = rgb2ind(im,256); 
%       % Write to the GIF File 
%       if n == 1 
%           imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
%       else 
%           imwrite(imind,cm,filename,'gif','WriteMode','append'); 
%       end 
end
  %%
img_v.node_data = vh.volt(:,1);
show_fem(img_v, [0, 1, 0]);  
%% ========================================================================
%===================== Generate training Data =============================
%==========================================================================

Simulation = struct('input', {}, 'output', {});

num_trials = 100;
p = ProgressBar(num_trials);

for  i = 1:num_trials
    
% Temp inverse model
inh_trial = hom_img;

    sel = randi([1, 3]);
    % Randomly select to model a [circle/triangle/square]
    switch(1)
         case 1
             % CIRCLE OBJECT
             t_circ = makeCircle();
             circle = @(x,y,z) (x-t_circ.x).^2 + (y-t_circ.y).^2 < t_circ.r.^2;
             inh_trial.elem_data = hom_img.elem_data + elem_select(inh_trial.fwd_model, circle);                                  
         case 2         
             t_tri = makeTriangle2('t_tri');
             triangle = eval(t_tri.fcn);
             inh_trial.elem_data = hom_img.elem_data + elem_select(inh_trial.fwd_model, triangle);
         case 3
             % SQUARE OBJECT
             t_square = makeSquare();
             square = @(x,y,z) ( y<= (x*t_square.l1(1)+t_square.l1(2)) ... 
                    & y >=(x*t_square.l2(1)+t_square.l2(2)) ...
                    & y >=(x*t_square.l3(1)+t_square.l3(2)) ... 
                    & y <=(x*t_square.l4(1)+t_square.l4(2)) );
             inh_trial.elem_data = hom_img.elem_data + elem_select(inh_trial.fwd_model, square);
    end


    % Fwd solve inv model
    inh_data = fwd_solve(inh_trial);        
    % Inv solve to get image
    img=inv_solve(inv_model, hom_data, inh_data);
        
    % Add noise
    noise = 0.1*std(inh_data.meas - hom_data.meas )*randn(size(hom_data.meas, 1),1);
    inh_data.meas = inh_data.meas + noise;

%     h1=subplot(131); show_fem(inh_trial);
%     h2=subplot(132); show_fem(img);    
%     h3=subplot(133); show_fem(inv_solve(inv_model, hom_data, inh_data));
                  
    %Store the inputs and outputs of the network
    % Input     = differential boundary voltages (208 or 256)
    % Output    = the element data in img
    Simulation(i).input = hom_data.meas - inh_data.meas;
%     Simulation(i).output = img.elem_data;        
    Simulation(i).output = inh_trial.elem_data;        
    
    p.progress;
end
p.stop;


nn_input = zeros(size(hom_data.meas,1), num_trials);
nn_output = zeros(size(img.elem_data,1), num_trials);


for j=1:num_trials
   nn_input(:,j) = Simulation(j).input;
   nn_output(:,j) = Simulation(j).output;
end

%%
% Write to CSV file

write_path = '/Users/matteokalogirou/Google Drive/Colab Notebooks/data/';
fname= [num2str(num_trials/1000) 'k_' ...    %NUMBER OF TRIALS
        num2str(n_elec) ...                  %NUMBER ELECTRODES
        '_c_' ...                   %CIRCLE
        'nonoise_'];                %TRAINED TO RECONSTRUCT NOISEY OUTPUT?

dlmwrite([write_path fname 'input_data.csv'], ...
          nn_input', 'delimiter', ',', '-append');

dlmwrite([write_path fname 'output_data.csv'], ...
        nn_output', 'delimiter', ',', '-append');




