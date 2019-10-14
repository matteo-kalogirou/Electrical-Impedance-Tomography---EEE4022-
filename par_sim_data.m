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
% 16-electrode, adjacent current injection (1mA) and measurement system.

n_electrode = 16;
n_rings = 1;
I = 0.001;

% c = mesh density
% 2 = 2D
% c/C = point electrode model/combined electrode model
model = mk_common_model('c2c', n_electrode);

% --- Homogenous Image
% Assume the internal impedance value                                       <-- CHECK Assumption
z_model = 1;
hom_img = mk_image(model, z_model);
%     figure(1);
    % 0/1 = on/off  for [colorBar, electrode_numbering, element_numbering]
    show_fem(hom_img, [0, 1, 1]);

% Don't meas current carrying electrodes, Don't rotate meas w stimulation
% {ad} = adjacent drive and measurement
options = {'no_meas_current', 'no_rotate_meas'};
[stim, meas_sel] = mk_stim_patterns(n_electrode, n_rings, '{ad}', ...
                                        '{ad}', options, I);
hom_img.stimulation = stim;
hom_img.meas_sel = meas_sel;
hom_img.fwd_solve.get_all_meas = 1; %(data.volt = all FEM nodes, not CEM)

% FWD SOLVE -> get the boundary voltages of the model
hom_data = fwd_solve(hom_img);

% === Generate training Data

% An end to end scheme is to be implemented. Therefore, the network
% requires the boundary voltages for each of the training sets (208
% inputs from inh_trial.meas)
% The output for the network will be the element_data associated with the
% model, ie the value of each element in the FEM.

%Create inverse model same as forward
inv_model = eidors_obj('inv_model', 'EIT_inverse');
inv_model.fwd_model=model.fwd_model;

% ===================================TODO: find the best solver
%Solver parameters
inv_model.reconst_type = 'difference'; %difference imaging
inv_model.jacobian_bkgnd.value = 1;
inv_model.solve = @inv_solve_diff_GN_one_step;
inv_model.hyperparameter.value = 3e-3;
inv_model.RtR_prior = @prior_tikhonov;


num_trials = 100;
Simulation = struct('input', {}, 'output', {});

for  i = 1:num_trials
    
% Make a new inverse model
inh_trial = hom_img;

sel = randi([1, 3]);
% Randomly select to model a [circle/triangle/square]
    switch(sel)
         case 1
             % CIRCLE OBJECT
    %          display('Circle');
             t_circ = makeCircle();
             circle = @(x,y,z) (x-t_circ.x).^2 + (y-t_circ.y).^2 < t_circ.r.^2;
             inh_trial.elem_data = hom_img.elem_data + elem_select(inh_trial.fwd_model, circle);                                  
         case 2    
    %          display('Triangle');
             t_tri = makeTriangle2('t_tri');
             triangle = eval(t_tri.fcn);
             inh_trial.elem_data = hom_img.elem_data + elem_select(inh_trial.fwd_model, triangle);
         case 3
             % SQUARE OBJECT
    %          display('Square');
             t_square = makeSquare();
             square = @(x,y,z) ( y<= (x*t_square.l1(1)+t_square.l1(2)) ... 
                    & y >=(x*t_square.l2(1)+t_square.l2(2)) ...
                    & y >=(x*t_square.l3(1)+t_square.l3(2)) ... 
                    & y <=(x*t_square.l4(1)+t_square.l4(2)) );
             inh_trial.elem_data = hom_img.elem_data + elem_select(inh_trial.fwd_model, square);
    end

    % Fwd solve inv model
    inh_data = fwd_solve(inh_trial);        
    % Add noise
    noise = 0.1*std( hom_data.meas - inh_data.meas )*randn(size(hom_data.meas, 1),1);
    inh_data.meas = inh_data.meas + noise;
    % Inv solve model
    img=inv_solve(inv_model, hom_data, inh_data);
    
%     show_fem(img, [0,0,0]);
%     hold on;              %triangle
%     plot(xx, t_tri.l1(1)*xx + t_tri.l1(2), '-r');
%     plot(xx, t_tri.l2(1)*xx + t_tri.l2(2), '-r');
%     plot(xx, t_tri.l3(1)*xx + t_tri.l3(2), '-r');
%     hold off;
%     axis([-1 1 -1 1]);    
%     grid on;
    

               
    %Store the inputs and outputs of the network
    % Input     = differential boundary voltages (208)
    % Output    = the element data in img
    Simulation(i).input = hom_data.meas - inh_data.meas;
    Simulation(i).output = img.elem_data;    
end



nn_input = zeros(size(hom_data.meas,1), num_trials);
nn_output = zeros(size(img.elem_data,1), num_trials);


for j=1:num_trials
   nn_input(:,j) = Simulation(j).input;
   nn_output(:,j) = Simulation(j).output;
end


%% Write to CSV file

%  dlmwrite('.csv', nn_input', 'delimiter', ',', '-append');
%  dlmwrite('.csv', nn_output', 'delimiter', ',', '-append');
