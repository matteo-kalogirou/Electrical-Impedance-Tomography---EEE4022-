%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Testing neural networks on the reconstruction process
%   208-10-x
%
%   Matteo Kalogirou
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Run Eidors
run /Users/matteokalogirou/Documents/MATLAB/eidors-v3.9.1-ng/eidors/startup.m

%%
path = '/Users/matteokalogirou/Google Drive/UCT/2019/EIT/Networks/Matlab_Networks/';
% b2c_path = 'MANN/b2c_mann/';
b2c_nonoise_path = 'MANN/b2c_nonoise_mann/';
c2c_noise_path = 'MANN/c2c_noise_mann/';

% MANN name
b2c_filename = 'b2c_nonoise_';
c2c_filename = 'mann_576_net_';

output_density= 256;

% Load the network
for i=1:output_density
    nonoise_output_net{i} = load([path b2c_nonoise_path 'b2c_nonoise_' num2str(i)  '.mat']);
%     noise_output_net{i} = load([path b2c_path 'mann_net_' num2str(i)  '.mat']);
%       noise_output_net{i} = load([path c2c_noise_path c2c_filename num2str(i) '.mat']);
end


%% Create an image

n_electrode = 16;
n_rings = 1;
I = 0.001;

% c = mesh density
% 2 = 2D
% c/C = point electrode model/combined electrode model
model = mk_common_model('b2c', n_electrode);

% --- Homogenous Image
% Assume the internal impedance value                                       <-- CHECK Assumption
z_model = 1;
hom_img = mk_image(model, z_model);
options = {'no_meas_current', 'no_rotate_meas'};
[stim, meas_sel] = mk_stim_patterns(n_electrode, n_rings, '{ad}', ...
                                        '{ad}', options, I);
hom_img.stimulation = stim;
hom_img.meas_sel = meas_sel;
hom_img.fwd_solve.get_all_meas = 1;
hom_data = fwd_solve(hom_img);

%Create inverse model same as forward
inv_model = eidors_obj('inv_model', 'EIT_inverse');
inv_model.fwd_model=model.fwd_model;

%Solver parameters
inv_model.reconst_type = 'difference';
inv_model.jacobian_bkgnd.value = 1;
inv_model.solve = @inv_solve_diff_GN_one_step;
inv_model.hyperparameter.value = 3e-3;
inv_model.RtR_prior = @prior_tikhonov;

%%
inh_trial = hom_img;

sel = randi([1, 3]);
switch(sel)
     case 1
         % CIRCLE OBJECT         
         t_circ = makeCircle();
         circle = @(x,y,z) (x-t_circ.x).^2 + (y-t_circ.y).^2 < t_circ.r.^2;
         inh_trial.elem_data = hom_img.elem_data + elem_select(inh_trial.fwd_model, circle);                                  
     case 2         
         % TRIANGULAR OBJECT
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
classical_no_noise_img = inv_solve(inv_model, hom_data, inh_data);
noise = 0.1*std( hom_data.meas - inh_data.meas )*randn(size(hom_data.meas, 1),1);
inh_data.meas = inh_data.meas + noise;
% Inv solve model
classical_img=inv_solve(inv_model, hom_data, inh_data);

Simulation.input = hom_data.meas - inh_data.meas;
Simulation.output = classical_img.elem_data;

% Reconsruct using the network
% nonoise_network_img = hom_img;
noise_network_img = hom_img;

t_start = tic;
for j=1:output_density   
    nonoise_network_img.elem_data(j) = nonoise_output_net{j}.a(Simulation.input);   %nonosie
%     noise_network_img.elem_data(j) = noise_output_net{j}.f(Simulation.input); %noise
end
time_network_recon = toc(t_start);

% PLotting

figure(1);
subplot(3,3,2);
[orig_x, orig_y] = makeCirclePlot(0,0,1);
switch(sel)    % Cirlce    
    case 1        
        [phant_x, phant_y] = makeCirclePlot(t_circ.x,t_circ.y,t_circ.r);
        plot(orig_x, orig_y, '-k', phant_x, phant_y, '-r'); axis equal; axis([-1 1 -1 1]);
        set(gca, 'xtick', []);
        set(gca, 'ytick', []);
        title('Simulated Phantom - Circle');
    case 2
        
        plot(orig_x, orig_y, '-k', ...
             t_tri.intX_l1(1,:), t_tri.intX_l1(2,:), '-r', ...
             t_tri.intX_l2(1,:), t_tri.intX_l2(2,:), '-r', ...
             t_tri.intX_l3(1,:), t_tri.intX_l3(2,:), '-r');
             axis equal; axis([-1 1 -1 1]);
             set(gca, 'xtick', []); set(gca, 'ytick', []);
             title('Simulated Phantom - Triangle');
    case 3        
        plot(orig_x, orig_y, '-k', ...
             t_square.intX_l1(1,:), t_square.intX_l1(2,:), '-r', ...
             t_square.intX_l2(1,:), t_square.intX_l2(2,:), '-r', ...
             t_square.intX_l3(1,:), t_square.intX_l3(2,:), '-r', ...
             t_square.intX_l3(1,:), t_square.intX_l3(2,:), '-r', ...
             t_square.intX_l4(1,:), t_square.intX_l4(2,:), '-r');
             axis equal; axis([-1 1 -1 1]);
             set(gca, 'xtick', []); set(gca, 'ytick', []);
             title('Simulated Phantom - Square');        
end
subplot(3,3,3); show_fem(classical_no_noise_img); axis off; title('Classical Reconstruction without noise');
subplot(3,3,1); show_fem(classical_img); axis off; title('Classical Reconstruction');
subplot(3,3,4); show_fem(noise_network_img); axis off; title('MANN Reconstruction (noise)');
% subplot(3,3,6); show_fem(nonoise_network_img); axis off; title('MANN Reconstruction (nonoise)');
subplot(3,3,7:9);
%     plot(classical_no_noise_img.elem_data-nonoise_network_img.elem_data, '-r');
% axis([0 size(nonoise_network_img.elem_data,1) -1.2 1.2])
% title('Difference: Classical Reconstruction w/o nosie and MANN Reconstruciton w/o noise')



