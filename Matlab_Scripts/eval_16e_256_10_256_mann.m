%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Testing neural networks on the reconstruction process
%   
%   16 Electrode model with adjacent current drive and measurement
%   (256-10-256)x256
%   Trained on 50k pure Element data
%
%   Matteo Kalogirou
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Run Eidors
run /Users/matteokalogirou/Documents/MATLAB/eidors-v3.9.1-ng/eidors/startup.m

%%
path = '/Users/matteokalogirou/Documents/MATLAB/Electrical Impedance Tomography - EEE4022/Networks/b2c_16_256_';

b2c_pureElem = 'pureElem';
b2c_noiseyElem = 'noiseyElem';


output_density= 256;

% Load the network
for i=1:output_density
    output_net{i} = load([path b2c_pureElem '_mann/b2c_256_' b2c_pureElem '_' num2str(i)  '.mat']);
end


%% Create an image

n_electrode = 16;
n_rings = 1;
I = 0.001;

model = mk_common_model('b2c', n_electrode);

% --- Homogenous Image
z_model = 1;
hom_img = mk_image(model, z_model);
options = {'meas_current', 'no_rotate_meas'};
[stim, meas_sel] = mk_stim_patterns(n_electrode, n_rings, '{ad}', ...
                                        '{ad}', options, I);
hom_img.fwd_model.stimulation = stim;
hom_img.fwd_model.meas_sel = meas_sel;
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

network_input = hom_data.meas - inh_data.meas;

%
% Reconsruct using the network
% nonoise_network_img = hom_img;
network_img = hom_img;

t_start = tic;
for j=1:output_density   
    network_img.elem_data(j) = output_net{j}.a(network_input); %noise
end
time_network_recon = toc(t_start);

% PLotting %
figure(1);
subplot(2,3,1);
[orig_x, orig_y] = makeCirclePlot(0,0,1);
switch(sel)    % Cirlce    
    case 1        
        [phant_x, phant_y] = makeCirclePlot(t_circ.x,t_circ.y,t_circ.r);
        plot(orig_x, orig_y, '-k', phant_x, phant_y, '-r'); axis equal; axis([-1 1 -1 1]);
        set(gca, 'xtick', []);set(gca, 'ytick', []);
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
subplot(2,3,3); show_fem(inh_trial); axis off; title('Pure Element Data');
subplot(2,3,4); show_fem(classical_img); axis off; title('Classical Reconstruction');
subplot(2,3,5); show_fem(classical_no_noise_img); axis off; title('Classical Reconstruction without noise');
subplot(2,3,6); show_fem(network_img); axis off; title('MANN Reconstruction (256-10-256)x256');

%% --- Reconsructing the experimental images

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

% --- Classical reconstruction
inh_left_cimg = inv_solve(inv_model, hom_exp_data_1k, inh_exp_left);
inh_right_cimg = inv_solve(inv_model, hom_exp_data_1k, inh_exp_right);
inh_top_cimg = inv_solve(inv_model, hom_exp_data_1k, inh_exp_top);
inh_bottom_cimg = inv_solve(inv_model, hom_exp_data_1k, inh_exp_bottom);

% --- Network Reconstruction
inh_left_nimg = hom_img;
inh_right_nimg = hom_img;
inh_top_nimg = hom_img;
inh_bottom_nimg = hom_img;


for j=1:output_density   
    inh_left_nimg.elem_data(j) = output_net{j}.a(hom_exp_data_1k + inh_exp_left);
    inh_right_nimg.elem_data(j) = output_net{j}.a(hom_exp_data_1k + inh_exp_right);
    inh_top_nimg.elem_data(j) = output_net{j}.a(hom_exp_data_1k + inh_exp_top);
    inh_bottom_nimg.elem_data(j) = output_net{j}.a(hom_exp_data_1k + inh_exp_bottom);
end

%% --- Plotting 1kHz
figure();
subplot(2,4,1); show_fem(inh_left_cimg); title('Classical Left');
subplot(2,4,2); show_fem(inh_right_cimg); title('Classical Right');
subplot(2,4,3); show_fem(inh_top_cimg); title('Classical Top');
subplot(2,4,4); show_fem(inh_bottom_cimg); title('Classical Bottom');

subplot(2,4,5); show_fem(inh_left_nimg); title('Network Left');
subplot(2,4,6); show_fem(inh_right_nimg); title('Network Right');
subplot(2,4,7); show_fem(inh_top_nimg); title('Network Top');
subplot(2,4,8); show_fem(inh_bottom_nimg); title('Network Bottom');

%% --- 300Hz

% Load inhomogenous data
inh_exp_c1 = importEITinput([p2 '/inh_test_cylinder/']);
inh_exp_c2 = importEITinput([p2 '/inh_test_cylinder2/']);
inh_exp_2c1 = importEITinput([p2 '/inh_test_2cylinders1/']);
inh_exp_2c2 = importEITinput([p2 '/inh_test_2cylinders2/']);
inh_exp_gb = importEITinput([p2 '/inh_test_glassbottle/']);

% --- Classical reconstruction
inh_c1_cimg = inv_solve(inv_model, hom_exp_data_300, inh_exp_c1);
inh_c2_cimg = inv_solve(inv_model, hom_exp_data_300, inh_exp_c2);
inh_2c1_cimg = inv_solve(inv_model, hom_exp_data_300, inh_exp_2c1);
inh_2c2_cimg = inv_solve(inv_model, hom_exp_data_300, inh_exp_2c2);
inh_gb_cimg = inv_solve(inv_model, hom_exp_data_300, inh_exp_gb);

% --- Network Reconstruction
inh_c1_nimg = hom_img;
inh_c2_nimg = hom_img;
inh_2c1_nimg = hom_img;
inh_2c2_nimg = hom_img;
inh_gb_nimg = hom_img;

for j=1:output_density   
    inh_c1_nimg.elem_data(j) = output_net{j}.a(-1*inh_exp_c1);
    inh_c2_nimg.elem_data(j) = output_net{j}.a(-1*inh_exp_c2);
    inh_2c1_nimg.elem_data(j) = output_net{j}.a(-1*inh_exp_2c1);
    inh_2c2_nimg.elem_data(j) = output_net{j}.a(-1*inh_exp_2c2);
    inh_gb_nimg.elem_data(j) = output_net{j}.a(-1*inh_exp_gb);
%     inh_c1_nimg.elem_data(j) = output_net{j}.a(hom_exp_data - inh_exp_c1);
%     inh_c2_nimg.elem_data(j) = output_net{j}.a(hom_exp_data - inh_exp_c2);
%     inh_2c1_nimg.elem_data(j) = output_net{j}.a(hom_exp_data - inh_exp_2c1);
%     inh_2c2_nimg.elem_data(j) = output_net{j}.a(hom_exp_data - inh_exp_2c2);
%     inh_gb_nimg.elem_data(j) = output_net{j}.a(hom_exp_data - inh_exp_gb);
end

%% --- Plotting 300Hz
figure();
subplot(2,5,1); show_fem(inh_c1_cimg); title('Classical Cylinder');
subplot(2,5,2); show_fem(inh_c2_cimg); title('Classical Cylinder');
subplot(2,5,3); show_fem(inh_2c1_cimg); title('Classical 2 Cylinders');
subplot(2,5,4); show_fem(inh_2c2_cimg); title('Classical 2 Cylinders');
subplot(2,5,5); show_fem(inh_gb_cimg); title('Classical Glass Bottle');

subplot(2,5,6); show_fem(inh_c1_nimg); title('Network Cylinder');
subplot(2,5,7); show_fem(inh_c2_nimg); title('Network Cylinder');
subplot(2,5,8); show_fem(inh_2c1_nimg); title('Network 2 Cylinders');
subplot(2,5,9); show_fem(inh_2c2_nimg); title('Network 2 Cylinders');
subplot(2,5,10); show_fem(inh_gb_nimg); title('Network Glass Bottle');

