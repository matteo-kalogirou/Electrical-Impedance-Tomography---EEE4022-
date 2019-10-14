%
%   Evaluate the performance of keras ANN networks
%   256-20-256
%   Adam
%

% Load the Keras Network

path = '';
model_name = '';

network = importKerasNetwork([path model_name]);
network = network.model;

%% Simulate an inh in a 32e model

n_elec = 32;
I = 0.001;
n_meas = 256;
load('sp_mp.mat');      % Custom stimulation and measurement pattern for tank

stim = stim_meas_list(sp_mp, n_elec, I, 1);

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
inv_model.reconst_type = 'difference';
inv_model.jacobian_bkgnd.value = 1;
inv_model.solve = @inv_solve_diff_GN_one_step;
inv_model.hyperparameter.value = 3e-3;
inv_model.RtR_prior = @prior_tikhonov;

%% Simulated Data

test_img = hom_img;
t_circ = makeCircle();
circle = @(x,y,z) (x-t_circ.x).^2 + (y-t_circ.y).^2 < t_circ.r.^2;

test_img.elem_data = hom_img.elem_data + elem_select(test_img.fwd_model,circle);
test_data = fwd_solve(test_img);

subplot(131);show_fem(test_img); title('Simulated Inhomogeneity');    
subplot(132);show_fem(inv_solve(inv_model, hom_data, test_data), [0,0,0]); title('Classical Reconstruction')

% Reconstruct the image - network
network_img = hom_img;
network_img.elem_data = network.predict(test_data);

subplot(133); show_fem(network_img); title('Network Reconstruction');

%% Load the experimental data and test the network on it

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

inh_left_nimg.elem_data = network.predict(hom_exp_data_1k + inh_exp_left);
inh_right_nimg.elem_data = network.predict(hom_exp_data_1k + inh_exp_right);
inh_top_nimg.elem_data = network.predict(hom_exp_data_1k + inh_exp_top);
inh_bottom_nimg.elem_data = network.predict(hom_exp_data_1k + inh_exp_bottom);

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
inh_c1_cimg     = inv_solve(inv_model, hom_exp_data_300, inh_exp_c1);
inh_c2_cimg     = inv_solve(inv_model, hom_exp_data_300, inh_exp_c2);
inh_2c1_cimg    = inv_solve(inv_model, hom_exp_data_300, inh_exp_2c1);
inh_2c2_cimg    = inv_solve(inv_model, hom_exp_data_300, inh_exp_2c2);
inh_gb_cimg     = inv_solve(inv_model, hom_exp_data_300, inh_exp_gb);

% --- Network Reconstruction
inh_c1_nimg = hom_img;
inh_c2_nimg = hom_img;
inh_2c1_nimg = hom_img;
inh_2c2_nimg = hom_img;
inh_gb_nimg = hom_img;

inh_c1_nimg.elem_data   = network.predict.a(-1*inh_exp_c1);
inh_c2_nimg.elem_data   = network.predict(-1*inh_exp_c2);
inh_2c1_nimg.elem_data  = network.predict(-1*inh_exp_2c1);
inh_2c2_nimg.elem_data  = network.predict(-1*inh_exp_2c2);
inh_gb_nimg.elem_data   = network.predict(-1*inh_exp_gb);

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






