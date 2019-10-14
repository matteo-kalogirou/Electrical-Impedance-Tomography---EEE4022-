%
%   Network Performance Metrics - Simulated Data
%   
%   Simulate a FEM of each shape
%   Reconstruct the simualtion with classical methods
%   For each network configuration, reconstruct the image
%   
%   Measure the MSE between the outputs and the pure element data.
%   Save a picture of each reconstruction and perform PSNR measurement
%   Correlation of the images


run /Users/matteokalogirou/Documents/MATLAB/eidors-v3.9.1-ng/eidors/startup.m

%% --- Define the Models

%--------------------------32 ELECTRODE SYSTEM-----------------------------
I = 0.001;
n_meas = 256;
load('sp_mp.mat');      % Custom stimulation and measurement pattern for tank

stim = stim_meas_list(sp_mp, 32, I, 1);

mdl_32 = mk_common_model('b2c', 32);
z_model = 1;
hom_img_32 = mk_image(mdl_32, z_model);

hom_img_32.fwd_model.stimulation = stim;
hom_img_32.fwd_model.meas_select = true(n_meas,1);
hom_img_32.fwd_solve.get_all_meas = 1;

% FWD SOLVE
hom_data_32 = fwd_solve(hom_img_32);

%Create inverse model same as forward
inv_model_32 = eidors_obj('inv_model', 'EIT_inverse');
f_mdl = mk_common_model('b2c', 16);
inv_model_32.fwd_model = f_mdl.fwd_model;

%Solver parameters
inv_model_32.reconst_type = 'difference'; %difference imaging
inv_model_32.jacobian_bkgnd.value = 1;
inv_model_32.solve = @inv_solve_diff_GN_one_step;
inv_model_32.hyperparameter.value = 3e-3;
inv_model_32.RtR_prior = @prior_tikhonov;

%%
%--------------------------16 ELECTRODE SYSTEM-----------------------------

mdl_16 = mk_common_model('b2c', 16);
z_model = 1;
hom_img = mk_image(mdl_16,z_model);

options = {'meas_current', 'no_rotate_meas'};
[stim, meas_sel] = mk_stim_patterns(16, 1, '{ad}', '{ad}', options, I);

hom_img_16.fwd_model.stimulation = stim;
hom_img_16.fwd_model.meas_select = meas_sel;
hom_img_16.fwd_solve.get_all_meas = 1;

%FWD SOLVE
hom_data_16 = fwd_solve(hom_img_16);

%Create inverse model same as forward
inv_model_16 = eidors_obj('inv_model', 'EIT_inverse');
inv_model_16.fwd_model = mdl_16.fwd_model;

%Solver parameters
inv_model_16.reconst_type = 'difference'; %difference imaging
inv_model_16.jacobian_bkgnd.value = 1;
inv_model_16.solve = @inv_solve_diff_GN_one_step;
inv_model_16.hyperparameter.value = 3e-3;
inv_model_16.RtR_prior = @prior_tikhonov;


%% --- Simulate Inhomogeneity

images_32.c = hom_img_32;
images_32.s = hom_img_32;
images_32.t = hom_img_32;

images_16.c = hom_img_16;
images_16.s = hom_img_16;
images_16.t = hom_img_16;

% CIRCLE
t_circ = makeCircle();
circle = @(x,y,z) (x-t_circ.x).^2 + (y-t_circ.y).^2 < t_circ.r.^2;
images_32.c.elem_data = hom_img_32.elem_data + elem_select(images_32.c.fwd_model, circle);
images_16.c.elem_data = hom_img_16.elem_data + elem_select(images_16.c.fwd_model, circle);
%TRIANGLE
t_tri = makeTriangle2('t_tri');
triangle = eval(t_tri.fcn);
images_32.t.elem_data = hom_img_32.elem_data + elem_select(images_32.t.fwd_model, triangle);
images_26.t.elem_data = hom_img_16.elem_data + elem_select(images_26.t.fwd_model, triangle);
% SQUARE
t_square = makeSquare();
square = @(x,y,z) ( y<= (x*t_square.l1(1)+t_square.l1(2)) ... 
    & y >=(x*t_square.l2(1)+t_square.l2(2)) ...
    & y >=(x*t_square.l3(1)+t_square.l3(2)) ... 
    & y <=(x*t_square.l4(1)+t_square.l4(2)) );
images_32.s.elem_data = hom_img_32.elem_data + elem_select(images_32.s.fwd_model, square);
images_16.s.elem_data = hom_img_16.elem_data + elem_select(images_16.s.fwd_model, square);

%% --- Solve The Reconstructions

% --- Classical

images_32.c.inh_data = fwd_solve(images_32.c);
images_32.c.classical_image = inv_solve(inv_model_32, hom_data_32, images_32.c.inh_data);
noise = 0.1*std(images_32.c.inh_data-hom_data_32.meas)*rand(size(hom_data_32.meas,1),1);
images_32.c.inh_data.meas = images_32.c.inh_data.meas+noise;
% wn - with noise
images_32.c.classical_image_wn = inv_solve(inv_model_32, hom_data_32, images_32.c.inh_data);

images_32.s.inh_data = fwd_solve(images_32.s);
images_32.s.classical_image = inv_solve(inv_model_32, hom_data_32, images_32.s.inh_data);
noise = 0.1*std(images_32.s.inh_data-hom_data_32.meas)*rand(size(hom_data_32.meas,1),1);
images_32.s.inh_data.meas = images_32.s.inh_data.meas+noise;
% wn - with noise
images_32.s.classical_image_wn = inv_solve(inv_model_32, hom_data_32, images_32.s.inh_data);

images_32.t.inh_data = fwd_solve(images_32.t);
images_32.t.classical_image = inv_solve(inv_model_32, hom_data_32, images_32.t.inh_data);
noise = 0.1*std(images_32.t.inh_data-hom_data_32.meas)*rand(size(hom_data_32.meas,1),1);
images_32.t.inh_data.meas = images_32.t.inh_data.meas+noise;
% wn - with noise
images_32.t.classical_image_wn = inv_solve(inv_model_32, hom_data_32, images_32.t.inh_data);


images_16.c.inh_data = fwd_solve(images_16.c);
images_16.c.classical_image = inv_solve(inv_model_16, hom_data_16, images_16.c.inh_data);
noise = 0.1*std(images_16.c.inh_data-hom_data_16.meas)*rand(size(hom_data_16.meas,1),1);
images_16.c.inh_data.meas = images_16.c.inh_data.meas+noise;
% wn - with noise
images_16.c.classical_image_wn = inv_solve(inv_model_16, hom_data_16, images_16.c.inh_data);

images_16.s.inh_data = fwd_solve(images_16.s);
images_16.s.classical_image = inv_solve(inv_model_16, hom_data_16, images_16.s.inh_data);
noise = 0.1*std(images_16.s.inh_data-hom_data_16.meas)*rand(size(hom_data_16.meas,1),1);
images_16.s.inh_data.meas = images_16.s.inh_data.meas+noise;
% wn - with noise
images_16.s.classical_image_wn = inv_solve(inv_model_16, hom_data_16, images_16.s.inh_data);

images_16.t.inh_data = fwd_solve(images_16.t);
images_16.t.classical_image = inv_solve(inv_model_16, hom_data_16, images_16.t.inh_data);
noise = 0.1*std(images_16.t.inh_data-hom_data_16.meas)*rand(size(hom_data_16.meas,1),1);
images_16.t.inh_data.meas = images_16.t.inh_data.meas+noise;
% wn - with noise
images_16.t.classical_image_wn = inv_solve(inv_model_16, hom_data_16, images_16.t.inh_data);

%% --- Network Reconstruction


