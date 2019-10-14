%
%
% Testing to see if I can reconstruct an image from the collected
% experimental data using the reshaping input hack.
%
%

run /Users/matteokalogirou/Documents/MATLAB/eidors-v3.9.1-ng/eidors/startup.m

%%
run import_hom_exp_data.m

%% Choose Inh data to use

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


%%
% --- Modelling

imdl = mk_common_model('b2c', 16);
himg = mk_image(imdl,1);

options = {'meas_current', 'no_rotate_meas'};
[stim, meas_sel] = mk_stim_patterns(16, 1, '{ad}', '{ad}', options, 0.001);

himg.fwd_model.stimulation = stim;
himg.fwd_model.meas_select = meas_sel;

imdl.jacobian_bkgnd.value = 1;
imdl.solve = @inv_solve_diff_GN_one_step;
imdl.hyperparameter.value = 3e-3;
imdl.RtR_prior = @prior_tikhonov;


img_left = inv_solve(imdl, hom_data, inh_data_left);
img_right = inv_solve(imdl, hom_data, inh_data_right);
img_top = inv_solve(imdl, hom_data, inh_data_top);
img_bottom = inv_solve(imdl, hom_data, inh_data_bottom);

img_cylinder = inv_solve(imdl, hom_data, inh_data_cylinder);
img_cylinder2 = inv_solve(imdl, hom_data, inh_data_cylinder2);
img_2cylinder1 = inv_solve(imdl, hom_data, inh_data_2cylinder1);
img_2cylinder2 = inv_solve(imdl, hom_data, inh_data_2cylinder2);
img_glassbottle = inv_solve(imdl, hom_data, inh_data_glassbottle);


%% Plot

figure(1); title('1 KHz');
h11 = subplot(4,1,1); show_fem(img_top, [1, 1, 0]); title('Top');
h12 = subplot(4,1,2); show_fem(img_left, [1, 1, 0]); title('Left');
h13 = subplot(4,1,3); show_fem(img_right, [1, 1, 0]); title('Right');
h14 = subplot(4,1,4); show_fem(img_bottom, [1, 1, 0]); title('Bottom');
img_bottom.calc_colours.cb_shrink_move = [.3,.8,-0.02];
common_colourbar([h11,h12,h13,h14]);
%%

x = linspace(1,256,256);
figure(2);
subplot(4,1,1); plot(x,hom_data, x, inh_data_left); legend('Hom', 'Inh');
title('Left');
subplot(4,1,2); plot(x,hom_data, x, inh_data_right); legend('Hom', 'Inh');
title('Right');
subplot(4,1,3); plot(x,hom_data, x, inh_data_top); legend('Hom', 'Inh');
title('Top');
subplot(4,1,4); plot(x,hom_data, x, inh_data_bottom); legend('Hom', 'Inh');
title('Bottom');

figure(7); title('300Hz');
subplot(3,3,2); show_fem(img_cylinder, [1, 1, 0]); title('Cylinder');
subplot(3,3,4); show_fem(img_cylinder2, [1, 1, 0]); title('Cylinder');
subplot(3,3,8); show_fem(img_2cylinder1, [1, 1, 0]); title('2 Cylinders');
subplot(3,3,6); show_fem(img_2cylinder2, [1, 1, 0]); title('2 Cylinders');
subplot(3,3,8); show_fem(img_glassbottle, [1, 1, 0]); title('Glass Bottle');

figure(4);
subplot(5,1,1); plot(x,hom_data, x, inh_data_cylinder); legend('Hom', 'Inh');
title('Cylinder');
subplot(5,1,2); plot(x,hom_data, x, inh_data_cylinder2); legend('Hom', 'Inh');
title('Cylinder');
subplot(5,1,3); plot(x,hom_data, x, inh_data_2cylinder1); legend('Hom', 'Inh');
title('2 Cylinders');
subplot(5,1,4); plot(x,hom_data, x, inh_data_2cylinder2); legend('Hom', 'Inh');
title('2 Cylinders');
subplot(5,1,5); plot(x,hom_data, x, inh_data_glassbottle); legend('Hom', 'Inh');
title('Glass Bottle');

%% Reconstrut using ANN

