
%%
run import_hom_exp_data.m

%%

model = mk_common_model('b2c', 16);

him = mk_image(model,1);

options = {'meas_current', 'no_rotate_meas'};
[stim, meas_sel] = mk_stim_patterns(16, 1, '{ad}', '{ad}', options, 0.001);

him.fwd_model.stimulation = stim;
him.fwd_model.meas_select = meas_sel;

model.jacobian_bkgnd.value = 1;
model.solve = @inv_solve_diff_GN_one_step;
model.hyperparameter.value = 3e-3;
model.RtR_prior = @prior_tikhonov;

rec_img = inv_solve(model, hom_data, inh_data_left);

show_fem(rec_img, [0,1,0]);

%%


n_elec = 16;
I = 0.001;      %Current = 1mA
n_rings = 1;
mdl_2d = mk_common_model('b2c', n_elec);

zhomg = 1;

homg_img = mk_image(mdl_2d, zhomg);   
options = {'meas_current', 'no_rotate_meas'};        
[stim, meas_sel] = mk_stim_patterns(n_elec, n_rings, '{ad}','{ad}', options, I);
homg_img.fwd_model.stimulation = stim;
homg_img.fwd_model.meas_select = meas_sel;

%Solve forward problem using EIDORS solver
homg_img.fwd_solve.get_all_meas = 1;      
homg_data = fwd_solve(homg_img);

inh_img=homg_img;
shape=@(x,y,z) (x-0.2).^2+(y-0.4).^2<0.25.^2; %circular inhomgeneity
inh_img.elem_data = homg_img.elem_data + elem_select(inh_img.fwd_model,shape);
inh_data = fwd_solve(inh_img);

figure(2);
show_fem(inh_img, [0, 0, 0]);

%Create inverse model same as forward
inv_model = eidors_obj('inv_model', 'EIT_inverse');
inv_model.fwd_model=mdl_2d.fwd_model;

%Solver parameters
inv_model.reconst_type = 'difference'; %difference imaging
inv_model.jacobian_bkgnd.value = 1;
inv_model.solve = @inv_solve_diff_GN_one_step;
inv_model.hyperparameter.value = 3e-3;
inv_model.RtR_prior = @prior_tikhonov;

%Solve inverse problem
img=inv_solve(inv_model, homg_data, inh_data);

%Display reconstruction
figure(3)
show_fem(img,[0,0,0]);
