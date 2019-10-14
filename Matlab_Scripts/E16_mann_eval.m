%
% Evaluate the MANNs
%


clear;
clc;
%%

x = csvread('/Users/matteokalogirou/Google Drive/Colab Notebooks/data/50k_16_c_pureElem_input_data.csv');
t = csvread('/Users/matteokalogirou/Google Drive/Colab Notebooks/data/50k_16_c_pureElem_output_data.csv');

x = x(10000:17499, :);
t = t(10000:17499, :);

%%  b2c_16_256 pureElem

net_16_path = './Networks/b2c_16_256_pureElem_mann/';
output_density = 256;

for i=1:output_density
    E16_mann{i} = load([net_16_path 'b2c_256_pureELem_' num2str(i)  '.mat']);
end


%% Collect performance for each network ?? hUH

N = 256;
T  = 5;
p = ProgressBar(N);

for j = 1:N
    y = zeros(T, 256);
    for i=1:T  
        E16_mann{j}.y(i, 1:256) = E16_mann{j}.a(x(i, :));
    end
    perf(j) = perform(E16_mann{j}.a, t(1:T,:), E16_mann{j}.y);
    p.progress;
end
p.stop;

av_perf = mean(perf);

%% CREATE E16

n_elec = 16;
n_rings = 1;
I = 0.001;
n_meas = 256;

mdl = mk_common_model('b2c', 16);
z_model = 1;
hom_img = mk_image(mdl,z_model);

options = {'meas_current', 'no_rotate_meas'};
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


%%

p = ProgressBar(50);
for K = 1:size(t_circ.x,2)

inh_trial = hom_img;

sel = randi([1, 3]);
sel =1;
switch(sel)
     case 1
         % CIRCLE OBJECT         
%          t_circ = makeCircle();
%              % FOR TESTING
%                t_circ.x = -0.16; t_circ.y = .28; t_circ.r = 0.21;
               circle = @(x,y,z) (x-t_circ.x(K)).^2 + (y-t_circ.y(K)).^2 < t_circ.r(K).^2;
%          circle = @(x,y,z) (x-t_circ.x).^2 + (y-t_circ.y).^2 < t_circ.r.^2;       
         inh_trial.elem_data = hom_img.elem_data + elem_select(inh_trial.fwd_model, circle);                                  
     case 2         
         % TRIANGULAR OBJECT
         t_tri = makeTriangle2('t_tri');
         triangle = eval(t_tri.fcn);
             % FOR TESTING
             inh_trial.elem_data([45, 75, 113, 59, 93, 44, 74, 58, 43 ])= 2;
%          inh_trial.elem_data = hom_img.elem_data + elem_select(inh_trial.fwd_model, triangle);
     case 3
         % SQUARE OBJECT
         t_square = makeSquare();
         square = @(x,y,z) ( y<= (x*t_square.l1(1)+t_square.l1(2)) ... 
                & y >=(x*t_square.l2(1)+t_square.l2(2)) ...
                & y >=(x*t_square.l3(1)+t_square.l3(2)) ... 
                & y <=(x*t_square.l4(1)+t_square.l4(2)) );
         inh_trial.elem_data = hom_img.elem_data + elem_select(inh_trial.fwd_model, square);
    case 4
        % Two Circles
         t_circ1.x = -0.3; t_circ1.y = .3; t_circ1.r = 0.21;
         circle1 = @(x,y,z) (x-t_circ1.x).^2 + (y-t_circ1.y).^2 < t_circ1.r.^2;
         t_circ2.x = 0.3; t_circ2.y = -.3; t_circ2.r = 0.21;
         circle2 = @(x,y,z) (x-t_circ2.x).^2 + (y-t_circ2.y).^2 < t_circ2.r.^2;
         inh_trial.elem_data = hom_img.elem_data + elem_select(inh_trial.fwd_model, circle1);
         inh_trial.elem_data = inh_trial.elem_data + elem_select(inh_trial.fwd_model, circle2);
end

% Fwd solve inv model
inh_data = fwd_solve(inh_trial);        
gn_clean_img = inv_solve(inv_model, hom_data, inh_data);

% SNR = 20db
% noise = 0.1*std( hom_data.meas - inh_data.meas )*randn(size(hom_data.meas, 1),1);
% SNR = 5db
noise = 10^(-5/20)*std( hom_data.meas - inh_data.meas )*randn(size(hom_data.meas, 1),1);
inh_data.meas = inh_data.meas + noise;

% Inv solve model
t_start = tic;
gn_img=inv_solve(inv_model, hom_data, inh_data);
E16GN.time_gn_recon(K) = toc(t_start);

E16GN.mse(K)            = immse(inh_trial.elem_data,gn_img.elem_data);
E16GN.imgCorrCoeff(K)   = corr2(inh_trial.elem_data,gn_img.elem_data); 
E16GN.spread(K)         = spread(gn_img, inh_trial, 0.4, 0.8);


Simulation.input = hom_data.meas - inh_data.meas;
Simulation.output = gn_img.elem_data;

% --- Reconsruct using the network
mann_img = hom_img;

t_start = tic;
for j=1:output_density   
    mann_img.elem_data(j) = E16_mann{j}.a(Simulation.input);   %nonosie
end
E16MANN.netReconstTime(K) = toc(t_start);


E16MANN.mse(K)            = immse(inh_trial.elem_data,mann_img.elem_data);
E16MANN.imgCorrCoeff(K)   = corr2(inh_trial.elem_data,mann_img.elem_data); 
E16MANN.spread(K)         = spread(mann_img, inh_trial,0.8, 0.8);

p.progress;
end
p.stop;
%% --- PLotting

figure(1);
subplot(1,3,1); set(gca, 'FontSize', 14);
show_fem(inh_trial);
[orig_x, orig_y] = makeCirclePlot(0,0,1);
hold on;
switch(sel)    % Cirlce    
    case 1        
        [phant_x, phant_y] = makeCirclePlot(t_circ.x,t_circ.y,t_circ.r);        
        plot(orig_x, orig_y, '-k', phant_x, phant_y, '-g', 'LineWidth', 2); axis equal; axis([-1 1 -1 1]);        
        axis off;
        title('Simulated Phantom - Circle');
    case 2
        
%         plot(orig_x, orig_y, '-k', ...
%              t_tri.intX_l1(1,:), t_tri.intX_l1(2,:), '-g', ...
%              t_tri.intX_l2(1,:), t_tri.intX_l2(2,:), '-g', ...
%              t_tri.intX_l3(1,:), t_tri.intX_l3(2,:), '-g');
%              axis equal; axis([-1 1 -1 1]);
             axis off;
             title('Simulated Phantom - Triangle');
    case 3        
        plot(orig_x, orig_y, '-k', ...
             t_square.intX_l1(1,:), t_square.intX_l1(2,:), '-g', ...
             t_square.intX_l2(1,:), t_square.intX_l2(2,:), '-g', ...
             t_square.intX_l3(1,:), t_square.intX_l3(2,:), '-g', ...
             t_square.intX_l3(1,:), t_square.intX_l3(2,:), '-g', ...
             t_square.intX_l4(1,:), t_square.intX_l4(2,:), '-g');
             axis equal; axis([-1 1 -1 1]);
             axis off;
             title('Simulated Phantom - Square');        
    case 4
        [phant_x, phant_y] = makeCirclePlot(t_circ1.x,t_circ1.y,t_circ1.r);        
        plot(orig_x, orig_y, '-k', phant_x, phant_y, '-g', 'LineWidth', 2); axis equal;
        [phant_x, phant_y] = makeCirclePlot(t_circ2.x,t_circ2.y,t_circ2.r);        
        plot(orig_x, orig_y, '-k', phant_x, phant_y, '-g', 'LineWidth', 2); axis equal;
        axis off;
        title('Simulated Phantom - Circles');
end
hold off;
% subplot(1,3,2); show_fem(gn_clean_img); axis off; title('GN one step (no noise)');
subplot(1,3,2); show_fem(gn_img); axis off; title('GN one step');
set(gca, 'FontSize', 14);
subplot(1,3,3); show_fem(mann_img); axis off; title('MANN');
set(gca, 'FontSize', 14);

figure(2)
% subplot(2,4,5:8);
plot(inh_trial.elem_data-mann_img.elem_data, '-r'); axis tight;
set(gca, 'FontSize', 14);
xlabel('Measurement number');
ylabel('Network Error [S/m]');


%% Performance Metrics

E16MANN.mse            = immse(inh_trial.elem_data,mann_img.elem_data);
E16MANN.imgCorrCoeff   = corr2(inh_trial.elem_data,mann_img.elem_data); 
E16MANN.spread         = spread(mann_img.elem_data, inh_trial.elem_data,...
                                0.8*max(inh_trial.elem_data));

resMSE = [resMSE E16MANN.mse];
resR = [resR E16MANN.imgCorrCoeff];
resSpread = [resSpread E16MANN.spread];

%% Performance Metrics


E16MANN.mse_ = mean(E16MANN.mse);
E16MANN.netReconstTime_ = mean(E16MANN.netReconstTime);
E16MANN.imgCorrCoeff_ = mean(E16MANN.imgCorrCoeff);
E16MANN.spread_ = mean(E16MANN.spread);

E16GN.mse_ = mean(E16GN.mse);
E16GN.time_gn_recon_ = mean(E16GN.time_gn_recon);
E16GN.imgCorrCoeff_ = mean(E16GN.imgCorrCoeff);
E16GN.spread_ = mean(E16GN.spread);




%% ========================================================================
%              --- RECONSTRUCT NEW EXPERIMENTAL DATA ---
%==========================================================================

% Load Experimental Data
Hom_path = './NewResults/hom_1k_1vpp_';
b6_path = './NewResults/banana_6oclk_1k_1vpp/';
b12_path = './NewResults/banana_12oclk_1k_1vpp/';

% Load Hom data
rmsAmp = [];
for j = 1:3
    for i=0:15
        comma2point_overwrite([Hom_path num2str(j) '/frame' num2str(i) '.txt']);
        Frames(:,:,i+1) =  dlmread([Hom_path num2str(j) '/frame' num2str(i) '.txt'], ';', 0,0);
        RmsFrames(:,:,i+1) = rms(Frames(:,:,i+1), 1);
    end
    x = reshape(RmsFrames, [1,256]);
    rmsAmp = [rmsAmp ; x];
end

Hom_data = mean(rmsAmp)';

% Load Inh data
b6 = importEITinput(b6_path);
b12 = importEITinput(b12_path);

% Reconstruct Image == 3e-4
inv_model.hyperparameter.value = 3e-4; 

% Classical Reconstr
img_b6_gn = inv_solve(inv_model, Hom_data, b6);
img_b12_gn = inv_solve(inv_model, Hom_data, b12);

img_b6 = hom_img;
img_b12 = hom_img;
for j=1:output_density   
    img_b6.elem_data(j)  = E16_mann{j}.a(Hom_data - b6);   
    img_b12.elem_data(j) = E16_mann{j}.a(Hom_data - b12);   
end
   
figure();
subplot(2,3,1); show_fem(img_b12_gn); axis off; title('Banana at 12');
set(gca, 'FontSize', 14);
subplot(2,3,2); show_fem(img_b12_gn); axis off; title('GN');
set(gca, 'FontSize', 14);
subplot(2,3,3); show_fem(img_b12); axis off; title('MANN');
set(gca, 'FontSize', 14);
subplot(2,3,4); show_fem(img_b6_gn); axis off; title('Banana at 6');
set(gca, 'FontSize', 14);
subplot(2,3,5); show_fem(img_b6_gn); axis off; title('GN');
set(gca, 'FontSize', 14);
subplot(2,3,6); show_fem(img_b6); axis off; title('MANN');
set(gca, 'FontSize', 14);

E16MANN.meas_cor(1) = corr2(img_b6.elem_data,img_b6_gn.elem_data);
E16MANN.meas_cor(2) = corr2(img_b12.elem_data,img_b12_gn.elem_data);


%%

inv_model.hyperparameter.value = 3e-4; % best results 3e-4
img_b6_gn = inv_solve(inv_model, Hom_data, b6_nm);
img_b12_gn = inv_solve(inv_model, Hom_data, b12_nm);

figure();

subplot(1,2,1); show_fem(img_b6_gn); axis off;
subplot(1,2,2); show_fem(img_b12_gn); axis off;

    
    
%% ========================================================================
%              --- RECONSTRUCT NEW EXPERIMENTAL DATA ---
%==========================================================================

Hom_path = './NewResults/no_mux_hom_';
b6_path = './NewResults/no_mux_ban_6oclk/';
b12_path = './NewResults/no_mux_ban_12oclk/';

% Load Hom data
rmsAmp = [];
for j = 1:2
    for i=0:15
        comma2point_overwrite([Hom_path num2str(j) '/frame' num2str(i) '.txt']);
        Frames(:,:,i+1) =  dlmread([Hom_path num2str(j) '/frame' num2str(i) '.txt'], ';', 0,0);
        RmsFrames(:,:,i+1) = rms(Frames(:,:,i+1), 1);
    end
    x = reshape(RmsFrames, [1,256]);
    rmsAmp = [rmsAmp ; x];
end
%%
Hom_data = mean(rmsAmp)';
%SCALEING

% Load Inh data
b6_nm = importEITinput(b6_path);
b12_nm = importEITinput(b12_path);

%  Reconstruct Image

inv_model.hyperparameter.value = 3e-4;

img_b6_gn = inv_solve(inv_model, Hom_data, b6_nm);
img_b12_gn = inv_solve(inv_model, Hom_data, b12_nm);

img_b6_nm = hom_img;
img_b12_nm = hom_img;
for j=1:output_density   
    img_b6_nm.elem_data(j)  = E16_mann{j}.a(Hom_data - b6_nm);   
    img_b12_nm.elem_data(j) = E16_mann{j}.a(Hom_data - b12_nm);   
end
   
figure();
subplot(2,3,1); show_fem(img_b12_gn); axis off; title('Banana at 12');
set(gca, 'FontSize', 14);
subplot(2,3,2); show_fem(img_b12_gn); axis off; title('GN');
set(gca, 'FontSize', 14);
subplot(2,3,3); show_fem(img_b12_nm); axis off; title('MANN');
set(gca, 'FontSize', 14);
subplot(2,3,4); show_fem(img_b6_gn); axis off; title('Banana at 6');
set(gca, 'FontSize', 14);
subplot(2,3,5); show_fem(img_b6_gn); axis off; title('GN');
set(gca, 'FontSize', 14);
subplot(2,3,6); show_fem(img_b6_nm); axis off; title('MANN');
set(gca, 'FontSize', 14);

%%






