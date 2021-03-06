%
% Evaluate the E32 MANNs
%


clear;
clc;
%%

x = csvread('/Users/matteokalogirou/Google Drive/Colab Notebooks/data/50k_16_c_pureElem_input_data.csv');
t = csvread('/Users/matteokalogirou/Google Drive/Colab Notebooks/data/50k_16_c_pureElem_output_data.csv');

x = x(10000:17499, :);
t = t(10000:17499, :);

%%  b2c_16_256 pureElem

net_32_path = './Networks/b2c_32_256_pureElem_mann/';
output_density = 256;

for i=1:output_density
    E32_mann{i} = load([net_32_path 'b2c_256_pureELem_' num2str(i)  '.mat']);
end


%% Collect performance for each network ?? hUH

N = 256;
T  = 5;
p = ProgressBar(N);

for j = 1:N
    y = zeros(T, 256);
    for i=1:T  
        E32_mann{j}.y(i, 1:256) = E32_mann{j}.a(x(i, :));
    end
    perf(j) = perform(E32_mann{j}.a, t(1:T,:), E32_mann{j}.y);
    p.progress;
end
p.stop;

av_perf = mean(perf);

%% CREATE E32 Model



n_elec = 32;
I = 0.001;


%---MEASUREMENT PATTERN -> 16x32 array with 1 and -1
meas_sel = zeros(16,32);
for i = 1:16
    meas_sel(i,i*2) = 1;
    if(i == 16)
       meas_sel(i,2) = -1;
    else 
       meas_sel(i,2*i+2) = -1;
    end
end

%---STIMUATION PATTERN
for i=1:16
    s = zeros(32,1);   
    if(i*2-1 < 31)
        s(i*2-1) = I;                        
        s(i*2+1) = -I;
    else
        s(i*2-1) = I;                        
        s(1) = -I;
    end

    stim(i).stimulation = 'A';
    stim(i).stim_pattern = s;
    stim(i).meas_pattern = meas_sel;
end


% --- FWD MODEL

fmdl = mk_circ_tank(8, [], 32);
fmdl.stimulation = stim;
fmdl.solve = @fwd_solve_1st_order;
fmdl.system_mat = @eidors_default;
fmdl.jacobian = @eidors_default;

if(~valid_fwd_model(fmdl))
    display('INVALID FWD MODEL');
end

% --Homogeneous Image

hom_img = mk_image(fmdl, 1);    %constant conductivity ~ 55 S/m for saline
hom_img.fwd_solve.get_all_meas = 1;

hom_data = fwd_solve(hom_img);

% Simulate Inhomogeneous Object
% t_circ = makeCircle(); %create random circle object
% t_circ.x = -0.16; t_circ.y = .28; t_circ.r = 0.21;
% 
% circle = @(x,y,z) (x-t_circ.x).^2 + (y-t_circ.y).^2 < t_circ.r.^2;
% inh_img = hom_img;
% inh_img.elem_data = hom_img.elem_data + elem_select(inh_img.fwd_model,circle);
% 
% inh_data = fwd_solve(inh_img);

% --- INV MODEL

imdl = eidors_obj('inv_model', '32e EIT Inverse model');
imdl.fwd_model = fmdl;
imdl.reconst_type = 'difference';
imdl.jacobian_bkgnd.value = 1;
imdl.solve = @inv_solve_diff_GN_one_step;
imdl.hyperparameter.value = 3e-6;
imdl.RtR_prior = @prior_tikhonov;


%%
for K = 1:size(t_circ.x,2)

inh_trial = hom_img;
    
sel = randi([1, 3]);
sel =1;
switch(sel)
     case 1
         % CIRCLE OBJECT         
%          t_circ = makeCircle();
%              % For Testing
%              t_circ.x = -0.16; t_circ.y = .28; t_circ.r = 0.21;
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
gn_clean_img = inv_solve(imdl, hom_data, inh_data);

% SNR = 20db
% noise = 0.1*std( hom_data.meas - inh_data.meas )*randn(size(hom_data.meas, 1),1);
% SNR = 5db
noise = 10^(-5/20)*std( hom_data.meas - inh_data.meas )*randn(size(hom_data.meas, 1),1);
inh_data.meas = inh_data.meas + noise;


% Inv solve model
t_start = tic;
gn_img=inv_solve(imdl, hom_data, inh_data);
    E32GN.time_gn_recon(K) = toc(t_start);
    E32GN.mse(K)            = immse(inh_trial.elem_data,gn_img.elem_data);
    E32GN.imgCorrCoeff(K)   = corr2(inh_trial.elem_data,gn_img.elem_data); 
    E32GN.spread(K)         = spread(gn_img, inh_trial, 0.4, 0.8);

Simulation.input = hom_data.meas - inh_data.meas;
Simulation.output = gn_img.elem_data;

% --- Reconsruct using the network
mann_img = hom_img;

t_start = tic;
for j=1:output_density   
    mann_img.elem_data(j) = E32_mann{j}.a(Simulation.input);   %nonosie
end
    E32MANN.netReconstTime(K) = toc(t_start);
    E32MANN.mse(K)            = immse(inh_trial.elem_data,mann_img.elem_data);
    E32MANN.imgCorrCoeff(K)   = corr2(inh_trial.elem_data,mann_img.elem_data); 
    E32MANN.spread(K)         = spread(mann_img, inh_trial,0.8, 0.8);

end
%% --- PLotting

figure();
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
subplot(1,3,2); show_fem(gn_clean_img); axis off; title('GN one step');
set(gca, 'FontSize', 14);
subplot(1,3,3); show_fem(mann_img); axis off; title('MANN');
set(gca, 'FontSize', 14);

% figure(4)
% % subplot(2,4,5:8);
% plot(inh_trial.elem_data-mann_img.elem_data, '-r'); axis tight;
% set(gca, 'FontSize', 14);
% xlabel('Measurement number');
% ylabel('Network Error [S/m]');


%% Performance Metrics


E32MANN.mse_ = mean(E32MANN.mse);
E32MANN.netReconstTime_ = mean(E32MANN.netReconstTime);
E32MANN.imgCorrCoeff_ = mean(E32MANN.imgCorrCoeff);
E32MANN.spread_ = mean(E32MANN.spread);

E32GN.mse_ = mean(E32GN.mse);
E32GN.time_gn_recon_ = mean(E32GN.time_gn_recon);
E32GN.imgCorrCoeff_ = mean(E32GN.imgCorrCoeff);
E32GN.spread_ = mean(E32GN.spread);



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

% Reconstruct Image == 3e-7
imdl.hyperparameter.value = 3e-6; 

% Classical Reconstr
img_b6_gn = inv_solve(imdl, Hom_data, b6);
img_b12_gn = inv_solve(imdl, Hom_data, b12);

img_b6 = hom_img;
img_b12 = hom_img;
for j=1:output_density   
    img_b6.elem_data(j)  = E32_mann{j}.a(Hom_data - b6);   
    img_b12.elem_data(j) = E32_mann{j}.a(Hom_data - b12);   
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

E32MANN.meas_cor(1) = corr2(img_b6.elem_data,img_b6_gn.elem_data);
E32MANN.meas_cor(2) = corr2(img_b12.elem_data,img_b12_gn.elem_data);

%%

imdl.hyperparameter.value = 3e-6; % best results 3e-7
img_b6_gn = inv_solve(imdl, Hom_data, b6);
img_b12_gn = inv_solve(imdl, Hom_data, b12);

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

Hom_data = mean(rmsAmp)';
%SCALEING

% Load Inh data
b6_nm = importEITinput(b6_path);
b12_nm = importEITinput(b12_path);

%  Reconstruct Image

imdl.hyperparameter.value = 3e-5;

img_b6_gn = inv_solve(imdl, Hom_data, b6_nm);
img_b12_gn = inv_solve(imdl, Hom_data, b12_nm);

img_b6_nm = hom_img;
img_b12_nm = hom_img;
for j=1:output_density   
    img_b6_nm.elem_data(j)  = E32_mann{j}.a(Hom_data - b6_nm);   
    img_b12_nm.elem_data(j) = E32_mann{j}.a(Hom_data - b12_nm);   
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

E32MANN.meas_cor(3) = corr2(img_b6_gn.elem_data,img_b6_nm.elem_data);
E32MANN.meas_cor(4) = corr2(img_b12_gn.elem_data,img_b12_nm.elem_data);

%%