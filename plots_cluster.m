clear
addpath(genpath('/users/PAS1622/nickmoise/.matlab/down'))

%% check 25 beat CV


scale_gj_loc_vec = [8 4 2 1 1/2 1/4 1/8];
scale_chan_loc_vec = [8 4 2 1 1/2 1/4 1/8];




%%% CV for all beat
% load_file = "/fs/scratch/PAS1622/nickmoise/ID_2025/hetg_25beats_c1000/" + ...
%             "msh1_baseline_msh2_p60_ip60_cycle1000_beats25_D01_gj_loc1_8_chan_loc1_8.mat";


%all except G_save - we don't need that for plots
% load(load_file,'phi_save', 'I_save','S_save', 'bcl', 'D', 'scale_gj_loc', 'scale_chan_loc', ...
%         'p', 'iEC', 'Nnodes', 'Ncell', 'Ncurrents', 'indices', 'Mdisc', ...
%         'phi_axial_all', 'Iind', 'ts', 'model', 'FEM_file_list', 'tissue_legend',...
%         'tup', 'trepol', 'ts_save', 'Nint');


data_dir = "/fs/scratch/PAS1622/nickmoise/ID_2025/hetg_25beats_c1000/";
data_list = dir(data_dir);
data_list = natsortfiles(data_list);
data_list(1:2) = [];

thresh_activation = -10;


for i_all = 1:length(data_list)

    fprintf("%d / %d \n", i_all, length(data_list))
    filename = data_dir + "/" +  data_list(i_all).name;
    load(filename, 'phi_axial_all','ts', 'scale_gj_loc', 'scale_chan_loc')




    x = ts;
    cell_index = [10,40]; % cells from which CV is comp
    activ_time = zeros(2,10);
    for i_cell = 1:length(cell_index)
        data_iso_int = phi_axial_all(cell_index(i_cell),:);
        data_iso1 = data_iso_int(:,1:end-1);
        data_iso2 = data_iso_int(:,2:end);
        ind_th = find(data_iso1<thresh_activation & data_iso2> thresh_activation); % find indices below/above threshold
        y1 = data_iso1(ind_th); y2 = data_iso2(ind_th);
        m_slope = (y2-y1)./(x(ind_th+1)-x(ind_th));  % linear slope
        activ_time_cell = x(ind_th) - (data_iso_int(ind_th)-thresh_activation)./m_slope;
        activ_time(i_cell,1:length(activ_time_cell)) = activ_time_cell;
    end
    
    cell_dist = (cell_index(2) - cell_index(1))*100; %100 = cell length in um
    CV = cell_dist./(activ_time(2,:) - activ_time(1,:)) ./10;  

    CV_list(i_all,:) = CV;

    
end

plot(CV_list')

save_dir = "/fs/ess/PAS1622/nickmoise/ID_2025/" + ...
           "processed_data/";

save_name = save_dir + "CV_25b_1000_all";
save(save_name,'CV_list', 'data_list')

%% single sim hetg CV
thresh_activation = -10;

x = ts;
activ_time = [];
for i = 1:50
    data_iso_int = phi_axial_all(i,:);
    data_iso1 = data_iso_int(:,1:end-1);
    data_iso2 = data_iso_int(:,2:end);
    ind_th = find(data_iso1<thresh_activation & data_iso2> thresh_activation); % find indices below/above threshold
    y1 = data_iso1(ind_th); y2 = data_iso2(ind_th);
    m_slope = (y2-y1)./(x(ind_th+1)-x(ind_th));  % linear slope
    for j = 1:length(x(ind_th))
        activ_vec = x(ind_th) - (data_iso_int(ind_th)-thresh_activation)./m_slope;
        activ_time(i,j) = activ_vec(j);
    end
end

activ_time(5:end,3:5) = activ_time(5:end,2:4);
activ_time(5:end,2) = NaN;
% 
% activ_time(11:end,5) = activ_time(11:end,4);
% activ_time(11:end,4) = NaN;

beat_no = 3;
local_cv = (100./diff(activ_time(:,beat_no)))./10; %cm/s
% local_cv = local_cv(5:45); %cm/s


figure
imagesc()
plot(local_cv, 'LineWidth',2)
box off
set(gca,'FontSize',12,'TickDir','out');
xlim([1,50])
xlabel('cell no')
ylabel("CV (cm/s)")
xticks([1,10:10:50])

%% calcium hetg
data_dir = "/fs/scratch/PAS1622/nickmoise/ID_2025/" + ...
           "gj_chan_loc_base_D01/";

save_dir = "/fs/scratch/PAS1622/nickmoise/ID_2025/" + ...
           "processed_data/";

save_name = save_dir + "Ca_hetg_mat_base_D01";


data_list = dir(data_dir);
data_list = natsortfiles(data_list);
data_list(1:2) = [];

scale_gj_loc_vec = [8 4 2 1 1/2 1/4 1/8];
scale_chan_loc_vec = [8 4 2 1 1/2 1/4 1/8];
[X,Y] = meshgrid(scale_gj_loc_vec,scale_chan_loc_vec);

Ca_mat = zeros(7,7);

Njunc = 25;
thresh_activation = 0;
cycle_vec = [200:5:300 350:50:1000];
for i_cycle = 1:length(data_list)
    fprintf("%d / %d \n", i_cycle, length(data_list))
    filename = data_dir + "/" +  data_list(i_cycle).name;
    load(filename, 'S_save', 'iEC', 'Nnodes', 'Mdisc', ...
         'scale_gj_loc', 'scale_chan_loc')

    Ca_cleft_all = S_save(iEC+2*Nnodes,:);
    Ca_cleft = Ca_cleft_all((Njunc-1)*Mdisc+1:Njunc*Mdisc,:); 

    max_Ca = max(Ca_cleft,[],'all');

    Ca_mat(find(scale_gj_loc_vec == scale_gj_loc), ...
       find(scale_chan_loc_vec == scale_chan_loc)) = max_Ca;

end

figure
surf(log2(X),log2(Y),Ca_mat)
xlabel('chan loc')
ylabel('gj loc')

save(save_name,'Ca_mat')


%% CV map scale
% msh1_baseline_cycle1000_beats10_D01_gj_loc0125_0125_chan_loc0125_0125.mat

data_dir = "/fs/scratch/PAS1622/nickmoise/ID_2025/" + ...
           "gj_chan_only_na/";

save_dir = "/fs/scratch/PAS1622/nickmoise/ID_2025/" + ...
           "processed_data/";

sim_root = "msh2_p60_ip60_cycle1000_beats10_D1_gj_loc";
save_name = save_dir + "CV_map7x7_na_" + sim_root + "_chan_loc";


% data_list = dir(data_dir);
% data_list = natsortfiles(data_list);
% data_list(1:2) = [];

scale_gj_loc_vec = [8 4 2 1 1/2 1/4 1/8];
scale_chan_loc_vec = [8 4 2 1 1/2 1/4 1/8];
CV_map = zeros(7,7);
[X,Y] = meshgrid(scale_gj_loc_vec,scale_chan_loc_vec);
cycle = 1000;

thresh_activation = 0;
for i_all = 1:length(scale_gj_loc_vec)
for j_all = 1:length(scale_chan_loc_vec)
    fprintf("%d / %d \n", i_all,j_all)
    file = sim_root ...
     + string(scale_gj_loc_vec(i_all)) + "_" + string(scale_gj_loc_vec(i_all)) ...
     + "_chan_loc" ...
     + string(scale_chan_loc_vec(j_all)) + "_" + string(scale_chan_loc_vec(j_all));
    file = strrep(file, '.', '');

    
    filename = data_dir + file + ".mat";

    if isfile(filename)
        load(filename, 'phi_axial_all','ts', 'scale_gj_loc', 'scale_chan_loc')
        scale_gj_loc = scale_gj_loc(1);
        scale_chan_loc = scale_chan_loc(1);
    else
        disp("missing: " + filename)
        CV_map(i_cycle) = NaN;
        continue
    end
   
   
    phi_axial_all(:,end) = phi_axial_all(:,end-1);
    ind_last_2beat = find(ts>(cycle*8-10));
    phi_axial_2beat = phi_axial_all(:,ind_last_2beat);
    x = ts(ind_last_2beat);
    cell_index = [10,40]; % cells from which CV is comp
    activ_time = zeros(2,20);
    for i_cell = 1:length(cell_index)
        data_iso_int = phi_axial_2beat(cell_index(i_cell),:);
        data_iso1 = data_iso_int(:,1:end-1);
        data_iso2 = data_iso_int(:,2:end);
        ind_th = find(data_iso1<thresh_activation & data_iso2> thresh_activation); % find indices below/above threshold
        y1 = data_iso1(ind_th); y2 = data_iso2(ind_th);
        m_slope = (y2-y1)./(x(ind_th+1)-x(ind_th));  % linear slope
        activ_time_cell = x(ind_th) - (data_iso_int(ind_th)-thresh_activation)./m_slope;
        activ_time(i_cell,1:length(activ_time_cell)) = activ_time_cell;
    end
    
    beats1 = activ_time(1,:) ~= 0;
    beats1(beats1==0) = [];
    beats2 = activ_time(2,:) ~= 0;
    beats2(beats2==0) = [];
    
    cell_dist = (cell_index(2) - cell_index(1))*100; %100 = cell length in um
    CV = cell_dist./(activ_time(2,length(beats2)) - activ_time(1,length(beats1))) ./10;    

    CV_map(find(scale_gj_loc_vec == scale_gj_loc), ...
       find(scale_chan_loc_vec == scale_chan_loc)) = CV;
    
    % figure
    % imagesc(ts,1:50,phi_axial_all)
    % ax = gca;
    % ax.DataAspectRatio = [100 1 1];
    
end
end
% plot(ts,phi_axial_all(5,:))

% plot(CV_cycles)




figure
surf(log2(X),log2(Y),CV_map)
xlabel('chan loc')
ylabel('gj loc')


save(save_name,'CV_map')


%% CV resti batch

%missing files:
% missing: /fs/scratch/PAS1622/nickmoise/ID_2025/homog_tissue_resti/msh1_baseline_cycle225_beats10_D01_gj_loc8_8_chan_loc8_8.mat
% msh2_p60_ip60
data_dir = "/fs/scratch/PAS1622/nickmoise/ID_2025/" + ...
           "homog_tissue_resti/";

save_dir = "/fs/scratch/PAS1622/nickmoise/ID_2025/" + ...
           "processed_data/";

gj_vec = [1];
chan_vec = [1];
comb_vec = combvec(gj_vec, chan_vec);

% sim_root = "msh2_p60_ip60_cycle%d_beats10_D01_gj_loc";
sim_root = "msh1_baseline_cycle%d_beats10_D01_gj_loc";
save_name = save_dir + "CV_resti_" + sprintf(sim_root,[]) + "_chan_loc";

% save_name = strrep(save_name, '.', ''); %remove dot to prevent file extension errors   
% cycle_vec = [200:5:300 350:50:1000];
cycle_vec = 225;
CV_cycles_all = zeros(length(comb_vec),length(cycle_vec));
beats_2_cycles_all = zeros(length(comb_vec),length(cycle_vec));

for i_all = 1:length(comb_vec)

file = sim_root ...
     + string(comb_vec(1,i_all)) + "_" + string(comb_vec(1,i_all)) ...
     + "_chan_loc" ...
     + string(comb_vec(2,i_all)) + "_" + string(comb_vec(2,i_all));
file = strrep(file, '.', '');


data_list = dir(data_dir);
data_list = natsortfiles(data_list);
data_list(1:2) = [];


thresh_activation = 0;


CV_cycles = zeros(1,length(cycle_vec));
beats_2_cycles = zeros(1,length(cycle_vec));
for i_cycle = 1:length(cycle_vec)
    % fprintf("%d / %d \n", i_file, length(data_list))
    % filename = data_dir + "/" +  data_list(i_file).name;

    fprintf("%d / %d \n", i_cycle, length(cycle_vec))
    load_file = sprintf(file, cycle_vec(i_cycle));
    filename = data_dir + load_file + ".mat";
    if isfile(filename)
        load(filename, 'phi_axial_all','ts')
    else
        disp("missing: " + filename)
        CV_cycles(i_cycle) = NaN;
        continue
    end

    phi_axial_all(:,end) = phi_axial_all(:,end-1);
    ind_last_2beat = find(ts>(cycle_vec(i_cycle)*8-10));
    phi_axial_2beat = phi_axial_all(:,ind_last_2beat);
    
    x = ts(ind_last_2beat);
    cell_index = [10,40]; % cells from which CV is comp
    activ_time = zeros(2,20);
    for i_cell = 1:length(cell_index)
        data_iso_int = phi_axial_2beat(cell_index(i_cell),:);
        data_iso1 = data_iso_int(:,1:end-1);
        data_iso2 = data_iso_int(:,2:end);
        ind_th = find(data_iso1<thresh_activation & data_iso2> thresh_activation); % find indices below/above threshold
        y1 = data_iso1(ind_th); y2 = data_iso2(ind_th);
        m_slope = (y2-y1)./(x(ind_th+1)-x(ind_th));  % linear slope
        activ_time_cell = x(ind_th) - (data_iso_int(ind_th)-thresh_activation)./m_slope;
        activ_time(i_cell,1:length(activ_time_cell)) = activ_time_cell;
    end
    
    beats1 = activ_time(1,:) ~= 0;
    beats1(beats1==0) = [];
    beats2 = activ_time(2,:) ~= 0;
    beats2(beats2==0) = [];
    
    cell_dist = (cell_index(2) - cell_index(1))*100; %100 = cell length in um

    
    
    if activ_time(2,2)>activ_time(1,2)
        beats_2_cycles(i_cycle) = 1;
        CV = cell_dist./(activ_time(2,2) - activ_time(1,2)) ./10;
        CV_cycles(i_cycle) = CV;
    else
        beats_2_cycles(i_cycle) = 0;
        CV = cell_dist./(activ_time(2,1) - activ_time(1,1)) ./10;
        CV_cycles(i_cycle) = CV;
    end
   
    
    
    % load(filename, 'scale_gj_loc', 'scale_chan_loc')
    % 
    % CV_mat(find(scale_gj_loc_vec == scale_gj_loc), ...
    %        find(scale_chan_loc_vec == scale_chan_loc)) = CV_cycles(i_file);


    % figure
    % imagesc(ts,1:50,phi_axial_all)
    % ax = gca;
    % ax.DataAspectRatio = [100 1 1];

end

% plot(ts(ind_last_2beat),phi_axial_all(5,:))


% plot(cycle_vec,CV_cycles)
% hold on


% 
% 
x_time = ts;
x_img = 0:x_time(end);
y_img = 1:50;
[X,Y] = meshgrid(x_time, y_img); %sample grid points
[Xq,Yq] = meshgrid(x_img, y_img);%querry grid points
phi_axial_interp = interp2(X,Y,phi_axial_all,Xq,Yq);

f = figure;
imagesc(x_img,y_img,phi_axial_interp)
colormap coolwarm
ax = gca;
ax.DataAspectRatio = [50 1 1];
% title("k_{gj} = " + string(comb_vec(1,i_all)) + ", k_{chan}"+ string(comb_vec(2,i_all)))
box off
set(gca,'FontSize',10,'TickDir','out','YDir','normal');
xlabel("time (ms)")
ylabel("cell no")
clim([-90,45])
% saveas(f,"/fs/scratch//PAS1622/nickmoise/ID_2025/processed_data/figures/" + string(i_all) + ".png")
% xlim([9000,9300])


CV_cycles_all(i_all,:) = CV_cycles;
beats_2_cycles_all(i_all,:) = beats_2_cycles;
end

% save(save_name,'cycle_vec','comb_vec','CV_cycles_all','beats_2_cycles_all')


%% CV resti/multi cycle hetg

%missing files:
% missing: /fs/scratch/PAS1622/nickmoise/ID_2025/homog_tissue_resti/msh1_baseline_cycle225_beats10_D01_gj_loc8_8_chan_loc8_8.mat



% msh1_baseline_cycle200_beats10_D01_gj_loc1_0125_chan_loc1_0125.mat

% msh1_baseline_msh2_p60_ip60_cycle250_beats10_D01_gj_loc1_0125_chan_loc1_8.mat


data_dir = "/fs/scratch/PAS1622/nickmoise/ID_2025/" + ...
           "hetg_tissue_restart_50b_all_clampall/";

save_dir = "/fs/scratch/PAS1622/nickmoise/ID_2025/" + ...
           "processed_data/";

gj_vec = [1 8];
chan_vec = [1 8];
comb_vec = combvec(gj_vec, chan_vec);

mesh_mid_vec = [1,2];
D_vec = [0.1, 1];

D_mesh_c = combvec(mesh_mid_vec,D_vec);

for i_m_D = 1:length(D_mesh_c)

mesh_mid = D_mesh_c(1,i_m_D);
D_str = D_mesh_c(2,i_m_D);

D_str = strrep(string(D_str), '.', '');

if mesh_mid == 1
    sim_root = "msh1_baseline_cycle%d_beats50_D" + D_str + "_gj_loc";
else
    sim_root = "msh1_baseline_msh2_p60_ip60_cycle%d_beats50_D" + D_str + "_gj_loc";
end

save_name = save_dir + "end_cycle_clampall_" + sprintf(sim_root,[]) + "_chan_loc";

% save_name = strrep(save_name, '.', ''); %remove dot to prevent file extension errors   
cycle_vec = [200:5:300];
% cycle_vec = [220];

% CV_cycles_all = zeros(length(comb_vec),length(cycle_vec));
end_cycles = zeros(length(comb_vec),length(cycle_vec));

% for i_all = [1,9]%1:length(comb_vec)
for i_all = 1:length(comb_vec)
file = sim_root ...
     + string(1) + "_" + string(comb_vec(1,i_all)) ...
     + "_chan_loc" ...
     + string(1) + "_" + string(comb_vec(2,i_all));
file = strrep(file, '.', '');

data_list = dir(data_dir);
data_list = natsortfiles(data_list);
data_list(1:2) = [];

thresh_activation = 0;
% CV_cycles_all = zeros(1,length(cycle_vec));
for i_cycle = 1:length(cycle_vec)
    % fprintf("%d / %d \n", i_file, length(data_list))
    % filename = data_dir + "/" +  data_list(i_file).name;

    fprintf("%d / %d \n", i_cycle, length(cycle_vec))
    load_file = sprintf(file, cycle_vec(i_cycle));
    filename = data_dir + load_file + ".mat";
    if isfile(filename)
        load(filename, 'phi_axial_all','ts','scale_chan_loc')
        disp(filename)
    else
        disp("missing: " + filename)
        CV_cycles(i_cycle) = NaN;
        continue
    end

    % plot(ts,phi_axial_all(50,:))

    [~,peak_locs] = findpeaks(phi_axial_all(end,:),'MinPeakHeight',5,'MinPeakDistance',150);

    end_cycles(i_all, i_cycle) = length(peak_locs);


    % ind_last_2beat = find(ts>(cycle_vec(i_file)*8-10));
    % phi_axial_2beat = phi_axial_all(:,ind_last_2beat);
    
    % x = ts(ind_last_2beat);
    % cell_index = [10,40]; % cells from which CV is comp
    % activ_time = zeros(2,20);
    % for i_cell = 1:length(cell_index)
    %     data_iso_int = phi_axial_2beat(cell_index(i_cell),:);
    %     data_iso1 = data_iso_int(:,1:end-1);
    %     data_iso2 = data_iso_int(:,2:end);
    %     ind_th = find(data_iso1<thresh_activation & data_iso2> thresh_activation); % find indices below/above threshold
    %     y1 = data_iso1(ind_th); y2 = data_iso2(ind_th);
    %     m_slope = (y2-y1)./(x(ind_th+1)-x(ind_th));  % linear slope
    %     activ_time_cell = x(ind_th) - (data_iso_int(ind_th)-thresh_activation)./m_slope;
    %     activ_time(i_cell,1:length(activ_time_cell)) = activ_time_cell;
    % end
    % 
    % beats1 = activ_time(1,:) ~= 0;
    % beats1(beats1==0) = [];
    % beats2 = activ_time(2,:) ~= 0;
    % beats2(beats2==0) = [];
    % 
    % cell_dist = (cell_index(2) - cell_index(1))*100; %100 = cell length in um
    % CV = cell_dist./(activ_time(2,length(beats2)) - activ_time(1,length(beats1))) ./10;
    % CV_cycles(i_file) = CV;
    

    % x_time = ts;
    % x_img = 0:x_time(end);
    % y_img = 1:50;
    % [X,Y] = meshgrid(x_time, y_img); %sample grid points
    % [Xq,Yq] = meshgrid(x_img, y_img);%querry grid points
    % phi_axial_interp = interp2(X,Y,phi_axial_all,Xq,Yq);
    % 
    % f = figure;
    % imagesc(x_img,y_img,phi_axial_interp)
    % colormap coolwarm
    % ax = gca;
    % ax.DataAspectRatio = [50 1 1];
    % title("k_{gj} = " + string(comb_vec(1,i_all)) + ", k_{chan}"+ string(comb_vec(2,i_all)))
    % box off
    % set(gca,'FontSize',10,'TickDir','out','YDir','normal');
    % xlabel("time (ms)")
    % ylabel("cell no")
    % clim([-90,45])
    % saveas(f,"/fs/scratch//PAS1622/nickmoise/ID_2025/processed_data/figures/" + string(i_all) + ".png")
    % xlim([9000,9300])
    

    % Compute CV
    % thresh_activation = 0;
    % 
    % ind_last_2beat = find(ts>(cycle_vec(i_cycle)*8-10));
    % phi_axial_2beat = phi_axial_all(:,ind_last_2beat);
    % x = ts(ind_last_2beat);
    % activ_time = zeros(50,2) .* NaN;
    % local_cv = [];
    % for i = 1:50
    %     data_iso_int = phi_axial_2beat(i,:);
    %     data_iso1 = data_iso_int(:,1:end-1);
    %     data_iso2 = data_iso_int(:,2:end);
    %     ind_th = find(data_iso1<thresh_activation & data_iso2> thresh_activation); % find indices below/above threshold
    %     y1 = data_iso1(ind_th); y2 = data_iso2(ind_th);
    %     m_slope = (y2-y1)./(x(ind_th+1)-x(ind_th));  % linear slope
    %     activ_time_local = x(ind_th) - (data_iso_int(ind_th)-thresh_activation)./m_slope;
    %     if length(activ_time_local) > 2
    %         activ_time_local = activ_time_local(1:2);
    %     end
    %     activ_time(i,:) = activ_time_local;
    % end
    % 
    % local_cv = (100./diff(activ_time(:,end)))./10; %cm/s
    % local_cv = local_cv(5:45); %cm/s
    
    % figure
    % plot(local_cv)
    % hold on

end



end

save(save_name,'cycle_vec','comb_vec','end_cycles')


end
%% CV block restart

data_dir = "/fs/scratch/PAS1622/nickmoise/ID_2025/" + ...
           "hetg_tissue_restart_50b_all_clampall/";

save_dir = "/fs/scratch/PAS1622/nickmoise/ID_2025/" + ...
           "processed_data/";

gj_vec = [8];
chan_vec = [8];
comb_vec = combvec(gj_vec, chan_vec);

mesh_mid = 2;

if mesh_mid == 1
    sim_root = "msh1_baseline_cycle%d_beats50_D01_gj_loc";
else
    sim_root = "msh1_baseline_msh2_p60_ip60_cycle%d_beats50_D01_gj_loc";
end

save_name = save_dir + "CV_resti_" + sprintf(sim_root,[]) + "_chan_loc";

% save_name = strrep(save_name, '.', ''); %remove dot to prevent file extension errors   
% cycle_vec = [200:5:300 350:50:1000];
cycle_vec = [205];
CV_cycles_all = zeros(length(comb_vec),length(cycle_vec));
beats_2_cycles_all = zeros(length(comb_vec),length(cycle_vec));

[~, loop_length] = size(comb_vec);

for i_all = 1:loop_length

file = sim_root ...
     + string(1) + "_" + string(comb_vec(1,i_all)) ...
     + "_chan_loc" ...
     + string(1) + "_" + string(comb_vec(2,i_all));
file = strrep(file, '.', '');


data_list = dir(data_dir);
data_list = natsortfiles(data_list);
data_list(1:2) = [];


thresh_activation = 0;


CV_cycles = zeros(1,length(cycle_vec));
beats_2_cycles = zeros(1,length(cycle_vec));
for i_cycle = 1:length(cycle_vec)
    % fprintf("%d / %d \n", i_file, length(data_list))
    % filename = data_dir + "/" +  data_list(i_file).name;

    fprintf("%d / %d \n", i_cycle, length(cycle_vec))
    load_file = sprintf(file, cycle_vec(i_cycle));
    filename = data_dir + load_file + ".mat";
    if isfile(filename)
        load(filename, 'Ca_cleft_mean_all','Na_cleft_mean_all','phi_axial_all','ts','ts_save')
    else
        disp("missing: " + filename)
        CV_cycles(i_cycle) = NaN;
        continue
    end

    phi_axial_all(:,end) = phi_axial_all(:,end-1);
    ind_last_2beat = find(ts>(cycle_vec(i_cycle)*8-10));
    phi_axial_2beat = phi_axial_all(:,ind_last_2beat);
    
    x = ts(ind_last_2beat);
    cell_index = [10,40]; % cells from which CV is comp
    activ_time = zeros(2,20);
    for i_cell = 1:length(cell_index)
        data_iso_int = phi_axial_2beat(cell_index(i_cell),:);
        data_iso1 = data_iso_int(:,1:end-1);
        data_iso2 = data_iso_int(:,2:end);
        ind_th = find(data_iso1<thresh_activation & data_iso2> thresh_activation); % find indices below/above threshold
        y1 = data_iso1(ind_th); y2 = data_iso2(ind_th);
        m_slope = (y2-y1)./(x(ind_th+1)-x(ind_th));  % linear slope
        activ_time_cell = x(ind_th) - (data_iso_int(ind_th)-thresh_activation)./m_slope;
        activ_time(i_cell,1:length(activ_time_cell)) = activ_time_cell;
    end
    
    beats1 = activ_time(1,:) ~= 0;
    beats1(beats1==0) = [];
    beats2 = activ_time(2,:) ~= 0;
    beats2(beats2==0) = [];
    
    cell_dist = (cell_index(2) - cell_index(1))*100; %100 = cell length in um

    
    
    if activ_time(2,2)>activ_time(1,2)
        beats_2_cycles(i_cycle) = 1;
        CV = cell_dist./(activ_time(2,2) - activ_time(1,2)) ./10;
        CV_cycles(i_cycle) = CV;
    else
        beats_2_cycles(i_cycle) = 0;
        CV = cell_dist./(activ_time(2,1) - activ_time(1,1)) ./10;
        CV_cycles(i_cycle) = CV;
    end
   
    
    
    % load(filename, 'scale_gj_loc', 'scale_chan_loc')
    % 
    % CV_mat(find(scale_gj_loc_vec == scale_gj_loc), ...
    %        find(scale_chan_loc_vec == scale_chan_loc)) = CV_cycles(i_file);


    % figure
    % imagesc(ts,1:50,Na_cleft_mean_all)
    % ax = gca;
    % ax.DataAspectRatio = [100 1 1];

    %%% phi axial
    x_time = ts;
    x_img = 0:x_time(end);
    y_img = 1:50;
    [X,Y] = meshgrid(x_time, y_img); %sample grid points
    [Xq,Yq] = meshgrid(x_img, y_img);%querry grid points
    phi_axial_interp = interp2(X,Y,phi_axial_all,Xq,Yq);

    if mesh_mid == 1
        title_mesh_string = "ID-baseline-";
    else
        title_mesh_string = "ID-wide-";
    end

    title_full = title_mesh_string + "k-{gj}-" + string(comb_vec(1,i_all)) ...
               + "-k-{na}-"+ string(comb_vec(2,i_all)) + "-cycle-" + string(cycle_vec(i_cycle));

    f = figure;
    imagesc(x_img,y_img,phi_axial_interp)
    colormap coolwarm
    ax = gca;
    ax.DataAspectRatio = [50 1 1];
    title(title_full)
    box off
    set(gca,'FontSize',10,'TickDir','out','YDir','normal');
    xlabel("time (ms)")
    ylabel("cell no")
    clim([-90,45])
    % saveas(f,"/fs/scratch//PAS1622/nickmoise/ID_2025/processed_data/figures/" + title_full + ".png")
    % close(f)
    % xlim([9000,9300])

    [~,peak_locs] = findpeaks(phi_axial_all(end,:),'MinPeakHeight',5,'MinPeakDistance',150);

    end_cycles(i_all, i_cycle) = length(peak_locs)
   

end

% plot(ts(ind_last_2beat),phi_axial_all(5,:))


% plot(cycle_vec,CV_cycles)
% hold on

% plot(Na_cleft_mean_all(25,:))

% 
% 



CV_cycles_all(i_all,:) = CV_cycles;
beats_2_cycles_all(i_all,:) = beats_2_cycles;
end

% save(save_name,'cycle_vec','comb_vec','CV_cycles_all','beats_2_cycles_all')


%%
% 
% %%
% figure
plot(ts_save,Ca_cleft(:,1:end),'color',[1 0 0 0.5])
ylim([1,3])
% %         
% figure
% plot(ts_save,Na_cleft(:,1:end),'color',[0 0 1 0.5])
% ylim([120,150])
% 
% figure
% plot(ts_save,phi_axial')



% for i = 1:25
%     phi = phi_axial(i,:);
%     phi(phi>-70) = 1;
%     phi(phi<=-70) = 0;
%     [~,peak_time] = findpeaks(phi);
% 
%     peak_list(i) = ts_save(peak_time(end));
% end

% figure
% hold on
% for i = 1:25
%     plot(ts_save,phi_axial(i,1:end))
% end


CV = 0.01./(diff(peak_list)./1000);
%         
%         figure
%         plot(ts_save,K_cleft(:,1:end-1),'color',[0 1 0 0.5])
%         ylim([4.5,8])
%         
%         figure
%         plot(ts_save,A_cleft(:,1:end-1),'color',[0 0 0 0.5])
%         ylim([140,150])
% 
figure
plot(ts_save,phi_axial(5,1:end))
% 
%         figure;
%         h = pcolor(ts, 1:length(iintra), phi_save(iintra,:));
%         set(gca,'ydir','reverse');
%         h.LineStyle = 'none';
%         set(gca,'ytick',1:length(iintra), 'yticklabel',iintra);
%         colorbar
%         caxis([-80,0])
% 
% subplot(2,2,3); plot(ts_save, phi_i,s); hold on;
% 
%         subplot(2,2,2);
%         h = pcolor(ts, 1:length(icleft), phi_save(icleft,:));
%         set(gca,'ydir','reverse');
%         h.LineStyle = 'none';
%         set(gca,'ytick',1:20:length(icleft), 'yticklabel',icleft(1:20:end));
% 
%         subplot(2,2,4); plot(ts, phi_cleft, s); hold on;

%         figure
%         plot(INa_all(ind_pre,:)'.*1e6)
%         ylim([-140,20])
%%



%%%%%%% start figure plots

%% figure 1 DATA
load("data/general_plots/" + ...
     "msh1_baseline_cycle1000_beats10_D1_gj_loc1_chan_loc1.mat")

Njunc = 25;

icleft = iEC(1:end-1);
iintra = setdiff(1:Nnodes-1,icleft);
phi_i = phi_save(iintra,:);
phi_cleft = phi_save(icleft,:);
% collect Vm / ionic currents
Vm = phi_save(Iind(:,1),:) - phi_save(Iind(:,2),:);

[~,ind] = sort(Iind(:,1));
Iind_cable = Iind(ind,:);
Vm_cable = Vm(ind,:);
tup = tup(ind,:); trepol = trepol(ind,:);

INa_all = I_save(p.iina:Ncurrents:end,:);

% sort axial / pre-/post-junctional membrane indices to match Vm order
[ind_axial, ~] = ind2sub([length(ind) length(indices.ind_axial)], find(ind == indices.ind_axial));
[ind_disc_pre, ~] = ind2sub([length(ind) length(indices.ind_disc_pre)], find(ind == indices.ind_disc_pre));
[ind_disc_post, ~] = ind2sub([length(ind) length(indices.ind_disc_post)], find(ind == indices.ind_disc_post));
INa_axial = INa_all(ind_axial,:);
INa_disc_pre = INa_all(ind_disc_pre,:);
INa_disc_post = INa_all(ind_disc_post,:);

ind_post = ind_disc_post((Njunc-1)*Mdisc+2:Njunc*Mdisc+1);
Vm_post = Vm_cable(ind_post,:); INa_post = INa_all(ind_post,:);
ind_pre = ind_disc_pre((Njunc-1)*Mdisc+1:Njunc*Mdisc);
Vm_pre = Vm_cable(ind_pre,:); INa_pre = INa_all(ind_pre,:);

Na_cleft_all = S_save(iEC,:);
Na_cleft = Na_cleft_all((Njunc-1)*Mdisc+1:Njunc*Mdisc,:);       
K_cleft_all = S_save(iEC+Nnodes,:);
K_cleft = K_cleft_all((Njunc-1)*Mdisc+1:Njunc*Mdisc,:);     
Ca_cleft_all = S_save(iEC+2*Nnodes,:);
Ca_cleft = Ca_cleft_all((Njunc-1)*Mdisc+1:Njunc*Mdisc,:); 
A_cleft_all = S_save(iEC+3*Nnodes,:);
A_cleft = A_cleft_all((Njunc-1)*Mdisc+1:Njunc*Mdisc,:);

tup_axial = tup(ind_axial(1:Nint:end),:);
i1 = round(.25*Ncell); i2 = round(.75*Ncell);
cv_est = 100*(i2-i1)*(p.L/1000)./(tup_axial(i2,:)-tup_axial(i1,:));  % mm/ms = m/s, 100*m/s = cm/s

phi_i = phi_save(iintra,:);
phi_axial = phi_i(ind_axial,:); 

%% Figure 1 FIGURE

% figure
% hold on
% for i = 1:50
%     plot(ts_save(1:end-1),phi_axial(i,:)+i*10)
% end


x_time = ts_save(1:end-1) - (ts_save(1));
y_loc = ones(1,length(x_time));


figure
for i = 1:50

    plot3(x_time,y_loc+i,phi_axial(i,:))
    hold on
    % plot3(x_time,y_loc+i+0.5,Ca_cleft_all.*50-40,'color',[0 0.447, 0.741,0.01])
    % plot3(x_time,y_loc+1,phi_axial(25,:))
end
daspect([1 0.05 1])
xlim([1,400])