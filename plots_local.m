clear

% α, β, γ, δ, ε, ζ, η, θ, ι, κ, λ, μ, ν, ξ, ο, π, ρ, σ, τ, υ, φ, χ, ψ, ω

%%%%%%% start figure plots
%% default colors for ions throughout
phi_color = [0.58, 0.403, 0.741];
Na_color = [.122, 0.467, 0.706];
Ca_color = [0.839, 0.153, 0.157];
K_color = [0.173, 0.627, 0.173];
A_color = [1.000, 0.498, 0.055];


%% figure baseline DATA
load_file = "data/hetg_tissue_block/" + ...
     "msh1_baseline_cycle1000_beats10_D1_gj_loc1_chan_loc1.mat";

%all except G_save - we don't need that for plots
load(load_file,'phi_save', 'I_save','S_save', 'bcl', 'D', 'scale_gj_loc', 'scale_chan_loc', ...
        'p', 'iEC', 'Nnodes', 'Ncell', 'Ncurrents', 'indices', 'Mdisc', ...
        'phi_axial_all', 'Iind', 'ts', 'model', 'FEM_file_list', 'tissue_legend',...
        'tup', 'trepol', 'ts_save', 'Nint');


icleft = iEC(1:end-1);
iintra = setdiff(1:Nnodes-1,icleft);
phi_i = phi_save(iintra,:);
phi_cleft_all = phi_save(icleft,:);
phi_axial = phi_save(Iind(1:Ncell,1),:); 
% collect Vm / ionic currents
Vm = phi_save(Iind(:,1),:) - phi_save(Iind(:,2),:);

INa_all = I_save(p.iina:Ncurrents:end,:);
ICa_all = I_save(p.iical:Ncurrents:end,:);

INa_axial = INa_all(1:Ncell,:);
ICa_axial = ICa_all(1:Ncell,:);

cleft1 = 25; 

ind_pre =  (Ncell+2 + (cleft1-1) * 2 * Mdisc )        : (Ncell+1 + (cleft1-1) * 2 * Mdisc + Mdisc);
ind_post = (Ncell+2 + (cleft1-1) * 2 * Mdisc + Mdisc) : (Ncell+1 + (cleft1-1) * 2 * Mdisc + 2*Mdisc);



Na_cleft_all = S_save(iEC,:);
K_cleft_all = S_save(iEC+Nnodes,:);
Ca_cleft_all = S_save(iEC+2*Nnodes,:);
A_cleft_all = S_save(iEC+3*Nnodes,:);

ion_cleft_index = (cleft1-1)*Mdisc+1:cleft1*Mdisc;

phi_cleft = phi_cleft_all((cleft1-1)*Mdisc+1:cleft1*Mdisc,:);
Na_cleft = Na_cleft_all((cleft1-1)*Mdisc+1:cleft1*Mdisc,:);
K_cleft = K_cleft_all((cleft1-1)*Mdisc+1:cleft1*Mdisc,:);
Ca_cleft = Ca_cleft_all((cleft1-1)*Mdisc+1:cleft1*Mdisc,:);
A_cleft = A_cleft_all((cleft1-1)*Mdisc+1:cleft1*Mdisc,:);


clear 'phi_save' 'I_save' 'S_save';
%% Figure baseline FIGURE

%D01
% phi_lim = [-107,8];
% Na_lim = [110,146];
% Ca_lim = [1.35,4];
% K_lim = [5,9];
% A_lim = [140,150];

%D1
phi_lim = [-84,8];
Na_lim = [115,146];
Ca_lim = [1.35,3.2];
K_lim = [5,7.8];
A_lim = [143,150];

%%% main figure
fig_save = figure;
% fig_save.Position = [0 0 1000 1200];
t = tiledlayout(5,2);
t.TileSpacing = 'compact';
t.Padding = 'compact';
fsize = 20;
legend_fsize = 12;

xlim_ap = [0,300];
xlim_upst = [11,15];
view_coords = [350,35];


x_time = ts_save - (ts_save(1))-45;
y_loc = zeros(1,length(x_time));


nexttile([5,1])
for i = flip(1:50) 
    hold on
    if i == 25 || i ==26
        plot3(x_time,y_loc+i,phi_axial(i,:),'LineWidth', 2, 'color', [0 0 0])
    else
        plot3(x_time,y_loc+i,phi_axial(i,:),'LineWidth', 2, 'color', [0.5 0.5 0.5])
    end
    % if i<50
    %     plot3(x_time,y_loc+i+0.5,mean(Ca_cleft_all((100*(i-1) + 1 : (100*i)),:)), ...
    %           'LineWidth', 2,'color','red')
    %     hold on
    % end
end
ax = gca;
ax.DataAspectRatio = [1 0.03 0.8];

% ax.XRuler.FirstCrossoverValue  = 0; % X crossover with Y axis
% ax.XRuler.SecondCrossoverValue  = 0; % X crossover with Z axis
% ax.YRuler.FirstCrossoverValue  = -50; % Y crossover with X axis
% ax.YRuler.SecondCrossoverValue  = -50; % Y crossover with Z axis
% ax.ZRuler.FirstCrossoverValue  = 0; % Z crossover with X axis
% ax.ZRuler.SecondCrossoverValue = 0; % Z crossover with Y axis
xticks(0:100:400)
yticks([1,10:10:50])
zticks([-80,-40,0,40])
xlim(xlim_ap)
ylim([1,50])
zlim([-90,50])
view(view_coords)
box off
set(gca,'FontSize',fsize,'TickDir','out');
xlabel("time (ms)")
ylabel("cell no")
zlabel("φ (mV)")




nexttile()
plot(x_time,phi_cleft,'LineWidth', 2, 'color', [phi_color 0.01])
hold on
plot(x_time,mean(phi_cleft,1),'LineWidth', 4, 'color', [phi_color]./1.5)
xlim(xlim_ap)
ylim(phi_lim);
box off
set(gca,'FontSize',fsize,'TickDir','out');
hold on
leg_color{1} = plot(nan,'Linewidth', 5, 'color', [phi_color]);
leg = legend([leg_color{:}], {'$\phi_{cleft}$'}, 'location', 'best','box','off');
legend boxoff 
set(leg, 'Interpreter','latex')
ax = gca;
ax.XAxis.Visible = 'off';
ylabel("φ (mV)")

nexttile()
plot(x_time,Na_cleft,'LineWidth', 2,'color',[Na_color 0.1])
hold on
plot(x_time,mean(Na_cleft,1),'LineWidth', 4,'color',[Na_color]./1.5)
xlim(xlim_ap)
ylim(Na_lim);
box off
set(gca,'FontSize',fsize,'TickDir','out');
hold on
leg_color{1} = plot(nan,'Linewidth', 5, 'color', [Na_color]);
leg = legend([leg_color{:}], {'$Na^+_{cleft}$'}, 'location', 'best','box','off');
legend boxoff 
set(leg, 'Interpreter','latex')
ax = gca;
ax.XAxis.Visible = 'off';
ylabel("Na^{+}(mM)")

% nexttile()
% plot(x_time,Na_cleft,'LineWidth', 2,'color',[Na_color 0.15])
% xlim(xlim_upst)
% ylim([118,146]);
% box off
% set(gca,'FontSize',fsize,'TickDir','out');
% ax = gca;
% ax.XAxis.Visible = 'off';
% ax.YAxis.Visible = 'off';

nexttile()
plot(x_time,Ca_cleft,'LineWidth', 2,'color',[Ca_color 0.1])
hold on
plot(x_time,mean(Ca_cleft,1),'LineWidth', 4,'color',[Ca_color]./1.5)
xlim(xlim_ap)
ylim(Ca_lim);
box off
set(gca,'FontSize',fsize,'TickDir','out');
hold on
leg_color{1} = plot(nan,'Linewidth', 5, 'color', [Ca_color]);
leg = legend([leg_color{:}], {'$Ca^{2+}_{cleft}$'}, 'location', 'best','box','off');
legend boxoff 
set(leg, 'Interpreter','latex')
ax = gca;
ax.XAxis.Visible = 'off';
ylabel("Ca^{2+}(mM)")

% nexttile()
% plot(x_time,Ca_cleft,'LineWidth', 2,'color',[Ca_color 0.15])
% xticks(0:100:400)
% xlim(xlim_upst)
% ylim([1.35,3.2]);
% box off
% set(gca,'FontSize',fsize,'TickDir','out');
% ax = gca;
% ax.XAxis.Visible = 'off';
% ax.YAxis.Visible = 'off';

nexttile()
plot(x_time,K_cleft,'LineWidth', 2,'color',[K_color 0.1])
hold on
plot(x_time,mean(K_cleft,1),'LineWidth', 4,'color',[K_color]./2)
xlim(xlim_ap)
ylim(K_lim);
box off
set(gca,'FontSize',fsize,'TickDir','out');
hold on
leg_color{1} = plot(nan,'Linewidth', 5, 'color', [K_color]);
leg = legend([leg_color{:}], {'$K^{+}_{cleft}$'}, 'location', 'best','box','off');
legend boxoff 
set(leg, 'Interpreter','latex')
ax = gca;
ax.XAxis.Visible = 'off';
ylabel("K^+(mM)")

% nexttile()
% plot(x_time,K_cleft,'LineWidth', 2,'color',[Ka_color 0.15])
% xlim(xlim_upst)
% ylim([5,8]);
% box off
% set(gca,'FontSize',fsize,'TickDir','out');
% ax = gca;
% ax.XAxis.Visible = 'off';
% ax.YAxis.Visible = 'off';

nexttile()
plot(x_time,A_cleft,'LineWidth', 2,'color',[A_color 0.15])
hold on
plot(x_time,mean(A_cleft,1),'LineWidth', 4,'color',[A_color] ./1.5)
xticks(0:100:400)
xlim(xlim_ap)
ylim(A_lim)
box off
set(gca,'FontSize',fsize,'TickDir','out');
hold on
leg_color{1} = plot(nan,'Linewidth', 5, 'color', [A_color]);
leg = legend([leg_color{:}], {'$A^{-}_{cleft}$'}, 'location', 'best','box','off');
legend boxoff 
set(leg, 'Interpreter','latex')
ylabel("A^-(mM)")
xlabel("time (ms)")

% nexttile()
% plot(x_time,A_cleft,'LineWidth', 2,'color',[A_color 0.2])
% xticks(0:20)
% xlim(xlim_upst)
% ylim([143.5,150])
% xticks([11,15,20])
% box off
% set(gca,'FontSize',fsize,'TickDir','out');
% ax = gca;
% ax.YAxis.Visible = 'off';
% 
% 


%%%phi_cleft inset
time_scale_bar_x = 13.5:0.1:14;
time_scale_bar_y = zeros(length(time_scale_bar_x),1) -60;
txt = string(time_scale_bar_x(end) - time_scale_bar_x(1)) + " ms";
txt_x = (time_scale_bar_x(1)*0.945 + time_scale_bar_x(end))/2;
txt_y = time_scale_bar_y(1)+8;

figure
plot(x_time,phi_cleft,'LineWidth', 2, 'color', [phi_color 0.15])
hold on
plot(x_time,mean(phi_cleft,1),'LineWidth', 4, 'color', [phi_color]./1.5)
plot(time_scale_bar_x,time_scale_bar_y,"LineWidth", 15, 'color', [0 0 0])
text(txt_x,txt_y,txt,'FontSize',fsize+10)

xticks(xlim_upst(1):1:xlim_upst(end))
xlim(xlim_upst)
ylim(phi_lim);
box off
set(gca,'FontSize',fsize,'TickDir','out');
ax = gca;
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';




%% figure hetg tissue D1 DATA
load_file = "data/hetg_tissue/" + ...
     "msh1_baseline_msh2_p60_ip60_cycle1000_beats10_D1_gj_loc1_chan_loc1";

%all except G_save - we don't need that for plots
load(load_file,'phi_save', 'I_save','S_save', 'bcl', 'D', 'scale_gj_loc', 'scale_chan_loc', ...
        'p', 'iEC', 'Nnodes', 'Ncell', 'Ncurrents', 'indices', 'Mdisc', ...
        'phi_axial_all', 'Na_cleft_mean_all', 'K_cleft_mean_all', 'Ca_cleft_mean_all', 'A_cleft_mean_all',...
        'Iind', 'ts', 'model', 'FEM_file_list', 'tissue_legend',...
        'tup', 'trepol', 'ts_save', 'Nint');

%current indices are according to Iind:
%first Ncell = axial
%one bulk, one bd
%100 pre/100 post etc

%phi indices are interior - membrane - cleft - membrane - interior etc

icleft = iEC(1:end-1);
iintra = setdiff(1:Nnodes-1,icleft);
phi_i = phi_save(iintra,:);
phi_cleft_all = phi_save(icleft,:);
phi_axial = phi_save(Iind(1:Ncell,1),:); 
% collect Vm / ionic currents
Vm = phi_save(Iind(:,1),:) - phi_save(Iind(:,2),:);

INa_all = I_save(p.iina:Ncurrents:end,:);
ICa_all = I_save(p.iical:Ncurrents:end,:);

INa_axial = INa_all(1:Ncell,:);
ICa_axial = ICa_all(1:Ncell,:);

cleft1 = 25; 

ind_pre =  (Ncell+2 + (cleft1-1) * 2 * Mdisc )        : (Ncell+1 + (cleft1-1) * 2 * Mdisc + Mdisc);
ind_post = (Ncell+2 + (cleft1-1) * 2 * Mdisc + Mdisc) : (Ncell+1 + (cleft1-1) * 2 * Mdisc + 2*Mdisc);



Na_cleft_all = S_save(iEC,:);
K_cleft_all = S_save(iEC+Nnodes,:);
Ca_cleft_all = S_save(iEC+2*Nnodes,:);
A_cleft_all = S_save(iEC+3*Nnodes,:);

ion_cleft_index = (cleft1-1)*Mdisc+1:cleft1*Mdisc;


clear 'phi_save' 'I_save' 'S_save';



%% Figure hetg tissue
%%% main figure

figure;
% fig_save.Position = [0 0 1000 1200];
t = tiledlayout(4,3);
t.TileSpacing = 'tight';
t.Padding = 'tight';
fsize = 20;
legend_fsize = 12;

cleft1 = 25;
cleft2 = 10;

xlim_ap = [0,300];
if D == 1
    xlim_upst_cleft1 = [0,3] + 12.5; 
    xlim_upst_cleft2 = [0,3] + 8.5;
elseif D == 0.1
    xlim_upst_cleft1 = [0,3] + 19.75; 
    xlim_upst_cleft2 = [0,3] + 10.25;
end
view_coords = [350,35];
tissue_yticks = [1,10:10:50];
voltage_ticks = [-80,-40,0,40];
voltage_lim = [-90,50];

ylim_INa = [-3.75e-4,1.5e-4];
ylim_Na = [110,145];
ticks_Na = [120:10:140];

ylim_ICa = [-18e-7,1.5e-7];
ylim_Ca = [1.5,4.0];
ticks_Ca = [2,3,4];

x_time = ts_save - (ts_save(1))-45;
y_loc = zeros(1,length(x_time));

%cable voltage
nexttile([4,1])
for i = flip(1:50) 
    hold on
    if i == cleft1 || i ==cleft1+1 || i==cleft2 || i == cleft2+1
        plot3(x_time,y_loc+i,phi_axial(i,:),'LineWidth', 2, 'color', [0 0 0])
    else
        plot3(x_time,y_loc+i,phi_axial(i,:),'LineWidth', 2, 'color', [0.5 0.5 0.5])
    end
    % if i<50
    %     plot3(x_time,y_loc+i+0.5,mean(Ca_cleft_all((100*(i-1) + 1 : (100*i)),:)), ...
    %           'LineWidth', 2,'color','red')
    %     hold on
    % end
end
ax = gca;
ax.DataAspectRatio = [1 0.03 0.8];

% ax.XRuler.FirstCrossoverValue  = 0; % X crossover with Y axis
% ax.XRuler.SecondCrossoverValue  = 0; % X crossover with Z axis
% ax.YRuler.FirstCrossoverValue  = -50; % Y crossover with X axis
% ax.YRuler.SecondCrossoverValue  = -50; % Y crossover with Z axis
% ax.ZRuler.FirstCrossoverValue  = 0; % Z crossover with X axis
% ax.ZRuler.SecondCrossoverValue = 0; % Z crossover with Y axis

yticks(tissue_yticks)
zticks(voltage_ticks)
xlim(xlim_ap)
ylim([1,50])
zlim(voltage_lim)
view(view_coords)
box off
set(gca,'FontSize',fsize,'TickDir','out');
% xlabel("time (ms)")
% ylabel("cell no")
% zlabel("φ (mV)")

%%%%% CLEFT 1
ind_pre =  (Ncell+2 + (cleft1-1) * 2 * Mdisc )        : (Ncell+1 + (cleft1-1) * 2 * Mdisc + Mdisc);
ind_post = (Ncell+2 + (cleft1-1) * 2 * Mdisc + Mdisc) : (Ncell+1 + (cleft1-1) * 2 * Mdisc + 2*Mdisc);


%Na conc
nexttile()
plot(x_time,phi_axial(cleft1,:),'LineWidth', 2, 'color', [0 0 0])
hold on
plot(x_time,phi_axial(cleft1+1,:),'LineWidth', 2, 'color', [0 0 0])
yticks(voltage_ticks)
ylim(voltage_lim)
ylabel("φ (mV)")

yyaxis right
plot(x_time,Na_cleft_all((cleft1-1)*Mdisc+1:cleft1*Mdisc,:),...
    '-','LineWidth', 2, 'color', [Na_color, 0.25]);
ylim(ylim_Na)
yticks(ticks_Na)
ylabel("Na^{+}(mM)")

xlim(xlim_upst_cleft1)
ax = gca;
box off
set(gca,'FontSize',fsize,'TickDir','out');
ax.YAxis(1).Color = [0 0 0];
ax.YAxis(2).Color = Na_color;
ax.XAxis.Visible = 'off';


leg_color = {};
leg_color{1} = plot(nan,'-','Linewidth', 5, 'color', [Na_color]);
legend([leg_color{:}], {'Na^+_{cleft}'}, 'location', 'best','box','off');


%Na currents
nexttile()
plot(x_time,phi_axial(cleft1,:),'LineWidth', 2, 'color', [0 0 0])
hold on
plot(x_time,phi_axial(cleft1+1,:),'LineWidth', 2, 'color', [0 0 0])
yticks(voltage_ticks)
ylim(voltage_lim)

yyaxis right
plot(x_time,INa_all(ind_pre,:),...
    '-','LineWidth', 2, 'color', [Na_color, 0.2]);
hold on
plot(x_time,INa_all(ind_post,:),...
    '-','LineWidth', 2, 'color', [1-Na_color, 0.2]);

plot(x_time,mean(INa_all(ind_pre,:),1),...
    '-','LineWidth', 4, 'color', [Na_color-0.1]);
plot(x_time,mean(INa_all(ind_post,:)),...
    '-','LineWidth', 4, 'color', [(1-Na_color).*0.65]);
ylim(ylim_INa)


xlim(xlim_upst_cleft1)
ax = gca;
box off
set(gca,'FontSize',fsize,'TickDir','out');
ax.YAxis(1).Color = [0 0 0];
ax.YAxis(2).Color = [0 0 0];
ax.XAxis.Visible = 'off';
ax.YAxis(1).Visible = 'off';
leg_color = {};
leg_color{1} = plot(nan,'-','Linewidth', 5, 'color', [Na_color]);
leg_color{2} = plot(nan,'-','Linewidth', 5, 'color', 1-[Na_color]);
legend([leg_color{:}], {'I_{Na}^{pre}','I_{Na}^{post}'}, 'location', 'best','box','off');
ylabel("I_{Na} (μA/μF)")


%Ca conc
nexttile()
plot(x_time,phi_axial(cleft1,:),'LineWidth', 2, 'color', [0 0 0])
hold on
plot(x_time,phi_axial(cleft1+1,:),'LineWidth', 2, 'color', [0 0 0])
yticks(voltage_ticks)
ylim(voltage_lim)
ylabel("φ (mV)")

yyaxis right
plot(x_time,Ca_cleft_all((cleft1-1)*Mdisc+1:cleft1*Mdisc,:),...
    '-','LineWidth', 2, 'color', [Ca_color, 0.25]);
ylim(ylim_Ca)
yticks(ticks_Ca)
ylabel("Ca^{2+}(mM)")



xlim(xlim_upst_cleft1)
ax = gca;
box off
set(gca,'FontSize',fsize,'TickDir','out');
ax.YAxis(1).Color = [0 0 0];
ax.YAxis(2).Color = Ca_color;
ax.XAxis.Visible = 'off';

leg_color = {};
leg_color{1} = plot(nan,'-','Linewidth', 5, 'color', [Ca_color]);
legend([leg_color{:}], {'Ca^{2+}_{cleft}'}, 'location', 'best','box','off');


%Ca current
nexttile()
plot(x_time,phi_axial(cleft1,:),'LineWidth', 2, 'color', [0 0 0])
hold on
plot(x_time,phi_axial(cleft1+1,:),'LineWidth', 2, 'color', [0 0 0])
yticks(voltage_ticks)
ylim(voltage_lim)

yyaxis right
plot(x_time,ICa_all(ind_pre,:),...
    '-','LineWidth', 2, 'color', [Ca_color, 0.15]);
hold on
plot(x_time,ICa_all(ind_post,:),...
    '-','LineWidth', 2, 'color', [1-Ca_color, 0.15]);

plot(x_time,mean(ICa_all(ind_pre,:),1),...
    '-','LineWidth', 4, 'color', [Ca_color-0.1]);
plot(x_time,mean(ICa_all(ind_post,:)),...
    '-','LineWidth', 4, 'color', [(1-Ca_color).*0.55]);
ylim(ylim_ICa)
ylabel("I_{Ca} (μA/μF)")

xlim(xlim_upst_cleft1)
ax = gca;
box off
set(gca,'FontSize',fsize,'TickDir','out');
ax.YAxis(1).Color = [0 0 0];
ax.YAxis(2).Color = [0 0 0];
ax.XAxis.Visible = 'off';
ax.YAxis(1).Visible = 'off';

leg_color = {};
leg_color{1} = plot(nan,'-','Linewidth', 5, 'color', [Ca_color]);
leg_color{2} = plot(nan,'-','Linewidth', 5, 'color', 1-[Ca_color]);
legend([leg_color{:}], {'I_{Ca}^{pre}','I_{Ca}^{post}'}, 'location', 'best','box','off');



%%%%% CLEFT2
ind_pre =  (Ncell+2 + (cleft2-1) * 2 * Mdisc )        : (Ncell+1 + (cleft2-1) * 2 * Mdisc + Mdisc);
ind_post = (Ncell+2 + (cleft2-1) * 2 * Mdisc + Mdisc) : (Ncell+1 + (cleft2-1) * 2 * Mdisc + 2*Mdisc);

% Na conc
nexttile()
plot(x_time,phi_axial(cleft2,:),'LineWidth', 2, 'color', [0 0 0])
hold on
plot(x_time,phi_axial(cleft2+1,:),'LineWidth', 2, 'color', [0 0 0])
yticks(voltage_ticks)
ylim(voltage_lim)
ylabel("φ (mV)")

yyaxis right
plot(x_time,Na_cleft_all((cleft2-1)*Mdisc+1:cleft2*Mdisc,:),...
    '-','LineWidth', 2, 'color', [Na_color, 0.25]);
ylim(ylim_Na)
yticks(ticks_Na)
ylabel("Na^{+}(mM)")



xlim(xlim_upst_cleft2)
ax = gca;
box off
set(gca,'FontSize',fsize,'TickDir','out');
ax.YAxis(1).Color = [0 0 0];
ax.YAxis(2).Color = Na_color;
ax.XAxis.Visible = 'off';

leg_color = {};
leg_color{1} = plot(nan,'-','Linewidth', 5, 'color', [Na_color]);
legend([leg_color{:}], {'Na^+_{cleft}'}, 'location', 'best','box','off');


%Na currents
nexttile()
time_scale_bar_x = (xlim_upst_cleft2(1)+2):0.1:(xlim_upst_cleft2(end)-0.5);
time_scale_bar_y = zeros(length(time_scale_bar_x),1) - 60;
txt_x = (time_scale_bar_x(1)*0.955 + time_scale_bar_x(end))/2;
txt_y = time_scale_bar_y(1)+14;

plot(x_time,phi_axial(cleft2,:),'LineWidth', 2, 'color', [0 0 0])
hold on
plot(x_time,phi_axial(cleft2+1,:),'LineWidth', 2, 'color', [0 0 0])
yticks(voltage_ticks)
ylim(voltage_lim)
plot(time_scale_bar_x, time_scale_bar_y,"LineWidth", 6, 'color', [0 0 0])
txt = string(time_scale_bar_x(end) - time_scale_bar_x(1)) + " ms";
text(txt_x,txt_y,txt,'FontSize',fsize)

yyaxis right
plot(x_time,INa_all(ind_pre,:),...
    '-','LineWidth', 2, 'color', [Na_color, 0.2]);
hold on
plot(x_time,INa_all(ind_post,:),...
    '-','LineWidth', 2, 'color', [1-Na_color, 0.2]);

plot(x_time,mean(INa_all(ind_pre,:),1),...
    '-','LineWidth', 4, 'color', [Na_color-0.1]);
plot(x_time,mean(INa_all(ind_post,:)),...
    '-','LineWidth', 4, 'color', [(1-Na_color).*0.65]);
ylim(ylim_INa)
ylabel("I_{Na} (μA/μF)")

xlim(xlim_upst_cleft2)
ax = gca;
box off
set(gca,'FontSize',fsize,'TickDir','out');
ax.YAxis(1).Color = [0 0 0];
ax.YAxis(2).Color = [0 0 0];
ax.XAxis.Visible = 'off';
ax.YAxis(1).Visible = 'off';

leg_color = {};
leg_color{1} = plot(nan,'-','Linewidth', 5, 'color', [Na_color]);
leg_color{2} = plot(nan,'-','Linewidth', 5, 'color', 1-[Na_color]);
legend([leg_color{:}], {'I_{Na}^{pre}','I_{Na}^{post}'}, 'location', 'best','box','off');

%Ca conc
nexttile()
plot(x_time,phi_axial(cleft2,:),'LineWidth', 2, 'color', [0 0 0])
hold on
plot(x_time,phi_axial(cleft2+1,:),'LineWidth', 2, 'color', [0 0 0])
yticks(voltage_ticks)
ylim(voltage_lim)
ylabel("φ (mV)")


yyaxis right
plot(x_time,Ca_cleft_all((cleft2-1)*Mdisc+1:cleft2*Mdisc,:),...
    '-','LineWidth', 2, 'color', [Ca_color, 0.25]);
ylim(ylim_Ca)
yticks(ticks_Ca)
ylabel("Ca^{2+}(mM)")

xlim(xlim_upst_cleft2)
ax = gca;
box off
set(gca,'FontSize',fsize,'TickDir','out');
ax.YAxis(1).Color = [0 0 0];
ax.YAxis(2).Color = Ca_color;
ax.XAxis.Visible = 'off';

leg_color = {};
leg_color{1} = plot(nan,'-','Linewidth', 5, 'color', [Ca_color]);
legend([leg_color{:}], {'Ca^{2+}_{cleft}'}, 'location', 'best','box','off');

%Ca current
nexttile()
plot(x_time,phi_axial(cleft2,:),'LineWidth', 2, 'color', [0 0 0])
hold on
plot(x_time,phi_axial(cleft2+1,:),'LineWidth', 2, 'color', [0 0 0])
yticks(voltage_ticks)
ylim(voltage_lim)

yyaxis right
plot(x_time,ICa_all(ind_pre,:),...
    '-','LineWidth', 2, 'color', [Ca_color, 0.15]);
hold on
plot(x_time,ICa_all(ind_post,:),...
    '-','LineWidth', 2, 'color', [1-Ca_color, 0.15]);

plot(x_time,mean(ICa_all(ind_pre,:),1),...
    '-','LineWidth', 4, 'color', [Ca_color-0.1]);
plot(x_time,mean(ICa_all(ind_post,:)),...
    '-','LineWidth', 4, 'color', [(1-Ca_color).*0.55]);
ylim(ylim_ICa)
ylabel("I_{Ca} (μA/μF)")

xlim(xlim_upst_cleft2)
ax = gca;
box off
set(gca,'FontSize',fsize,'TickDir','out');
ax.YAxis(1).Color = [0 0 0];
ax.YAxis(2).Color = [0 0 0];
ax.XAxis.Visible = 'off';
ax.YAxis(1).Visible = 'off';

leg_color = {};
leg_color{1} = plot(nan,'-','Linewidth', 5, 'color', [Ca_color]);
leg_color{2} = plot(nan,'-','Linewidth', 5, 'color', 1-[Ca_color]);
legend([leg_color{:}], {'I_{Ca}^{pre}','I_{Ca}^{post}'}, 'location', 'best','box','off');

%
%%%% CV along tissue
thresh_activation = -10;

x = ts_save;
activ_time = [];
for i = 1:50
    data_iso_int = phi_axial(i,:);
    data_iso1 = data_iso_int(:,1:end-1);
    data_iso2 = data_iso_int(:,2:end);
    ind_th = find(data_iso1<thresh_activation & data_iso2> thresh_activation); % find indices below/above threshold
    y1 = data_iso1(ind_th); y2 = data_iso2(ind_th);
    m_slope = (y2-y1)./(x(ind_th+1)-x(ind_th));  % linear slope
    activ_time(i) = x(ind_th) - (data_iso_int(ind_th)-thresh_activation)./m_slope;
end

local_cv = (100./diff(activ_time))./10; %cm/s
local_cv = local_cv(5:45); %cm/s


figure
fsize = 30;
plot(5:45,local_cv, 'LineWidth',8)
box off
set(gca,'FontSize',fsize,'TickDir','out');
xlim([1,50])
xlabel('cell no')
ylabel("CV (cm/s)")
xticks([1,10:10:50])


%% current comp




cleft1 = 25;
cleft2 = 10;
% xlim_ap = [0,300];
% if D == 1
%     xlim_upst_cleft1 = [0,3] + 12.5; 
%     xlim_upst_cleft2 = [0,3] + 8.5;
% elseif D == 0.1
%     xlim_upst_cleft1 = [0,3] + 19.75; 
%     xlim_upst_cleft2 = [0,3] + 10.25;
% end


figure;
t = tiledlayout(2,2);
t.TileSpacing = 'tight';
t.Padding = 'tight';
fsize = 20;
legend_fsize = 12;

load_file = "data/hetg_tissue/" + ...
     "msh1_baseline_msh2_p60_ip60_cycle1000_beats10_D1_gj_loc1_chan_loc1";

%all except G_save - we don't need that for plots
load(load_file,'I_save','D', 'Ncell', 'Ncurrents', 'Mdisc', 'ts_save','p');
ICa_all = I_save(p.iical:Ncurrents:end,:);
INa_all = I_save(p.iina:Ncurrents:end,:);
x_time = ts_save - (ts_save(1))-45 - 8;


%%% D1 INa
nexttile(1)
title("Normal GJ coupling")
hold on

ind_pre =  (Ncell+2 + (cleft1-1) * 2 * Mdisc )        : (Ncell+1 + (cleft1-1) * 2 * Mdisc + Mdisc);
ind_post = (Ncell+2 + (cleft1-1) * 2 * Mdisc + Mdisc) : (Ncell+1 + (cleft1-1) * 2 * Mdisc + 2*Mdisc);
INa_total = INa_all(ind_pre,:) + INa_all(ind_post,:);
plot(x_time,mean(INa_total,1),'Linewidth',2,'color',Na_color);
min1 = min(mean(INa_total,1),[],'all');


ind_pre =  (Ncell+2 + (cleft2-1) * 2 * Mdisc )        : (Ncell+1 + (cleft2-1) * 2 * Mdisc + Mdisc);
ind_post = (Ncell+2 + (cleft2-1) * 2 * Mdisc + Mdisc) : (Ncell+1 + (cleft2-1) * 2 * Mdisc + 2*Mdisc);
INa_total = INa_all(ind_pre,:) + INa_all(ind_post,:);
plot(x_time,mean(INa_total,1),':','Linewidth',2,'color',Na_color);
min2 = min(mean(INa_total,1),[],'all');

legend("I_{Na}, baseline", "I_{Na}, wide ID",'location','best')

xlabel("time (ms)")
ylabel("I_{Na} (μA/μF)")
box off
set(gca,'FontSize',fsize,'TickDir','out');
xlim([0,20])
ylim([-1.2e-4,0.1e-4])
% xlim([0,250])
% ylim([-1.3e-6,0])


%%% D1 ICa
nexttile(3)
hold on

ind_pre =  (Ncell+2 + (cleft1-1) * 2 * Mdisc )        : (Ncell+1 + (cleft1-1) * 2 * Mdisc + Mdisc);
ind_post = (Ncell+2 + (cleft1-1) * 2 * Mdisc + Mdisc) : (Ncell+1 + (cleft1-1) * 2 * Mdisc + 2*Mdisc);
ICa_total = ICa_all(ind_pre,:) + ICa_all(ind_post,:);
plot(x_time,mean(ICa_total,1),'Linewidth',2,'color',Ca_color);
min1 = min(mean(ICa_total,1),[],'all');


ind_pre =  (Ncell+2 + (cleft2-1) * 2 * Mdisc )        : (Ncell+1 + (cleft2-1) * 2 * Mdisc + Mdisc);
ind_post = (Ncell+2 + (cleft2-1) * 2 * Mdisc + Mdisc) : (Ncell+1 + (cleft2-1) * 2 * Mdisc + 2*Mdisc);
ICa_total = ICa_all(ind_pre,:) + ICa_all(ind_post,:);
plot(x_time,mean(ICa_total,1),':','Linewidth',2,'color',Ca_color);
min2 = min(mean(ICa_total,1),[],'all');

legend("I_{Ca}, baseline", "I_{Ca}, wide ID",'location','best')

xlabel("time (ms)")
ylabel("I_{Ca} (μA/μF)")
box off
set(gca,'FontSize',fsize,'TickDir','out');
% xlim([0,50])
% ylim([-1.2e-4,0])
xlim([0,250])
ylim([-1.3e-6,0.1e-6])



load_file = "data/hetg_tissue/" + ...
     "msh1_baseline_msh2_p60_ip60_cycle1000_beats10_D01_gj_loc1_chan_loc1";

%all except G_save - we don't need that for plots
load(load_file,'I_save','D', 'Ncell', 'Ncurrents', 'Mdisc', 'ts_save','p');
ICa_all = I_save(p.iical:Ncurrents:end,:);
INa_all = I_save(p.iina:Ncurrents:end,:);
x_time = ts_save - (ts_save(1))-45 - 10;

%%% D01 INa
nexttile(2)
title("Low GJ coupling")
hold on

ind_pre =  (Ncell+2 + (cleft1-1) * 2 * Mdisc )        : (Ncell+1 + (cleft1-1) * 2 * Mdisc + Mdisc);
ind_post = (Ncell+2 + (cleft1-1) * 2 * Mdisc + Mdisc) : (Ncell+1 + (cleft1-1) * 2 * Mdisc + 2*Mdisc);
INa_total = INa_all(ind_pre,:) + INa_all(ind_post,:);
plot(x_time,mean(INa_total,1),'Linewidth',2,'color',Na_color);
min1 = min(mean(INa_total,1),[],'all');


ind_pre =  (Ncell+2 + (cleft2-1) * 2 * Mdisc )        : (Ncell+1 + (cleft2-1) * 2 * Mdisc + Mdisc);
ind_post = (Ncell+2 + (cleft2-1) * 2 * Mdisc + Mdisc) : (Ncell+1 + (cleft2-1) * 2 * Mdisc + 2*Mdisc);
INa_total = INa_all(ind_pre,:) + INa_all(ind_post,:);
plot(x_time,mean(INa_total,1),':','Linewidth',2,'color',Na_color);
min2 = min(mean(INa_total,1),[],'all');

legend("I_{Na}, baseline", "I_{Na}, wide ID",'location','best')

xlabel("time (ms)")
ylabel("I_{Na} (μA/μF)")
box off
set(gca,'FontSize',fsize,'TickDir','out');
ax = gca;
ax.YAxis.Visible = 'off';
xlim([0,20])
ylim([-1.2e-4,0.1e-4])
% xlim([0,250])
% ylim([-1.3e-6,0])



%%% D01 ICa
nexttile(4)
hold on

ind_pre =  (Ncell+2 + (cleft1-1) * 2 * Mdisc )        : (Ncell+1 + (cleft1-1) * 2 * Mdisc + Mdisc);
ind_post = (Ncell+2 + (cleft1-1) * 2 * Mdisc + Mdisc) : (Ncell+1 + (cleft1-1) * 2 * Mdisc + 2*Mdisc);
ICa_total = ICa_all(ind_pre,:) + ICa_all(ind_post,:);
plot(x_time,mean(ICa_total,1),'Linewidth',2,'color',Ca_color);
min1 = min(mean(ICa_total,1),[],'all');


ind_pre =  (Ncell+2 + (cleft2-1) * 2 * Mdisc )        : (Ncell+1 + (cleft2-1) * 2 * Mdisc + Mdisc);
ind_post = (Ncell+2 + (cleft2-1) * 2 * Mdisc + Mdisc) : (Ncell+1 + (cleft2-1) * 2 * Mdisc + 2*Mdisc);
ICa_total = ICa_all(ind_pre,:) + ICa_all(ind_post,:);
plot(x_time,mean(ICa_total,1),':','Linewidth',2,'color',Ca_color);
min2 = min(mean(ICa_total,1),[],'all');

legend("I_{Ca}, baseline", "I_{Ca}, wide ID",'location','best')

xlabel("time (ms)")
ylabel("I_{Ca} (μA/μF)")
box off
set(gca,'FontSize',fsize,'TickDir','out');
ax = gca;
ax.YAxis.Visible = 'off';
% xlim([0,50])
% ylim([-1.2e-4,0])
xlim([0,250])
ylim([-1.3e-6,0.1e-6])



ratio = min2/min1
%% test all mean cleft currents
figure
for i = (1:49)
    ind_pre = (Ncell+2 + (i-1) * 2 * Mdisc ) : (Ncell+1 + (i-1) * 2 * Mdisc + Mdisc);
    ind_post = (Ncell+2 + (i-1) * 2 * Mdisc + Mdisc) : (Ncell+1 + (i-1) * 2 * Mdisc + 2*Mdisc);

    plot3(x_time,y_loc+i,mean(ICa_all(ind_pre,:),1) + mean(ICa_all(ind_post,:),1), ...
        'LineWidth', 2,'color',[[Ca_color]./1.25, 0.5])
    hold on
    % plot3(x_time,y_loc+i,mean(INa_all(ind_pre,:),1),'LineWidth', 2,'color',[[Ca_color], 1])
    % hold on
    % plot3(x_time,y_loc+i+0.5,mean(INa_all(ind_post,:),1),'LineWidth', 2,'color',[1-[Ca_color]./1.25,1])
    % plot3(x_time,y_loc+i,cumsum(mean(ICa_all(ind_pre,:)+mean(ICa_all(ind_post,:)),1)),...
    %       'LineWidth', 2,'color',[[Ca_color]./1.25, 0.5])

    hold on

end


% mesh(INa_all)


ax = gca;
% ax.DataAspectRatio = [1 0.03 1e-9];
% ax.DataAspectRatio = [1 0.03 5e-7];

% ax.XRuler.FirstCrossoverValue  = 0; % X crossover with Y axis
% ax.XRuler.SecondCrossoverValue  = 0; % X crossover with Z axis
% ax.YRuler.FirstCrossoverValue  = -50; % Y crossover with X axis
% ax.YRuler.SecondCrossoverValue  = -50; % Y crossover with Z axis
% ax.ZRuler.FirstCrossoverValue  = 0; % Z crossover with X axis
% ax.ZRuler.SecondCrossoverValue = 0; % Z crossover with Y axis
xticks(0:100:400)
yticks(tissue_yticks)
% zticks([-80,-40,0,40])
xlim(xlim_ap)
ylim([1,50])
% zlim([-1e-4,0])
% zlim([-1e-6,0])
% zlim([-1.3e-2,0])
% zlim([-0.15e-2,0])
% zlim([-3e-4,0])

view([300 35])
box off
set(gca,'FontSize',fsize,'TickDir','out');
% title("integral INa cleft D01")


%% figure CV 7x7 maps
fsize = 16;

load("data/processed_data/CV_map7x7_msh1_baseline_cycle1000_beats10_D1_gj_loc_chan_loc")
CV_map_msh1 = CV_map;

load("data/processed_data/CV_map7x7_na_msh1_baseline_cycle1000_beats10_D1_gj_loc_chan_loc")
CV_map_msh2 = CV_map;

clear CV_map

scale_gj_loc_vec = [8 4 2 1 1/2 1/4 1/8];
scale_chan_loc_vec = [8 4 2 1 1/2 1/4 1/8];
[X,Y] = meshgrid(scale_gj_loc_vec,scale_chan_loc_vec);
xtick_values = flip(scale_gj_loc_vec); %

% figure with 2 colormaps
nColors = 256;
cmap_colors_msh1 = get_2dcolormap_value(CV_map_msh1(1:end-1,1:end-1), viridis(nColors));
cmap_colors_msh2 = get_2dcolormap_value(CV_map_msh1(1:end-1,1:end-1), inferno(nColors));

figure;
s = surf(X, Y, CV_map_msh1);
set(s, 'FaceColor', 'flat', ...
       'CData', cmap_colors_msh1, ...
       'CDataMapping', 'direct', ...
       'EdgeColor', 'k');
hold on
s = surf(X, Y, CV_map_msh2);
set(s, 'FaceColor', 'flat', ...
       'CData', cmap_colors_msh2, ...
       'CDataMapping', 'direct', ...
       'EdgeColor', 'k');
xticks(xtick_values);
yticks(xtick_values);
xticklabels(arrayfun(@(v) sprintf('2^{%d}', log2(v)), xtick_values, 'UniformOutput', false));
yticklabels(arrayfun(@(v) sprintf('2^{%d}', log2(v)), xtick_values, 'UniformOutput', false));
xlabel('k_{chan}')
ylabel('k_{gj}')
zlabel('CV (cm/s)')
box off
set(gca,'FontSize',fsize,'TickDir','out','xscale','log','yscale','log')

%% figure chan/gj contrast + CV maps

%%% main figure
figure;
% fig_save.Position = [0 0 1000 1200];
t = tiledlayout(2,2);
t.TileSpacing = 'compact';
t.Padding = 'compact';
fsize = 18;
legend_fsize = 12;

load("mesh_data/FEMDATA_baseline.mat")
gj_area = FEM_data.gj_area_norm;
chan_area = FEM_data.Na_area_norm;

load("data/processed_data/CV_map7x7_msh1_baseline_cycle1000_beats10_D1_gj_loc_chan_loc")
CV_map_msh_base_D1 = CV_map;
load("data/processed_data/CV_map7x7_msh2_p60_ip60_cycle1000_beats10_D1_gj_loc_chan_loc")
CV_map_msh_6060_D1 = CV_map;

load("data/processed_data/CV_map7x7_msh1_baseline_cycle1000_beats10_D01_gj_loc_chan_loc")
CV_map_msh_base_D01 = CV_map;
load("data/processed_data/CV_map7x7_msh2_p60_ip60_cycle1000_beats10_D01_gj_loc_chan_loc")
CV_map_msh_6060_D01 = CV_map;
clear CV_map


nexttile()
plot((chan_area.^1)./(sum(chan_area.^1)),"LineWidth",2)
hold on
plot((chan_area.^0.0125)./(sum(chan_area.^0.0125)),"LineWidth",2)
plot((chan_area.^8)./(sum(chan_area.^8)),"LineWidth",2)
legend("k_{Na} = 1", "k_{Na} = 2^{-3}", "k_{Na} = 2^{3}",'location','northwest')
box off
set(gca,'FontSize',fsize,'TickDir','out');
xlabel("partition no")
ylabel("relative density of Na channels")

nexttile()
plot((gj_area.^1)./(sum(gj_area.^1)),"LineWidth",2)
hold on
plot((gj_area.^0.0125)./(sum(gj_area.^0.0125)),"LineWidth",2)
plot((gj_area.^8)./(sum(gj_area.^8)),"LineWidth",2)
legend("k_{gj} = 1", "k_{gj} = 2^{-3}", "k_{gj} = 2^{3}",'location','northwest')
box off
set(gca,'FontSize',fsize,'TickDir','out');
xlabel("partition no")
ylabel("relative density of GJs")


%%%%% CV maps 
scale_gj_loc_vec = [8 4 2 1 1/2 1/4 1/8];
scale_chan_loc_vec = [8 4 2 1 1/2 1/4 1/8];
[X,Y] = meshgrid(scale_gj_loc_vec,scale_chan_loc_vec);
xtick_values = flip(scale_gj_loc_vec); %

view1 = [110,7];
view2 = [150,10];

nexttile()
% figure with 2 colormaps
nColors = 256;
cmap_colors_msh1 = get_2dcolormap_value(CV_map_msh_base_D1(1:end-1,1:end-1), viridis(nColors));
cmap_colors_msh2 = get_2dcolormap_value(CV_map_msh_6060_D1(1:end-1,1:end-1), inferno(nColors));
s = surf(X, Y, CV_map_msh_base_D1);
set(s, 'FaceColor', 'flat', ...
       'CData', cmap_colors_msh1, ...
       'CDataMapping', 'direct', ...
       'EdgeColor', 'k');
hold on
s = surf(X, Y, CV_map_msh_6060_D1);
set(s, 'FaceColor', 'flat', ...
       'CData', cmap_colors_msh2, ...
       'CDataMapping', 'direct', ...
       'EdgeColor', 'k');
xticks(xtick_values);
yticks(xtick_values);
xticklabels(arrayfun(@(v) sprintf('2^{%d}', log2(v)), xtick_values, 'UniformOutput', false));
yticklabels(arrayfun(@(v) sprintf('2^{%d}', log2(v)), xtick_values, 'UniformOutput', false));
xlabel('k_{chan}')
ylabel('k_{gj}')
zlabel('CV (cm/s)')
view(view1)
box off
set(gca,'FontSize',fsize,'TickDir','out','xscale','log','yscale','log')


nexttile()
% figure with 2 colormaps
nColors = 256;
cmap_colors_msh1 = get_2dcolormap_value(CV_map_msh_base_D01(1:end-1,1:end-1), viridis(nColors));
cmap_colors_msh2 = get_2dcolormap_value(CV_map_msh_6060_D01(1:end-1,1:end-1), inferno(nColors));
s = surf(X, Y, CV_map_msh_base_D01);
set(s, 'FaceColor', 'flat', ...
       'CData', cmap_colors_msh1, ...
       'CDataMapping', 'direct', ...
       'EdgeColor', 'k');
hold on
s = surf(X, Y, CV_map_msh_6060_D01);
set(s, 'FaceColor', 'flat', ...
       'CData', cmap_colors_msh2, ...
       'CDataMapping', 'direct', ...
       'EdgeColor', 'k');
xticks(xtick_values);
yticks(xtick_values);
xticklabels(arrayfun(@(v) sprintf('2^{%d}', log2(v)), xtick_values, 'UniformOutput', false));
yticklabels(arrayfun(@(v) sprintf('2^{%d}', log2(v)), xtick_values, 'UniformOutput', false));
xlabel('k_{chan}')
ylabel('k_{gj}')
zlabel('CV (cm/s)')
view(view2)
box off
set(gca,'FontSize',fsize,'TickDir','out','xscale','log','yscale','log')

%% figure CV resti 

% load("data/processed_data/CV_resti_msh1_baseline_cycle_beats10_D01_gj_loc_chan_loc")
load("data/processed_data/CV_map7x7_msh1_baseline_cycle1000_beats10_D1_gj_loc_chan_loc")
% CV_cycles_all(beats_2_cycles_all == 0) = NaN;
figure
plot(cycle_vec,CV_cycles_all(5,:),'LineWidth', 2)
hold on
plot(cycle_vec,CV_cycles_all(1,:),'LineWidth', 2)
plot(cycle_vec,CV_cycles_all(7,:),'LineWidth', 2)
plot(cycle_vec,CV_cycles_all(3,:),'LineWidth', 2)
plot(cycle_vec,CV_cycles_all(9,:),'LineWidth', 2)

set(gca,'ColorOrderIndex',1)
load("data/processed_data/CV_resti_msh2_p60_ip60_cycle_beats10_D01_gj_loc_chan_loc")
% CV_cycles_all(beats_2_cycles_all == 0) = NaN;

plot(cycle_vec,CV_cycles_all(5,:),'--','LineWidth', 2)
plot(cycle_vec,CV_cycles_all(1,:),'--','LineWidth', 2)
plot(cycle_vec,CV_cycles_all(7,:),'--','LineWidth', 2)
plot(cycle_vec,CV_cycles_all(3,:),'--','LineWidth', 2)
plot(cycle_vec,CV_cycles_all(9,:),'--','LineWidth', 2)

legend("gj 1 / chan 1", "0.125 / 0.125", "0.125/8","8/0.125", "8/8", ...
       "gj 1 / chan 1", "0.125 / 0.125", "0.125/8","8/0.125", "8/8",'fontsize',14)


%% figure block DATA

load_file = "data/hetg_tissue_block/" + ...
     "msh1_baseline_msh2_p60_ip60_cycle220_beats25_D01_gj_loc1_1_chan_loc1_1";

%all except G_save - we don't need that for plots
load(load_file, 'bcl', 'D', 'scale_gj_loc', 'scale_chan_loc', ...
        'p', 'iEC', 'Nnodes', 'Ncell', 'Ncurrents', 'indices', 'Mdisc', ...
        'phi_axial_all', 'Na_cleft_mean_all', 'K_cleft_mean_all', 'Ca_cleft_mean_all', 'A_cleft_mean_all',...
        'Iind', 'ts', 'model', 'FEM_file_list', 'tissue_legend',...
        'tup', 'trepol', 'ts_save', 'Nint');

ts_25b = ts;
phi_axial_all_25b = phi_axial_all;
Na_cleft_mean_all_25b = Na_cleft_mean_all;
K_cleft_mean_all_25b = K_cleft_mean_all;
Ca_cleft_mean_all_25b = Ca_cleft_mean_all;
A_cleft_mean_all_25b = A_cleft_mean_all;

load_file = "data/hetg_tissue_block/" + ...
     "msh1_baseline_msh2_p60_ip60_cycle220_beats8_D01_gj_loc1_1_chan_loc1_1_last4";

%all except G_save - we don't need that for plots
load(load_file,'phi_save', 'I_save','S_save', 'bcl', 'D', 'scale_gj_loc', 'scale_chan_loc', ...
        'p', 'iEC', 'Nnodes', 'Ncell', 'Ncurrents', 'indices', 'Mdisc', ...
        'phi_axial_all', 'Na_cleft_mean_all', 'K_cleft_mean_all', 'Ca_cleft_mean_all', 'A_cleft_mean_all',...
        'Iind', 'ts', 'model', 'FEM_file_list', 'tissue_legend',...
        'tup', 'trepol', 'ts_save', 'Nint');

%current indices are according to Iind:
%first Ncell = axial
%one bulk, one bd
%100 pre/100 post etc

%phi indices are interior - membrane - cleft - membrane - interior etc

icleft = iEC(1:end-1);
iintra = setdiff(1:Nnodes-1,icleft);
phi_i = phi_save(iintra,:);
phi_cleft_all = phi_save(icleft,:);
phi_axial = phi_save(Iind(1:Ncell,1),:); 
phi_axial_norm = (phi_axial - min(phi_axial,[],'all')) ./ (max(phi_axial,[],'all') - min(phi_axial,[],'all'));

% collect Vm / ionic currents
Vm = phi_save(Iind(:,1),:) - phi_save(Iind(:,2),:);

INa_all = I_save(p.iina:Ncurrents:end,:);
ICa_all = I_save(p.iical:Ncurrents:end,:);

INa_axial = INa_all(1:Ncell,:);
ICa_axial = ICa_all(1:Ncell,:);

cleft1 = 25; 

ind_pre =  (Ncell+2 + (cleft1-1) * 2 * Mdisc )        : (Ncell+1 + (cleft1-1) * 2 * Mdisc + Mdisc);
ind_post = (Ncell+2 + (cleft1-1) * 2 * Mdisc + Mdisc) : (Ncell+1 + (cleft1-1) * 2 * Mdisc + 2*Mdisc);



Na_cleft_all = S_save(iEC,:);
K_cleft_all = S_save(iEC+Nnodes,:);
Ca_cleft_all = S_save(iEC+2*Nnodes,:);
A_cleft_all = S_save(iEC+3*Nnodes,:);

ion_cleft_index = (cleft1-1)*Mdisc+1:cleft1*Mdisc;

phi_cleft = phi_cleft_all((cleft1-1)*Mdisc+1:cleft1*Mdisc,:);
Na_cleft = Na_cleft_all((cleft1-1)*Mdisc+1:cleft1*Mdisc,:);
K_cleft = K_cleft_all((cleft1-1)*Mdisc+1:cleft1*Mdisc,:);
Ca_cleft = Ca_cleft_all((cleft1-1)*Mdisc+1:cleft1*Mdisc,:);
A_cleft = A_cleft_all((cleft1-1)*Mdisc+1:cleft1*Mdisc,:);


clear 'phi_save' 'I_save' 'S_save';

%% figure block FIGURE




%%% main figure
fig_save = figure;
% % fig_save.Position = [0 0 1000 1200];
t = tiledlayout(3,4);
t.TileSpacing = 'compact';
t.Padding = 'compact';
fsize = 20;
legend_fsize = 12;


x_time = ts_25b;
y_loc = zeros(1,length(x_time));

xlim_ap = [0,2195]; % first 10 beats out of 25
xlim_upst = [11,15];
view_coords = [355,75];
tissue_yticks = [1,10:10:50];
voltage_ticks = [-80,-40,0,40];
voltage_lim = [-90,50];



ylim_INa = [-3.75e-4,1.5e-4];
ylim_Na = [125,142];
ticks_Na = [135 140];

ylim_ICa = [-18e-7,1.5e-7];
ylim_Ca = [1.65,2.75];
ticks_Ca = [1.7 1.9 2.1];

ylim_phi = [-15,2];

xlim_upst1 = [0,15] + 1125; 
% xlim_upst1 = [0,15] + 1555; 
xlim_upst2 = [0,3] + 8.5;

phi_axial_norm_Ca = 1.7 + phi_axial_norm * (3.3 - 1.7);
phi_axial_norm_phi = -15 + phi_axial_norm * (10 + 15);

%%%% tile APs
nexttile([1,4])
for i = 1:49 
    hold on
    if i == 25 || i ==26
        plot3(x_time,y_loc+i,phi_axial_all_25b(i,:),'LineWidth', 2, 'color', [0 0 0])
    elseif i == 21 || i ==30
        plot3(x_time,y_loc+i,phi_axial_all_25b(i,:),'LineWidth', 2, 'color', [0.5 0.5 0.5])
    else
        plot3(x_time,y_loc+i,phi_axial_all_25b(i,:),'LineWidth', 2, 'color', [0.75 0.75 0.75])
    end
    % if i<50
    %     plot3(x_time,y_loc+i+0.5,mean(Ca_cleft_all((100*(i-1) + 1 : (100*i)),:)), ...
    %           'LineWidth', 2,'color','red')
    %     hold on
    % end
end
ax = gca;
% ax.DataAspectRatio = [1 0.5 1];

% ax.XRuler.FirstCrossoverValue  = 0; % X crossover with Y axis
% ax.XRuler.SecondCrossoverValue  = 0; % X crossover with Z axis
% ax.YRuler.FirstCrossoverValue  = -50; % Y crossover with X axis
% ax.YRuler.SecondCrossoverValue  = -50; % Y crossover with Z axis
% ax.ZRuler.FirstCrossoverValue  = 0; % Z crossover with X axis
% ax.ZRuler.SecondCrossoverValue = 0; % Z crossover with Y axis
xticks(0:200:xlim_ap(end))
yticks([1,10:10:50])
zticks([-80,-40,0,40])
xlim(xlim_ap)
ylim([1,50])
zlim([-90,50])
view(view_coords)
box off
set(gca,'FontSize',fsize,'TickDir','out');
xlabel("time (ms)")
ylabel("cell no")
zlabel("φ (mV)")

%%%%% TIME 1
x_time = ts_save;

%%% T1 Ca/Na conc
nexttile()
% plot(x_time,phi_axial(cleft1,:),'LineWidth', 2, 'color', [0 0 0])
% hold on
% plot(x_time,phi_axial(cleft1+1,:),'LineWidth', 2, 'color', [0 0 0])
% yticks(voltage_ticks)
% ylim(voltage_lim)
% ylabel("φ (mV)")

plot(x_time,Ca_cleft, '-','LineWidth', 2, 'color', [Ca_color, 0.25]);
hold on
plot(x_time,phi_axial_norm_Ca(cleft1,:),'LineWidth', 2, 'color', [0 0 0])
plot(x_time,phi_axial_norm_Ca(cleft1+1,:),'LineWidth', 2, 'color', [0 0 0])

ylim(ylim_Ca)
yticks(ticks_Ca)
ylabel("Ca^{2+}(mM)")


yyaxis right
plot(x_time,Na_cleft,'-','LineWidth', 2, 'color', [Na_color, 0.25]);
ylim(ylim_Na)
yticks(ticks_Na)
ylabel("Na^{+}(mM)")

xlim(xlim_upst1)
ax = gca;
box off
set(gca,'FontSize',fsize,'TickDir','out');
ax.YAxis(1).Color = [0 0 0];
ax.YAxis(2).Color = Na_color;
ax.XAxis.Visible = 'off';

hold on
leg_color = {};
leg_color{1} = plot(nan,'-','Linewidth', 5, 'color', [Na_color]);
leg_color{2} = plot(nan,'-','Linewidth', 5, 'color', [Ca_color]);
legend([leg_color{:}], {'Na^{+}_{cleft}','Ca^{2+}_{cleft}'}, 'location', 'best','box','off');



%%%% T1 phi cleft 
nexttile()
plot(x_time,phi_cleft, '-','LineWidth', 2, 'color', [phi_color, 0.25]);
hold on
plot(x_time,phi_axial_norm_phi(cleft1,:),'LineWidth', 2, 'color', [0 0 0])
plot(x_time,phi_axial_norm_phi(cleft1+1,:),'LineWidth', 2, 'color', [0 0 0])
ax = gca;
box off
set(gca,'FontSize',fsize,'TickDir','out');
ax.YAxis(1).Color = [0 0 0];
ax.XAxis.Visible = 'off';

ylim(ylim_phi)
ylabel("φ_{cleft} (mV)")
xlim(xlim_upst1)

hold on
leg_color = {};
leg_color{1} = plot(nan,'-','Linewidth', 5, 'color', [phi_color]);
legend([leg_color{:}], {'φ_{cleft}'}, 'location', 'best','box','off');

% %%%phi_cleft inset
% time_scale_bar_x = 13.5:0.1:14;
% time_scale_bar_y = zeros(length(time_scale_bar_x),1) -60;
% txt = string(time_scale_bar_x(end) - time_scale_bar_x(1)) + " ms";
% txt_x = (time_scale_bar_x(1)*0.945 + time_scale_bar_x(end))/2;
% txt_y = time_scale_bar_y(1)+8;
% 
% figure
% plot(x_time,phi_cleft,'LineWidth', 2, 'color', [phi_color 0.15])
% hold on
% plot(x_time,mean(phi_cleft,1),'LineWidth', 4, 'color', [phi_color]./1.5)
% plot(time_scale_bar_x,time_scale_bar_y,"LineWidth", 15, 'color', [0 0 0])
% text(txt_x,txt_y,txt,'FontSize',fsize+10)
% 
% xticks(xlim_upst(1):1:xlim_upst(end))
% xlim(xlim_upst)
% ylim(phi_lim);
% box off
% set(gca,'FontSize',fsize,'TickDir','out');
% ax = gca;
% ax.XAxis.Visible = 'off';
% ax.YAxis.Visible = 'off';



%Ca conc
% nexttile()
% plot(x_time,phi_axial(cleft1,:),'LineWidth', 2, 'color', [0 0 0])
% hold on
% plot(x_time,phi_axial(cleft1+1,:),'LineWidth', 2, 'color', [0 0 0])
% yticks(voltage_ticks)
% ylim(voltage_lim)
% ylabel("φ (mV)")
% 
% 
% yyaxis right
% plot(x_time,Ca_cleft_all((cleft1-1)*Mdisc+1:cleft1*Mdisc,:),...
%     '-','LineWidth', 2, 'color', [Ca_color, 0.25]);
% ylim(ylim_Ca)
% yticks(ticks_Ca)
% ylabel("Ca^{2+}(mM)")
% 
% xlim(xlim_upst_cleft2)
% ax = gca;
% box off
% set(gca,'FontSize',fsize,'TickDir','out');
% ax.YAxis(1).Color = [0 0 0];
% ax.YAxis(2).Color = Ca_color;
% ax.XAxis.Visible = 'off';
% 
% leg_color = {};
% leg_color{1} = plot(nan,'-','Linewidth', 5, 'color', [Ca_color]);
% legend([leg_color{:}], {'Ca^{2+}_{cleft}'}, 'location', 'best','box','off');

%% figure block clamp DATA



load_file = "data/clamp_exp_Ca/" + ...
     "msh1_baseline_msh2_p60_ip60_cycle220_beats10_D01_gj_loc1_1_chan_loc1_1";

%all except G_save - we don't need that for plots
load(load_file,'phi_save', 'I_save','S_save', 'bcl', 'D', 'scale_gj_loc', 'scale_chan_loc', ...
        'p', 'iEC', 'Nnodes', 'Ncell', 'Ncurrents', 'indices', 'Mdisc', ...
        'phi_axial_all', 'Na_cleft_mean_all', 'K_cleft_mean_all', 'Ca_cleft_mean_all', 'A_cleft_mean_all',...
        'Iind', 'ts', 'model', 'FEM_file_list', 'tissue_legend',...
        'tup', 'trepol', 'ts_save', 'Nint');

%current indices are according to Iind:
%first Ncell = axial
%one bulk, one bd
%100 pre/100 post etc

%phi indices are interior - membrane - cleft - membrane - interior etc

icleft = iEC(1:end-1);
iintra = setdiff(1:Nnodes-1,icleft);
phi_i = phi_save(iintra,:);
phi_cleft_all = phi_save(icleft,:);
phi_axial = phi_save(Iind(1:Ncell,1),:); 
% collect Vm / ionic currents
Vm = phi_save(Iind(:,1),:) - phi_save(Iind(:,2),:);

INa_all = I_save(p.iina:Ncurrents:end,:);
ICa_all = I_save(p.iical:Ncurrents:end,:);

INa_axial = INa_all(1:Ncell,:);
ICa_axial = ICa_all(1:Ncell,:);

cleft1 = 25; 

ind_pre =  (Ncell+2 + (cleft1-1) * 2 * Mdisc )        : (Ncell+1 + (cleft1-1) * 2 * Mdisc + Mdisc);
ind_post = (Ncell+2 + (cleft1-1) * 2 * Mdisc + Mdisc) : (Ncell+1 + (cleft1-1) * 2 * Mdisc + 2*Mdisc);



Na_cleft_all = S_save(iEC,:);
K_cleft_all = S_save(iEC+Nnodes,:);
Ca_cleft_all = S_save(iEC+2*Nnodes,:);
A_cleft_all = S_save(iEC+3*Nnodes,:);

ion_cleft_index = (cleft1-1)*Mdisc+1:cleft1*Mdisc;


clear 'phi_save' 'I_save' 'S_save';

%%
plot(ts,Na_cleft_mean_all(6,:))


%% Compute CV
thresh_activation = -10;

x = ts_save;
activ_time = [];
for i = 1:50
    data_iso_int = phi_axial(i,:);
    data_iso1 = data_iso_int(:,1:end-1);
    data_iso2 = data_iso_int(:,2:end);
    ind_th = find(data_iso1<thresh_activation & data_iso2> thresh_activation); % find indices below/above threshold
    y1 = data_iso1(ind_th); y2 = data_iso2(ind_th);
    m_slope = (y2-y1)./(x(ind_th+1)-x(ind_th));  % linear slope
    activ_time(i) = x(ind_th) - (data_iso_int(ind_th)-thresh_activation)./m_slope;
end

local_cv = (100./diff(activ_time))./10; %cm/s
local_cv = local_cv(5:45); %cm/s


figure
plot(5:45,local_cv, 'LineWidth',2)
box off
set(gca,'FontSize',fsize,'TickDir','out');
xlim([1,50])
xlabel('cell no')
ylabel("CV (cm/s)")
xticks([1,10:10:50])


%% Figure hetg tissue ALL IONS - SUPP?? 

%%% main figure
fig_save = figure;
% fig_save.Position = [0 0 1000 1200];
t = tiledlayout(5,6);
t.TileSpacing = 'compact';
t.Padding = 'compact';
fsize = 20;
legend_fsize = 12;

xlim_ap = [0,300];
xlim_upst = [11,15];
view_coords = [330,35];
tissue_yticks = [1,10:10:50];


x_time = ts_save - (ts_save(1))-45;
y_loc = zeros(1,length(x_time));

%cable voltage
nexttile([5,1])
for i = flip(1:50) 
    hold on
    if i == 25 || i ==26
        plot3(x_time,y_loc+i,phi_axial(i,:),'LineWidth', 2, 'color', [0 0 0])
    else
        plot3(x_time,y_loc+i,phi_axial(i,:),'LineWidth', 2, 'color', [0.5 0.5 0.5])
    end
    % if i<50
    %     plot3(x_time,y_loc+i+0.5,mean(Ca_cleft_all((100*(i-1) + 1 : (100*i)),:)), ...
    %           'LineWidth', 2,'color','red')
    %     hold on
    % end
end
ax = gca;
ax.DataAspectRatio = [1 0.03 0.8];

% ax.XRuler.FirstCrossoverValue  = 0; % X crossover with Y axis
% ax.XRuler.SecondCrossoverValue  = 0; % X crossover with Z axis
% ax.YRuler.FirstCrossoverValue  = -50; % Y crossover with X axis
% ax.YRuler.SecondCrossoverValue  = -50; % Y crossover with Z axis
% ax.ZRuler.FirstCrossoverValue  = 0; % Z crossover with X axis
% ax.ZRuler.SecondCrossoverValue = 0; % Z crossover with Y axis
xticks(0:100:400)
yticks(tissue_yticks)
zticks([-80,-40,0,40])
xlim(xlim_ap)
ylim([1,50])
zlim([-90,50])
view(view_coords)
box off
set(gca,'FontSize',fsize,'TickDir','out');


%cable ion concs
nexttile([5,1])
for i = flip(1:49) 
    Njunc = i;
    ion_cleft = phi_cleft_all((Njunc-1)*Mdisc+1:Njunc*Mdisc,:);
    hold on
    plot3(x_time,y_loc+i,mean(ion_cleft,1),'LineWidth', 2,'color',[[phi_color]./1.25, 0.5])
    % if i == 25 || i ==26
    %     plot3(x_time,y_loc+i,mean(ion_cleft,1),'LineWidth', 2,'color',[[phi_color]/1.5, 0.5])
    % else
    %     plot3(x_time,y_loc+i,mean(ion_cleft,1),'LineWidth', 2,'color',[[Ca_color]./1.25, 0.5])
    % end
    % if i<50
    %     plot3(x_time,y_loc+i+0.5,mean(Ca_cleft_all((100*(i-1) + 1 : (100*i)),:)), ...
    %           'LineWidth', 2,'color','red')
    %     hold on
    % end
end
ax = gca;
ax.DataAspectRatio = [1 0.03 0.05];

% ax.XRuler.FirstCrossoverValue  = 0; % X crossover with Y axis
% ax.XRuler.SecondCrossoverValue  = 0; % X crossover with Z axis
% ax.YRuler.FirstCrossoverValue  = -50; % Y crossover with X axis
% ax.YRuler.SecondCrossoverValue  = -50; % Y crossover with Z axis
% ax.ZRuler.FirstCrossoverValue  = 0; % Z crossover with X axis
% ax.ZRuler.SecondCrossoverValue = 0; % Z crossover with Y axis
xticks(0:100:400)
yticks(tissue_yticks)
% zticks([-80,-40,0,40])
xlim(xlim_ap)
ylim([1,50])
% zlim([-90,50])
view(view_coords)
box off
set(gca,'FontSize',fsize,'TickDir','out');
title("phi cleft")



%cable ion concs
nexttile([5,1])
for i = flip(1:49) 
    Njunc = i;
    ion_cleft = Na_cleft_all((Njunc-1)*Mdisc+1:Njunc*Mdisc,:);
    hold on
    plot3(x_time,y_loc+i,mean(ion_cleft,1),'LineWidth', 2,'color',[[Na_color]./1.25, 0.5])
    % if i == 25 || i ==26
    %     plot3(x_time,y_loc+i,mean(ion_cleft,1),'LineWidth', 2,'color',[[phi_color]/1.5, 0.5])
    % else
    %     plot3(x_time,y_loc+i,mean(ion_cleft,1),'LineWidth', 2,'color',[[Ca_color]./1.25, 0.5])
    % end
    % if i<50
    %     plot3(x_time,y_loc+i+0.5,mean(Ca_cleft_all((100*(i-1) + 1 : (100*i)),:)), ...
    %           'LineWidth', 2,'color','red')
    %     hold on
    % end
end
ax = gca;
ax.DataAspectRatio = [1 0.03 0.01];

% ax.XRuler.FirstCrossoverValue  = 0; % X crossover with Y axis
% ax.XRuler.SecondCrossoverValue  = 0; % X crossover with Z axis
% ax.YRuler.FirstCrossoverValue  = -50; % Y crossover with X axis
% ax.YRuler.SecondCrossoverValue  = -50; % Y crossover with Z axis
% ax.ZRuler.FirstCrossoverValue  = 0; % Z crossover with X axis
% ax.ZRuler.SecondCrossoverValue = 0; % Z crossover with Y axis
xticks(0:100:400)
yticks(tissue_yticks)
% zticks([-80,-40,0,40])
xlim(xlim_ap)
ylim([1,50])
% zlim([-90,50])
view(view_coords)
box off
set(gca,'FontSize',fsize,'TickDir','out');
title("Na cleft")


%cable ion concs
nexttile([5,1])
for i = flip(1:49) 
    Njunc = i;
    ion_cleft = Ca_cleft_all((Njunc-1)*Mdisc+1:Njunc*Mdisc,:);
    hold on
    plot3(x_time,y_loc+i,mean(ion_cleft,1),'LineWidth', 2,'color',[[Ca_color]./1.25, 0.5])
    % if i == 25 || i ==26
    %     plot3(x_time,y_loc+i,mean(ion_cleft,1),'LineWidth', 2,'color',[[phi_color]/1.5, 0.5])
    % else
    %     plot3(x_time,y_loc+i,mean(ion_cleft,1),'LineWidth', 2,'color',[[Ca_color]./1.25, 0.5])
    % end
    % if i<50
    %     plot3(x_time,y_loc+i+0.5,mean(Ca_cleft_all((100*(i-1) + 1 : (100*i)),:)), ...
    %           'LineWidth', 2,'color','red')
    %     hold on
    % end
end
ax = gca;
ax.DataAspectRatio = [1 0.03 0.001];

% ax.XRuler.FirstCrossoverValue  = 0; % X crossover with Y axis
% ax.XRuler.SecondCrossoverValue  = 0; % X crossover with Z axis
% ax.YRuler.FirstCrossoverValue  = -50; % Y crossover with X axis
% ax.YRuler.SecondCrossoverValue  = -50; % Y crossover with Z axis
% ax.ZRuler.FirstCrossoverValue  = 0; % Z crossover with X axis
% ax.ZRuler.SecondCrossoverValue = 0; % Z crossover with Y axis
xticks(0:100:400)
yticks(tissue_yticks)
% zticks([-80,-40,0,40])
xlim(xlim_ap)
ylim([1,50])
% zlim([-90,50])
view(view_coords)
box off
set(gca,'FontSize',fsize,'TickDir','out');
title("Ca cleft")



%cable ion concs
nexttile([5,1])
for i = flip(1:49) 
    Njunc = i;
    ion_cleft = K_cleft_all((Njunc-1)*Mdisc+1:Njunc*Mdisc,:);
    hold on
    plot3(x_time,y_loc+i,mean(ion_cleft,1),'LineWidth', 2,'color',[[K_color]./1.25, 0.5])
    % if i == 25 || i ==26
    %     plot3(x_time,y_loc+i,mean(ion_cleft,1),'LineWidth', 2,'color',[[phi_color]/1.5, 0.5])
    % else
    %     plot3(x_time,y_loc+i,mean(ion_cleft,1),'LineWidth', 2,'color',[[Ca_color]./1.25, 0.5])
    % end
    % if i<50
    %     plot3(x_time,y_loc+i+0.5,mean(Ca_cleft_all((100*(i-1) + 1 : (100*i)),:)), ...
    %           'LineWidth', 2,'color','red')
    %     hold on
    % end
end
ax = gca;
ax.DataAspectRatio = [1 0.03 0.001];

% ax.XRuler.FirstCrossoverValue  = 0; % X crossover with Y axis
% ax.XRuler.SecondCrossoverValue  = 0; % X crossover with Z axis
% ax.YRuler.FirstCrossoverValue  = -50; % Y crossover with X axis
% ax.YRuler.SecondCrossoverValue  = -50; % Y crossover with Z axis
% ax.ZRuler.FirstCrossoverValue  = 0; % Z crossover with X axis
% ax.ZRuler.SecondCrossoverValue = 0; % Z crossover with Y axis
xticks(0:100:400)
yticks(tissue_yticks)
% zticks([-80,-40,0,40])
xlim(xlim_ap)
ylim([1,50])
% zlim([-90,50])
view(view_coords)
box off
set(gca,'FontSize',fsize,'TickDir','out');
title("K cleft")





%cable ion concs
nexttile([5,1])
for i = flip(1:49) 
    Njunc = i;
    ion_cleft = A_cleft_all((Njunc-1)*Mdisc+1:Njunc*Mdisc,:);
    hold on
    plot3(x_time,y_loc+i,mean(ion_cleft,1),'LineWidth', 2,'color',[[A_color]./1.25, 0.5])
    % if i == 25 || i ==26
    %     plot3(x_time,y_loc+i,mean(ion_cleft,1),'LineWidth', 2,'color',[[phi_color]/1.5, 0.5])
    % else
    %     plot3(x_time,y_loc+i,mean(ion_cleft,1),'LineWidth', 2,'color',[[Ca_color]./1.25, 0.5])
    % end
    % if i<50
    %     plot3(x_time,y_loc+i+0.5,mean(Ca_cleft_all((100*(i-1) + 1 : (100*i)),:)), ...
    %           'LineWidth', 2,'color','red')
    %     hold on
    % end
end
ax = gca;
ax.DataAspectRatio = [1 0.03 0.005];

% ax.XRuler.FirstCrossoverValue  = 0; % X crossover with Y axis
% ax.XRuler.SecondCrossoverValue  = 0; % X crossover with Z axis
% ax.YRuler.FirstCrossoverValue  = -50; % Y crossover with X axis
% ax.YRuler.SecondCrossoverValue  = -50; % Y crossover with Z axis
% ax.ZRuler.FirstCrossoverValue  = 0; % Z crossover with X axis
% ax.ZRuler.SecondCrossoverValue = 0; % Z crossover with Y axis
xticks(0:100:400)
yticks(tissue_yticks)
% zticks([-80,-40,0,40])
xlim(xlim_ap)
ylim([1,50])
% zlim([-90,50])
view(view_coords)
box off
set(gca,'FontSize',fsize,'TickDir','out');

title("A cleft")



%% mesh plot


load('mesh_data/FEMDATA_baseline.mat')



A_plot = zeros(100,100);
A_plot(FEM_data.cleft_adjacency_matrix ~=0)  = 1;



centroid_list = FEM_data.partition_centers;

% sum(A_plot,1)
figure
G = graph(FEM_data.cleft_adjacency_matrix);
lwidths = G.Edges.Weight/max(FEM_data.cleft_adjacency_matrix, [],'all').*5;
plot(G,'XData',centroid_list(:,1),'YData',centroid_list(:,2),'ZData',centroid_list(:,3), ...
    'NodeColor', 'k','NodeLabel',[],'LineWidth', lwidths,'NodeLabel',[]);    
hold on
% plot3(centroid_list(286,1),centroid_list(286,2),centroid_list(286,3), '.', 'MarkerSize',24, 'color','red')
% plot3(centroid_list(357,1),centroid_list(357,2),centroid_list(357,3), '.', 'MarkerSize',24, 'color','red')
% plot3(centroid_list(467,1),centroid_list(467,2),centroid_list(467,3), '.', 'MarkerSize',24, 'color','red')
ax = gca;
ax.XTick = [];
ax.YTick = [];
ax.ZTick = [];
ax.DataAspectRatio = [1 1 1];


plot(FEM_data.Na_area_norm)%node 78

z_mesh = []; 
figure
view(3); 
plot(G,'XData',centroid_list(:,1),'YData',centroid_list(:,2),'ZData',centroid_list(:,3), ...
    'NodeColor', 'k','NodeLabel',[],'LineWidth', lwidths,'NodeLabel',[]);    
hold on
for i = 1:100
    mesh_3d = stlread("mesh_data/3d_mesh/baseline/partitions/parts_" + string(i) +  ".stl");
    z_mesh.vertices = mesh_3d.Points;
    z_mesh.faces = mesh_3d.ConnectivityList;
    drawMesh(z_mesh.vertices, z_mesh.faces,'FaceColor',rand(1,3) .* [0.5,0.5,0.5],'FaceAlpha',.25,'EdgeAlpha', 0)
end
ax = gca;
ax.DataAspectRatio = [1 1 1];

%% single part plot MESH DATA
% fem and mesh data
load('mesh_data/FEMDATA_baseline.mat')

for i = 1:100
    mesh_3d = stlread("mesh_data/3d_mesh/baseline/partitions/parts_" + string(i) +  ".stl");
    z(i).vertices = mesh_3d.Points;
    z(i).faces = mesh_3d.ConnectivityList;
end

mesh_3d = stlread("mesh_data/3d_mesh/baseline/chan_all/gj_all.stl");
gj_mesh.vertices = mesh_3d.Points;
gj_mesh.faces = mesh_3d.ConnectivityList;

mesh_3d = stlread("mesh_data/3d_mesh/baseline/chan_all/Nav_chan_interplicate.stl");
Nav_ip_mesh.vertices = mesh_3d.Points;
Nav_ip_mesh.faces = mesh_3d.ConnectivityList;
mesh_3d = stlread("mesh_data/3d_mesh/baseline/chan_all/Nav_chan_plicate.stl");
Nav_p_mesh.vertices = mesh_3d.Points;
Nav_p_mesh.faces = mesh_3d.ConnectivityList;

[Nav_mesh.vertices, Nav_mesh.faces] = ...
    concatenateMeshes(Nav_ip_mesh.vertices, Nav_ip_mesh.faces, Nav_p_mesh.vertices, Nav_p_mesh.faces);

mesh_3d = stlread("mesh_data/3d_mesh/baseline/chan_all/Kir21_chan_interplicate.stl");
Kir21_ip_mesh.vertices = mesh_3d.Points;
Kir21_ip_mesh.faces = mesh_3d.ConnectivityList;
mesh_3d = stlread("mesh_data/3d_mesh/baseline/chan_all/Kir21_chan_plicate.stl");
Kir21_p_mesh.vertices = mesh_3d.Points;
Kir21_p_mesh.faces = mesh_3d.ConnectivityList;
[Kir21_mesh.vertices, Kir21_mesh.faces] = ...
    concatenateMeshes(Kir21_ip_mesh.vertices, Kir21_ip_mesh.faces, Kir21_p_mesh.vertices, Kir21_p_mesh.faces);

mesh_3d = stlread("mesh_data/3d_mesh/baseline/chan_all/NKA_chan_interplicate.stl");
NKA_ip_mesh.vertices = mesh_3d.Points;
NKA_ip_mesh.faces = mesh_3d.ConnectivityList;
mesh_3d = stlread("mesh_data/3d_mesh/baseline/chan_all/NKA_chan_plicate.stl");
NKA_p_mesh.vertices = mesh_3d.Points;
NKA_p_mesh.faces = mesh_3d.ConnectivityList;
[NKA_mesh.vertices, NKA_mesh.faces] = ...
    concatenateMeshes(NKA_ip_mesh.vertices, NKA_ip_mesh.faces, NKA_p_mesh.vertices, NKA_p_mesh.faces);

% partition 6 is representative - use rand to find a good on
part_no = 6; %randi([1 100],1,1); %6
part_xlim = [5,7.5];
part_ylim = [9,11];
part_zlim = [-0.1,0.3];

%remove chan vertices outside figure limits for easier plotting
gj_mesh = remove_vertices_outside_figlim(gj_mesh,part_xlim,part_ylim,part_zlim);
Nav_mesh = remove_vertices_outside_figlim(Nav_mesh,part_xlim,part_ylim,part_zlim);
Kir21_mesh = remove_vertices_outside_figlim(Kir21_mesh,part_xlim,part_ylim,part_zlim);
NKA_mesh = remove_vertices_outside_figlim(NKA_mesh,part_xlim,part_ylim,part_zlim);


%% single part FIGURE
fsize = 14;
figure
view(3); 
drawMesh(z(part_no).vertices, z(part_no).faces, ...
         'FaceColor', [0.961, 0.961, 0.863],'FaceAlpha',0,'EdgeAlpha', 1)
hold on
drawMesh(Nav_mesh.vertices, Nav_mesh.faces, ...
         'FaceColor', [0.467, 0.867, 0.467],'FaceAlpha',.25,'EdgeAlpha', 0)
drawMesh(gj_mesh.vertices, gj_mesh.faces, ...
         'FaceColor', [1.0, 0.702, 0.729],'FaceAlpha',1,'EdgeAlpha', 0)
% drawMesh(Kir21_mesh.vertices, Kir21_mesh.faces, ...
%          'FaceColor', [0.992, 0.992, 0.588],'FaceAlpha',.75,'EdgeAlpha', 0)
% drawMesh(NKA_mesh.vertices, NKA_mesh.faces, ...
%          'FaceColor', [0.792, 0.749, 1.000],'FaceAlpha',.25,'EdgeAlpha', 0)
ax = gca;
ax.DataAspectRatio = [1 1 1];
xlim(part_xlim)
ylim(part_ylim)
zlim(part_zlim)
xlabel('x (μm)')
ylabel('y (μm)')
zlabel('z (μm)')
set(gca,'fontsize',fsize)





%% video - data
% fem and mesh data
load('mesh_data/FEMDATA_baseline.mat')

for i = 1:100
    mesh_3d = stlread("mesh_data/3d_mesh/p60_ip60/partitions/parts_" + string(i) +  ".stl");
    z(i).vertices = mesh_3d.Points;
    z(i).faces = mesh_3d.ConnectivityList;
end

% sim data
load_folder = "data/hetg_tissue/";
load_file = "msh1_baseline_cycle1000_beats10_D1_gj_loc1_8_chan_loc1_1";

load_path = load_folder + load_file;
%all except G_save - we don't need that for plots
load(load_path,'phi_save', 'I_save','S_save', 'bcl', 'D', 'scale_gj_loc', 'scale_chan_loc', ...
        'p', 'iEC', 'Nnodes', 'Ncell', 'Ncurrents', 'indices', 'Mdisc', ...
        'phi_axial_all', 'Na_cleft_mean_all', 'K_cleft_mean_all', 'Ca_cleft_mean_all', 'A_cleft_mean_all',...
        'Iind', 'ts', 'model', 'FEM_file_list', 'tissue_legend',...
        'tup', 'trepol', 'ts_save', 'Nint');

%current indices are according to Iind:
%first Ncell = axial
%one bulk, one bd
%100 pre/100 post etc

%phi indices are interior - membrane - cleft - membrane - interior etc

icleft = iEC(1:end-1);
iintra = setdiff(1:Nnodes-1,icleft);
phi_i = phi_save(iintra,:);
phi_cleft_all = phi_save(icleft,:);
phi_axial = phi_save(Iind(1:Ncell,1),:); 
% collect Vm / ionic currents
Vm = phi_save(Iind(:,1),:) - phi_save(Iind(:,2),:);

INa_all = I_save(p.iina:Ncurrents:end,:);
ICa_all = I_save(p.iical:Ncurrents:end,:);

INa_axial = INa_all(1:Ncell,:);
ICa_axial = ICa_all(1:Ncell,:);

Na_cleft_all = S_save(iEC,:);
K_cleft_all = S_save(iEC+Nnodes,:);
Ca_cleft_all = S_save(iEC+2*Nnodes,:);
A_cleft_all = S_save(iEC+3*Nnodes,:);

cleft = 25; 
ind_pre =  (Ncell+2 + (cleft-1) * 2 * Mdisc )        : (Ncell+1 + (cleft-1) * 2 * Mdisc + Mdisc);
ind_post = (Ncell+2 + (cleft-1) * 2 * Mdisc + Mdisc) : (Ncell+1 + (cleft-1) * 2 * Mdisc + 2*Mdisc);
ion_cleft_index = (cleft-1)*Mdisc+1:cleft*Mdisc;

phi_cleft = phi_cleft_all(ion_cleft_index,:);
Na_cleft = Na_cleft_all(ion_cleft_index,:);
K_cleft = K_cleft_all(ion_cleft_index,:);
Ca_cleft = Ca_cleft_all(ion_cleft_index,:);
A_cleft = A_cleft_all(ion_cleft_index,:);

INa_cleft_pre = INa_all(ind_pre,:);
ICa_cleft_pre = ICa_all(ind_pre,:);

INa_cleft_post = INa_all(ind_post,:);
ICa_cleft_post = ICa_all(ind_post,:);


clear 'phi_save' 'I_save' 'S_save';

%% video
I_Na_both = -(INa_cleft_pre + INa_cleft_post);


phi_min = min(phi_cleft,[],'all');
phi_max = max(phi_cleft,[],'all');
y_limits_phi = [phi_min-3, phi_max+3];

Na_min = min(Na_cleft,[],'all');
Na_max = max(Na_cleft,[],'all');
y_limits_Na = [Na_min-3, Na_max+3];

Ca_min = min(Ca_cleft,[],'all');
Ca_max = max(Ca_cleft,[],'all');
y_limits_Ca = [Ca_min-3, Ca_max+3];

K_min = min(K_cleft,[],'all');
K_max = max(K_cleft,[],'all');
y_limits_K = [K_min-3, K_max+3];

A_min = min(A_cleft,[],'all');
A_max = max(A_cleft,[],'all');
y_limits_A = [A_min-3, A_max+3];


I_Na_min = min(I_Na_both,[],'all');
I_Na_max = max(I_Na_both,[],'all');
y_limits_I_Na = [I_Na_min-3, I_Na_max+3];

cmap_viridis = viridis(256);


clear('F_vid')
k = 1;
vid_range = 250:450; %D1
% vid_range = 450:650;
time_plot_xlim = [ts_save(vid_range(1)), ts_save(vid_range(end))];
for i_time = vid_range

f = figure;
f.Position = [100 100 800 900];
t = tiledlayout(5,2);
t.TileSpacing = 'tight';
t.Padding = 'none';
disp(i_time)

%%%%% phi 
% phi cleft mesh

ax = nexttile();
for i_parts = 1:100
    cmap_val = get_colormap_value(phi_cleft(i_parts,i_time), phi_min, phi_max, cmap_viridis);
    drawMesh(z(i_parts).vertices, z(i_parts).faces,'FaceColor',cmap_val,'FaceAlpha',.5,'EdgeAlpha', 0)
    hold on
end

ax.XTick = [];
ax.YTick = [];
ax.ZTick = [];
ax.DataAspectRatio = [1 1 1];
title('φ cleft dynamics')

% phi over time
ax = nexttile();
plot(ts_save(vid_range),phi_cleft(:,vid_range),'Color',phi_color)
box off
xlim(time_plot_xlim)
ylim(y_limits_phi)
set(gca,'FontSize',14,'TickDir','out')
xline(ts_save(i_time))
ax.XAxis.Visible = 'off';


%%%%% Na
% Na cleft mesh

ax = nexttile();
for i_parts = 1:100
    cmap_val = get_colormap_value(Na_cleft(i_parts,i_time), Na_min, Na_max, cmap_viridis);
    drawMesh(z(i_parts).vertices, z(i_parts).faces,'FaceColor',cmap_val,'FaceAlpha',.5,'EdgeAlpha', 0)
    hold on
end
ax.XTick = [];
ax.YTick = [];
ax.ZTick = [];
ax.DataAspectRatio = [1 1 1];
title('Na^{+} cleft dynamics')

% Na over time
ax = nexttile();
plot(ts_save(vid_range),Na_cleft(:,vid_range),'Color',Na_color)
box off
xlim(time_plot_xlim)
ylim(y_limits_Na)
set(gca,'FontSize',14,'TickDir','out')
xline(ts_save(i_time))
ax.XAxis.Visible = 'off';

%%%%% Ca
% Ca cleft mesh

ax = nexttile();
for i_parts = 1:100
    cmap_val = get_colormap_value(Ca_cleft(i_parts,i_time), Ca_min, Ca_max, cmap_viridis);
    drawMesh(z(i_parts).vertices, z(i_parts).faces,'FaceColor',cmap_val,'FaceAlpha',.5,'EdgeAlpha', 0)
    hold on
end
ax.XTick = [];
ax.YTick = [];
ax.ZTick = [];
ax.DataAspectRatio = [1 1 1];
title('Ca^{2+} cleft dynamics')

% Ca over time
ax = nexttile();
plot(ts_save(vid_range),Ca_cleft(:,vid_range),'Color',Ca_color)
box off
xlim(time_plot_xlim)
ylim(y_limits_Ca)
set(gca,'FontSize',14,'TickDir','out')
xline(ts_save(i_time))
ax.XAxis.Visible = 'off';

%%%%% K
% K cleft mesh

ax = nexttile();
for i_parts = 1:100
    cmap_val = get_colormap_value(K_cleft(i_parts,i_time), K_min, K_max, cmap_viridis);
    drawMesh(z(i_parts).vertices, z(i_parts).faces,'FaceColor',cmap_val,'FaceAlpha',.5,'EdgeAlpha', 0)
    hold on
end
ax.XTick = [];
ax.YTick = [];
ax.ZTick = [];
ax.DataAspectRatio = [1 1 1];
title('K^{+} cleft dynamics')

% K over time
ax = nexttile();
plot(ts_save(vid_range),K_cleft(:,vid_range),'Color',K_color)
box off
xlim(time_plot_xlim)
ylim(y_limits_K)
set(gca,'FontSize',14,'TickDir','out')
xline(ts_save(i_time))
ax.XAxis.Visible = 'off';

%%%%% A
% A cleft mesh

ax = nexttile();
for i_parts = 1:100
    cmap_val = get_colormap_value(A_cleft(i_parts,i_time), A_min, A_max, cmap_viridis);
    drawMesh(z(i_parts).vertices, z(i_parts).faces,'FaceColor',cmap_val,'FaceAlpha',.5,'EdgeAlpha', 0)
    hold on
end
ax.XTick = [];
ax.YTick = [];
ax.ZTick = [];
ax.DataAspectRatio = [1 1 1];
title('A^{-} cleft dynamics')

% A over time
ax = nexttile();
plot(ts_save(vid_range),A_cleft(:,vid_range),'Color',A_color)
box off
xlim(time_plot_xlim)
ylim(y_limits_A)
set(gca,'FontSize',14,'TickDir','out')
xline(ts_save(i_time))
% ax.XAxis.Visible = 'off';


F_vid(k) = getframe(f);
k = k+1;
close(f)
end

vid_name = "vid_" + load_file + "_cleft_" + string(cleft);
q = VideoWriter(vid_name, 'MPEG-4');
q.Quality = 100;
open(q)
writeVideo(q,F_vid)
close(q)







%% MISC

scale_gj_loc_vec = [8 1];
scale_chan_loc_vec = [8 1];
cycle_vec = [200:5:300];

gj_chan_vec = combvec(scale_gj_loc_vec, scale_chan_loc_vec,cycle_vec);
gj_chan_vec(:,gj_chan_vec(1,:) == 1 & gj_chan_vec(2,:) == 1) = [];


%interpolated imagesc for corrcet time stretch
x_img = 0:x_time(end);
y_img = 1:50;
[X,Y] = meshgrid(x_time, y_img); %sample grid points
[Xq,Yq] = meshgrid(x_img, y_img);%querry grid points
phi_axial_interp = interp2(X,Y,phi_axial,Xq,Yq);

figure
imagesc(x_img,y_img,phi_axial_interp)
ax = gca;
ax.DataAspectRatio = [50 1 1];



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

%% old indexing 

% cleft = 25;
% ind_disc_pre = ((Ncell+2) : (Ncell+1+Mdisc)) + (cleft-1).*Mdisc;
% ind_disc_post = ((Ncell+2+Mdisc) : (Ncell+1+2*Mdisc)) + (cleft-1).*Mdisc;
% ICa_pre = ICa_all(ind_disc_pre,:);
% ICa_post = ICa_all(ind_disc_post,:);

% [~,ind] = sort(Iind(:,1));
% Iind_cable = Iind(ind,:);
% Vm_cable = Vm(ind,:);
% tup = tup(ind,:); trepol = trepol(ind,:);

% % sort axial / pre-/post-junctional membrane indices to match Vm order
% [ind_axial, ~] = ind2sub([length(ind) length(indices.ind_axial)], find(ind == indices.ind_axial));
% [ind_disc_pre, ~] = ind2sub([length(ind) length(indices.ind_disc_pre)], find(ind == indices.ind_disc_pre));
% [ind_disc_post, ~] = ind2sub([length(ind) length(indices.ind_disc_post)], find(ind == indices.ind_disc_post));
% 
% INa_axial = INa_all(ind_axial,:);
% INa_disc_pre = INa_all(ind_disc_pre,:);
% INa_disc_post = INa_all(ind_disc_post,:);
% 
% ICa_axial = ICa_all(ind_axial,:);
% ICa_disc_pre = ICa_all(ind_disc_pre,:);
% ICa_disc_post = ICa_all(ind_disc_post,:);
% 
% phi_i = phi_save(iintra,:);
% phi_axial = phi_i(ind_axial,:); 
% 


% %junction specific
% Njunc = 25;
% ind_post = ind_disc_post((Njunc-1)*Mdisc+2:Njunc*Mdisc+1); 
% Vm_post = Vm_cable(ind_post,:); 
% INa_post = INa_all(ind_post,:);
% ICa_post = ICa_all(ind_post,:);
% 
% ind_pre = ind_disc_pre((Njunc-1)*Mdisc+1:Njunc*Mdisc);
% Vm_pre = Vm_cable(ind_pre,:); 
% INa_pre = INa_all(ind_pre,:);
% ICa_pre = ICa_all(ind_pre,:);

%% functions

function cmap_val = get_colormap_value(data, data_min, data_max, cmap)
    norm_val = (data - data_min)./(data_max - data_min);
    color_index = int32(norm_val.*255)+1;
    cmap_val = cmap(color_index,:);
end

function cmap_val_2D = get_2dcolormap_value(data, cmap)
    [nColors,~] = size(cmap);
    % Normalize data to colormap indices
    data_norm = (data - min(data(:))) / (max(data(:)) - min(data(:)));
    color_idx = round(1 + (nColors - 1) * data_norm);  % Map to [1, nColors]
    color_idx = min(max(color_idx, 1), nColors);       % Ensure within bounds
    % Convert indices to RGB
    cmap_val_2D = reshape(cmap(color_idx, :), [size(data,1), size(data,2), 3]);
end

function z = remove_vertices_outside_figlim(z,xlim,ylim,zlim)

    lim_ind = zeros(size(z.vertices));
    lim_ind(:,1) = (z.vertices(:,1)<xlim(1) | z.vertices(:,1)>xlim(2));
    lim_ind(:,2) = (z.vertices(:,2)<ylim(1) | z.vertices(:,2)>ylim(2));
    lim_ind(:,3) = (z.vertices(:,3)<zlim(1) | z.vertices(:,3)>zlim(2));
    
    [z.vertices, z.faces] = removeMeshVertices(z.vertices, z.faces, sum(lim_ind,2)>0);

end