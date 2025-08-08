%% ====================================================================
%% R-SOP系统分布式电源和负荷数据处理
%% 对应原case33_R_SOP.m的动态数据部分
%% ====================================================================

run ieee_33_node_system_R_SOP.m                % 加载基础系统数据

%% 1. 光伏和负荷原始时变数据（保持与原DG_Load.m相同的时变特性）
%------------------------------------光伏数据------------------------------
Solar_origin_data=[113.778,113.778,113.778,113.778,113.778,113.778,...
                    113.778,123.667,157.333,233.333,320.556,366.778,...
                    375.778,366.667,320.667,260.111,160.111,133.222,...
                    113.778,113.778,113.778,113.778,113.778,113.778];

%------------------------------------负荷数据------------------------------
Load_origin_data=[139.778 135.333 129.222 133.222 149.222 197.778,...
                  224.556 204.889 193.556,173.889 174.667 173.556,...
                  169.222 168.667 169.222 183.556 196.111 219.333,...
                  224.778 239.778 226.222 190.667 165 131];

%% 2. 数据标准化处理
alpha=2.92;                                     % 光伏缩放因子
Solar_data=(Solar_origin_data-113.778)/(765.778-113.778)*4*alpha*1.45*7/7;
Load_data=(Load_origin_data-57.2222)/(241.222-57.2222)*6.5*1;

%% 3. 生成三相不平衡负荷矩阵
% 对应mpc.load_3ph_kW
load_3ph_kW = zeros(32, 7);                     % [节点号, Pa, Qa, Pb, Qb, Pc, Qc]
load_3ph_kW(:,1) = (2:33)';                     % 节点编号

for i = 1:32
    bus_idx = load_3ph_kW(i,1);                 % 当前节点号
    P_total_kW = load_data_scaled(i, 1);        % 节点总有功负荷
    Q_total_kvar = load_data_scaled(i, 2);      % 节点总无功负荷
    
    % 获取该节点的三相分配系数
    phase_factors = three_phase_imbalance(bus_idx, :);
    
    % 分配到三相
    load_3ph_kW(i, 2) = P_total_kW * phase_factors(1);   % A相有功
    load_3ph_kW(i, 3) = Q_total_kvar * phase_factors(1); % A相无功
    load_3ph_kW(i, 4) = P_total_kW * phase_factors(2);   % B相有功
    load_3ph_kW(i, 5) = Q_total_kvar * phase_factors(2); % B相无功
    load_3ph_kW(i, 6) = P_total_kW * phase_factors(3);   % C相有功
    load_3ph_kW(i, 7) = Q_total_kvar * phase_factors(3); % C相无功
end

%% 4. 光伏空间分布系数计算
% 基于pv_3ph_config生成分布系数
Solar_radio = zeros(33, 1);                     % 光伏分布系数矩阵
total_pv_capacity_MW = sum(pv_3ph_config(:,3)) / 1000; % 总光伏容量(MW)

% 按节点统计光伏容量并计算分布系数
unique_pv_nodes = unique(pv_3ph_config(:,1));   % 获取光伏节点
for i = 1:length(unique_pv_nodes)
    node = unique_pv_nodes(i);
    node_capacity = sum(pv_3ph_config(pv_3ph_config(:,1)==node, 3)); % 该节点总容量
    Solar_radio(node) = (node_capacity/1000) / total_pv_capacity_MW; % 分布系数
end

%% 5. 负荷空间分布系数
Load_radio = zeros(33, 1);                      % 负荷分布系数矩阵
p_load = Bus(:, 2) / 1000;                      % 各节点基础有功负荷(MW)
q_load = Bus(:, 3) / 1000;                      % 各节点基础无功负荷(MVar)

for a = 2 : 33                                  % 计算负荷分布系数
    Load_radio(a) = p_load(a) / sum(p_load);
end

%% 6. 时空耦合数据生成
N = 24;                                         % 24小时
p_Solar = zeros(33, N);                         % 各节点各时段光伏出力
p_Load = zeros(33, N);                          % 各节点各时段有功负荷
q_Load = zeros(33, N);                          % 各节点各时段无功负荷

for a = 1 : N
    p_Solar(:, a) = Solar_radio * Solar_data(a); % 光伏出力分布
    p_Load(:, a) = Load_radio * Load_data(a);   % 负荷分布
    q_Load(:, a) = (p_Load(:, a) ./ p_load) .* (q_load); % 无功负荷
end

q_Load(1, :) = 0;                               % 平衡节点无功为0
q_Load_data = sum(q_Load);                      % 系统总无功负荷

%% 7. R-SOP专用数据结构
% VSC详细规格定义
vsc_spec = struct('id', {}, 'capacity_kVA', {}, 'loss_coeffs', {});
for k = 1:r_sop_vsc_num
    vsc_spec(k).id = k;
    vsc_spec(k).capacity_kVA = r_sop_vsc_capacity;
    vsc_spec(k).loss_coeffs = r_sop_loss_coeffs;
end

% 储能系统结构
ess_struct = struct();
ess_struct.id = ess_id;
ess_struct.bus = ess_bus;
ess_struct.phase = [1, 1, 1];                   % 三相可控
ess_struct.capacity_kWh = ess_capacity_kWh;
ess_struct.power_kW = ess_power_kW;
ess_struct.soc_min = ess_soc_min;
ess_struct.soc_max = ess_soc_max;
ess_struct.soc_init = ess_soc_init;
ess_struct.eff_ch = ess_eff_ch;
ess_struct.eff_dis = ess_eff_dis;

% 负荷投切结构
ls_struct = struct();
ls_struct.id = ls_id;
ls_struct.bus = ls_bus;
ls_struct.priority = ls_priority;
ls_struct.max_shed_pct = ls_max_shed_pct;

%% 8. 兼容性变量生成
% 生成与stage1_switch_WithVSM.m兼容的变量
tan_theta = q_load ./ p_load;                   % 功率因数正切值

%% 9. 验证和可视化
total_P_kW_3ph = sum(sum(load_3ph_kW(:, [2,4,6]))); % 三相总有功负荷
total_Q_kvar_3ph = sum(sum(load_3ph_kW(:, [3,5,7]))); % 三相总无功负荷

fprintf('\n========== R-SOP系统时变数据验证 ==========\n');
fprintf('三相总有功负荷: %.1f kW\n', total_P_kW_3ph);
fprintf('三相总无功负荷: %.1f kVar\n', total_Q_kvar_3ph);
fprintf('光伏峰值出力: %.2f MW\n', max(Solar_data));
fprintf('负荷峰值需求: %.2f MW\n', max(Load_data));
fprintf('VSC数量: %d, 候选节点: %s\n', r_sop_vsc_num, mat2str(r_sop_candidate_nodes));

% 时变曲线绘制
figure('Name', 'R-SOP系统光伏出力和负荷需求');
x = 1:1:24;
plot(x, Solar_data, 'r-', 'LineWidth', 1.5);
hold on
plot(x, Load_data, 'b--', 'LineWidth', 1.5);
xlabel('Time(h)'); ylabel('Active power(MW)');
legend('PV output', 'Load demand');
title('R-SOP系统24小时光伏出力与负荷需求');
grid on;

fprintf('=========================================\n');