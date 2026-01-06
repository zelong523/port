%% 港区调度-配电协同 MILP（改进版：提升可视化差异 & 可行性）
% 方案0：规则基线（设备轨迹不移位 + 充电尽早完成）
% 方案1：削峰优先（允许设备轨迹小范围错峰移位 + 智能充电，min P_peak）
% 方案2：平滑优先（更大平滑权重，min TV=Σ|ΔP|，并兼顾峰值）
%
% 关键改动（与原稿对照）：
% 1) 两阶段目标：先最小化 slack_end（保障可行），再在该最小值下优化 Ppeak/TV/移位等（防止权重冲突）。
% 2) 移位惩罚适中：shiftCostWeight = 5e-3，既允许错峰又避免过度拖延导致“平移整个曲线”。
% 3) 充电斜率放宽：Ramp_ch = 60 kW/步，减少因可充窗短导致的锯齿与不可行。
% 4) 方案0 提前充电系数进一步减小：1e-4 -> 5e-5，使其“高峰且波动大”但不过度极端。
% 5) TV 乘步长归一化：TVh = Δt_h * Σ|ΔP|，权重与 kW 级目标可比；方案2 增大平滑权重。
% 6) 可充窗口比例默认 0.70（先确保可行性，再通过权重体现差异），若想更紧张可调低。
% 7) SHIFT_MAX 适度减小到 20（10min），进一步降低整数复杂度，提升可行/最优概率。
% 8) 增加 SOC_min 软约束与母线容量软约束，防止随机可充窗口导致中途掉到 SOC_min 之下或轻微越限时直接判不可行。
% 8) Pgrid 直接用表达式，清晰化模型。

clear; clc;
yalmip('clear');
rng(20251125);

%% ========== 1) 场景参数 ==========
Delta_t_min = 0.5;            % 30 s = 0.5 min
Delta_t_h   = Delta_t_min/60; % 小时
T = 720;                      % 6h / 0.5min = 720步
time_min = (0:T-1)*Delta_t_min;

Q = 3; R = 3;                 % 3台岸桥、3台场桥
E = 9; S = 4;                 % 9辆车、4个桩

P_bus_max  = 2600;            % kW
P_pile_max = 150;             % kW（单桩上限）
eta_ch     = 0.95;            % 充电效率
Ramp_ch    = 60;              % kW/步（进一步放宽，便于短窗补能）

E_bat   = 300;                % kWh
SOC_min = 0.20;
SOC_max = 0.90;
SOC_tar = 0.70;

P_aux = 200;                  % 站内其它固定负荷 kW

% 箱型比例
prob_type = [0.20 0.30 0.15 0.35]; % [20E,20F,40E,40F]

%% ========== 2) 分轴体功率常数（kW） ==========
% 岸桥(QC)
P_QC_H = [150 350 180 450];   % Hoist
P_QC_T = [ 80 100  90 120];   % Trolley
P_QC_G = 50;                  % Gantry

% 场桥(RTG)
P_RTG_H = [ 80 180 100 220];
P_RTG_T = [ 30  45  35  50];
P_RTG_G = 40;

%% ========== 3) 位置->等待/额外步数 ==========
shipRow_extraWait = [0 1 2];
shipRow_prob      = [0.50 0.30 0.20];

yard_extraWait = [0 2 4];
yard_prob      = [0.40 0.35 0.25];

%% ========== 4) 生成“分轴体功率轨迹” ==========
SHIFT_MAX = 20;              % 允许最大延后：20步=10min，降低整数复杂度
T_eff = T - SHIFT_MAX;

[stateQC_base, P_QC_base]   = genCraneTrace_QC(T_eff, prob_type, shipRow_prob, shipRow_extraWait, P_QC_H, P_QC_T, P_QC_G);
[stateRTG_base, P_RTG_base] = genCraneTrace_RTG(T_eff, prob_type, yard_prob, yard_extraWait, P_RTG_H, P_RTG_T, P_RTG_G);

P_QC0  = repmat([P_QC_base; zeros(SHIFT_MAX,1)], 1, Q);   % T x Q
P_RTG0 = repmat([P_RTG_base; zeros(SHIFT_MAX,1)], 1, R);  % T x R

stateQC0  = repmat([stateQC_base; zeros(SHIFT_MAX,1)], 1, Q);
stateRTG0 = repmat([stateRTG_base; zeros(SHIFT_MAX,1)], 1, R);

M_QC  = cell(Q,1);
M_RTG = cell(R,1);
for q=1:Q,  M_QC{q}  = buildShiftMatrix(P_QC0(:,q), SHIFT_MAX);  end
for r=1:R,  M_RTG{r} = buildShiftMatrix(P_RTG0(:,r), SHIFT_MAX); end

%% ========== 5) 车辆SOC与可充窗口 ==========
SOC0 = 0.25 + 0.55*rand(E,1);
d_step = 0.10 + 0.05*rand(E,1);

avail_ratio = 0.70;            % 更宽可充比例，提升可行性；需要更紧张再下调
A = buildAvailability(E, T, avail_ratio);

%% ========== 6) 三个方案求解 ==========
Pgrid_all = zeros(3,T);
Pch_all   = zeros(3,T);
shiftQC_all  = zeros(3,Q);
shiftRTG_all = zeros(3,R);
Ych_all   = cell(3,1);
stateQC_all  = cell(3,1);
stateRTG_all = cell(3,1);

metrics = struct('Ppeak',zeros(3,1),'VarP',zeros(3,1),'TV',zeros(3,1), ...
                 'Ech',zeros(3,1),'minSOC',zeros(3,1),'Cmax_min',zeros(3,1));

ops = sdpsettings('solver','cplex','verbose',1);
ops.cplex.timelimit = 600;                % 增加时间上限，避免早停
ops.cplex.mip.tolerances.mipgap = 0.10;   % 放宽 MIP gap，优先拿到可行解
ops.cplex.threads = 0;
% 可调的日志级别，便于排查
ops.cplex.display = 2;
ops.cplex.emphasis.mip = 1;               % 可行性优先

for scheme = 0:2

    %% --- 变量 ---
    Ns = SHIFT_MAX + 1;

    uQC  = binvar(Q, Ns, 'full');
    uRTG = binvar(R, Ns, 'full');

    ych  = binvar(E, T, 'full');
    pch  = sdpvar(E, T, 'full');
    soc  = sdpvar(E, T+1, 'full');
    slack_end = sdpvar(E,1);              % 末端 SOC 软约束
    slack_socmin = sdpvar(E, T+1, 'full'); % 过程 SOC_min 软约束

    slack_bus = sdpvar(T,1);              % 母线越限软约束

    Ppeak = sdpvar(1,1);
    dplus  = sdpvar(T,1);
    dminus = sdpvar(T,1);

    %% --- 表达式：设备功率 ---
    P_QC = zeros(T,1);
    for q=1:Q
        P_QC = P_QC + M_QC{q} * uQC(q,:)';
    end

    P_RTG = zeros(T,1);
    for r=1:R
        P_RTG = P_RTG + M_RTG{r} * uRTG(r,:)';
    end

    P_ch = sum(pch,1)';  % T x1
    Pgrid = P_aux + P_QC + P_RTG + P_ch;

    %% --- 约束 ---
    cons = [];

    cons = [cons, sum(uQC,2)==1, sum(uRTG,2)==1];

    if scheme==0
        cons = [cons, uQC(:,1)==1, uRTG(:,1)==1];
    end

    cons = [cons, 0 <= pch <= P_pile_max*ych];
    cons = [cons, ych <= A];
    cons = [cons, sum(ych,1) <= S];

    for t=2:T
        cons = [cons, pch(:,t)-pch(:,t-1) <= Ramp_ch];
        cons = [cons, pch(:,t-1)-pch(:,t) <= Ramp_ch];
    end

    cons = [cons, soc(:,1) == SOC0];
    for t=1:T
        cons = [cons, soc(:,t+1) == soc(:,t) + eta_ch*pch(:,t)*Delta_t_h/E_bat - d_step/E_bat];
    end
    cons = [cons, soc <= SOC_max];
    cons = [cons, soc + slack_socmin >= SOC_min, slack_socmin >= 0];

    cons = [cons, slack_end >= 0];
    cons = [cons, soc(:,T+1) + slack_end >= SOC_tar];

    cons = [cons, slack_bus >= 0];
    cons = [cons, Pgrid <= P_bus_max + slack_bus, Pgrid >= 0];

    cons = [cons, Ppeak >= Pgrid, Ppeak >= 0];

    cons = [cons, dplus >= 0, dminus >= 0];
    cons = [cons, dplus(1)==0, dminus(1)==0];
    for t=2:T
        cons = [cons, Pgrid(t)-Pgrid(t-1) == dplus(t)-dminus(t)];
    end
    TVh = Delta_t_h * sum(dplus + dminus); % 乘以步长的总变差

    %% --- 目标函数（两阶段：先最小化 slack，再优化峰值/平滑/移位） ---
    BIG = 1e4; % 基准权重
    BIG_bus = 5e3; % 母线越限惩罚

    shiftCost = 0;
    idx = (0:SHIFT_MAX)'; % 步数
    for q=1:Q, shiftCost = shiftCost + idx'*uQC(q,:)'; end
    for r=1:R, shiftCost = shiftCost + idx'*uRTG(r,:)'; end
    shiftCostWeight = 5e-3; % 适中惩罚，既能错峰又不会把曲线平移过多

    % 阶段1：最小化 slack_end 总和，保障可行
    obj_stage1 = sum(slack_end) + sum(slack_bus) + sum(slack_socmin(:));
    fprintf('\n=========== 求解方案 %d (阶段1：最小化 slack) ===========\n', scheme);
    sol = optimize(cons, obj_stage1, ops);
    if sol.problem ~= 0
        fprintf('方案 %d 阶段1 未达最优/失败：%s\n', scheme, yalmiperror(sol.problem));
        if scheme>0
            Pgrid_all(scheme+1,:) = Pgrid_all(1,:);
            Pch_all(scheme+1,:)   = Pch_all(1,:);
            continue;
        else
            error('方案0阶段1失败，基础可行性未满足。');
        end
    end
    best_slack = value(obj_stage1);
    % 将 slack 固定为最小值（允许极小数值容差）
    cons = [cons, obj_stage1 <= best_slack + 1e-6];

    % 阶段2：在 slack 最小的基础上优化各自的二级目标
    if scheme==0
        timeWeight = (1:T)'; % 越早越好（系数很小）
        obj = BIG*obj_stage1 + BIG_bus*sum(slack_bus) + 5e-5*sum(timeWeight.*P_ch) + shiftCostWeight*shiftCost;

    elseif scheme==1
        lambda_peak = 1.2;   % 加强削峰权重
        lambda_tv   = 0.04;  % 适度平滑
        obj = BIG*obj_stage1 + BIG_bus*sum(slack_bus) + lambda_peak*Ppeak + lambda_tv*TVh + shiftCostWeight*shiftCost;

    else
        lambda_peak = 0.25;  % 保持峰值关注但次于平滑
        lambda_tv   = 0.35;  % 平滑权重大，曲线更柔和
        obj = BIG*obj_stage1 + BIG_bus*sum(slack_bus) + lambda_peak*Ppeak + lambda_tv*TVh + shiftCostWeight*shiftCost;
    end

    fprintf('=========== 求解方案 %d (阶段2：多目标) ===========\n', scheme);
    sol = optimize(cons, obj, ops);

    if sol.problem ~= 0
        fprintf('方案 %d 阶段2 未达最优/失败：%s\n', scheme, yalmiperror(sol.problem));
        % 备选求解：放松 TV 和移位权重 + 更宽 gap/time
        ops_relax = ops;
        ops_relax.cplex.mip.tolerances.mipgap = 0.15;
        ops_relax.cplex.timelimit = 900;
        lambda_tv_relax = 0.2 * (scheme==2) * lambda_tv + 0.2 * (scheme==1) * lambda_tv; % 方案1/2都减弱
        if scheme==0
            obj_relax = BIG*obj_stage1 + BIG_bus*sum(slack_bus) + 5e-5*sum(timeWeight.*P_ch) + shiftCostWeight*0.5*shiftCost;
        elseif scheme==1
            obj_relax = BIG*obj_stage1 + BIG_bus*sum(slack_bus) + lambda_peak*Ppeak + lambda_tv_relax*TVh + shiftCostWeight*0.5*shiftCost;
        else
            obj_relax = BIG*obj_stage1 + BIG_bus*sum(slack_bus) + lambda_peak*Ppeak + lambda_tv_relax*TVh + shiftCostWeight*0.5*shiftCost;
        end
        fprintf('方案 %d 阶段2 重试（放松TV/移位，扩大gap/time）...\n', scheme);
        sol = optimize(cons, obj_relax, ops_relax);
        if sol.problem ~= 0
            fprintf('方案 %d 阶段2 仍失败：%s，回退方案0曲线。\n', scheme, yalmiperror(sol.problem));
            if scheme>0
                Pgrid_all(scheme+1,:) = Pgrid_all(1,:);
                Pch_all(scheme+1,:)   = Pch_all(1,:);
                continue;
            else
                error('方案0阶段2失败，无法生成基线。');
            end
        end
    end

    %% --- 提取结果 ---
    Pgrid_val = value(Pgrid)';
    pch_val   = value(pch);
    ych_val   = value(ych);
    soc_val   = value(soc);

    Pgrid_all(scheme+1,:) = Pgrid_val;
    Pch_all(scheme+1,:)   = sum(pch_val,1);
    Ych_all{scheme+1}     = ych_val;

    shiftQC = zeros(1,Q);
    for q=1:Q
        uu = value(uQC(q,:));
        [~,k] = max(uu);
        shiftQC(q) = k-1;
    end
    shiftRTG = zeros(1,R);
    for r=1:R
        uu = value(uRTG(r,:));
        [~,k] = max(uu);
        shiftRTG(r) = k-1;
    end
    shiftQC_all(scheme+1,:)  = shiftQC;
    shiftRTG_all(scheme+1,:) = shiftRTG;

    stateQC_shifted  = zeros(Q,T);
    stateRTG_shifted = zeros(R,T);
    for q=1:Q
        s = shiftQC(q);
        base = stateQC0(:,q)';
        stateQC_shifted(q,:) = [zeros(1,s), base(1:end-s)];
    end
    for r=1:R
        s = shiftRTG(r);
        base = stateRTG0(:,r)';
        stateRTG_shifted(r,:) = [zeros(1,s), base(1:end-s)];
    end
    stateQC_all{scheme+1}  = stateQC_shifted;
    stateRTG_all{scheme+1} = stateRTG_shifted;

    metrics.Ppeak(scheme+1)   = max(Pgrid_val);
    metrics.VarP(scheme+1)    = var(Pgrid_val);
    metrics.TV(scheme+1)      = sum(abs(diff(Pgrid_val)));
    metrics.Ech(scheme+1)     = sum(sum(pch_val))*Delta_t_h; % kWh
    metrics.minSOC(scheme+1)  = min(soc_val(:));

    lastIdx = find(P_QC0(:,1)'>0 | P_RTG0(:,1)'>0, 1, 'last');
    Cmax_step = min(T, lastIdx + max([shiftQC,shiftRTG]));
    metrics.Cmax_min(scheme+1) = (Cmax_step-1)*Delta_t_min;

end

%% ========== 7) 输出表（中文） ==========
schemeNames = {'方案0(规则基线)','方案1(削峰优先)','方案2(平滑优先)'};
fprintf('\n=== 宏观指标汇总（中文）===\n');
fprintf('%-14s %-10s %-12s %-10s %-10s %-12s %-10s\n', ...
    '方案','峰值P_peak','波动Var(P)','TV=Σ|ΔP|','完工(min)','累计充电量(kWh)','最低SOC');
for i=1:3
    fprintf('%-14s %-10.1f %-12.1f %-10.1f %-10.1f %-12.1f %-10.3f\n', ...
        schemeNames{i}, metrics.Ppeak(i), metrics.VarP(i), metrics.TV(i), ...
        metrics.Cmax_min(i), metrics.Ech(i), metrics.minSOC(i));
end

%% ========== 8) 图：总功率对比 ==========
figure('Name','系统总功率对比');
plot(time_min, Pgrid_all(1,:), 'LineWidth',1.2); hold on;
plot(time_min, Pgrid_all(2,:), 'LineWidth',1.2);
plot(time_min, Pgrid_all(3,:), 'LineWidth',1.2);
yline(P_bus_max,'--','P_{bus,max}','LineWidth',1.2);
xlabel('时间 / min'); ylabel('系统总功率 / kW');
title('不同方案下系统总功率曲线对比');
legend(schemeNames,'Location','best');
grid on;

%% ========== 9) 图：充电功率时序 ==========
figure('Name','充电功率对比');
plot(time_min, Pch_all(1,:), 'LineWidth',1.2); hold on;
plot(time_min, Pch_all(2,:), 'LineWidth',1.2);
plot(time_min, Pch_all(3,:), 'LineWidth',1.2);
xlabel('时间 / min'); ylabel('充电功率 / kW');
title('不同方案下充电功率时序对比');
legend(schemeNames,'Location','best');
grid on;

%% ========== 10) 图：|ΔP|对比 ==========
figure('Name','功率变化幅度对比');
for i=1:3
    dP = [0 abs(diff(Pgrid_all(i,:)))];
    plot(time_min, dP, 'LineWidth',1.0); hold on;
end
xlabel('时间 / min'); ylabel('|P(t)-P(t-1)| / kW');
title('不同方案下功率变化幅度对比（越小越平滑）');
legend(schemeNames,'Location','best');
grid on;

%% ========== 11) 甘特图：岸桥轴状态 ==========
plotGanttAxis(stateQC_all, '岸桥QC轴体作业甘特图', schemeNames, Delta_t_min);

%% ========== 12) 甘特图：场桥轴状态 ==========
plotGanttAxis(stateRTG_all, '场桥RTG轴体作业甘特图', schemeNames, Delta_t_min);

%% ========== 13) 甘特图：车辆充电(0/1) ==========
plotGanttCharge(Ych_all, '车辆充电甘特图', schemeNames, Delta_t_min);

disp('=== 运行结束：已输出表+曲线+甘特图 ===');

%% ========================= 本地函数 =========================
function [state, P] = genCraneTrace_QC(T_eff, prob_type, rowProb, rowExtraWait, P_H, P_T, P_G)
    state = zeros(T_eff,1); P = zeros(T_eff,1);
    t = 1;
    while t <= T_eff
        c = drawDiscrete(prob_type);
        g = drawDiscrete(rowProb);
        extraIdle = rowExtraWait(g);
        seq = [ones(1,3), 2*ones(1,2), 3, zeros(1,1+extraIdle)]; % HHH TT G idle(+)
        for s = seq
            if t > T_eff, break; end
            state(t) = s;
            if s==1, P(t)=P_H(c);
            elseif s==2, P(t)=P_T(c);
            elseif s==3, P(t)=P_G;
            else, P(t)=0;
            end
            t = t + 1;
        end
    end
end

function [state, P] = genCraneTrace_RTG(T_eff, prob_type, yardProb, yardExtraWait, P_H, P_T, P_G)
    state = zeros(T_eff,1); P = zeros(T_eff,1);
    t = 1;
    while t <= T_eff
        c = drawDiscrete(prob_type);
        h = drawDiscrete(yardProb);
        extraIdle = yardExtraWait(h);
        seq = [ones(1,2), 2, 3, zeros(1,1+extraIdle)]; % HH T G idle(+)
        for s = seq
            if t > T_eff, break; end
            state(t) = s;
            if s==1, P(t)=P_H(c);
            elseif s==2, P(t)=P_T(c);
            elseif s==3, P(t)=P_G;
            else, P(t)=0;
            end
            t = t + 1;
        end
    end
end

function M = buildShiftMatrix(baseVec, SHIFT_MAX)
    T = length(baseVec);
    Ns = SHIFT_MAX + 1;
    M = zeros(T, Ns);
    for s = 0:SHIFT_MAX
        M(:,s+1) = [zeros(s,1); baseVec(1:T-s)];
    end
end

function A = buildAvailability(E, T, ratio)
    A = zeros(E,T);
    for e=1:E
        pos = 1;
        while pos <= T
            cycle = randi([16, 28]);                % 8~14 min
            idle  = max(4, round(ratio*cycle));     % 可充窗口长度
            start = pos + randi([0, max(0,cycle-idle)]);
            A(e, start:min(T,start+idle-1)) = 1;
            pos = pos + cycle;
        end
    end
end

function idx = drawDiscrete(prob)
    r = rand;
    cs = cumsum(prob(:));
    idx = find(r <= cs, 1, 'first');
    if isempty(idx), idx = length(prob); end
end

function plotGanttAxis(state_all, figTitle, schemeNames, dt_min)
    cmap = [1 1 1; 0.90 0.40 0.40; 0.40 0.60 0.90; 0.40 0.80 0.50]; % 白/红/蓝/绿
    for i=1:3
        figure('Name',[figTitle,'-',schemeNames{i}]);
        S = state_all{i};
        imagesc((0:size(S,2)-1)*dt_min, 1:size(S,1), S);
        colormap(cmap);
        caxis([0 3]);
        xlabel('时间 / min'); ylabel('设备编号');
        title([figTitle,' - ',schemeNames{i}]);
        grid on;
    end
end

function plotGanttCharge(Ych_all, figTitle, schemeNames, dt_min)
    cmap = [0.95 0.95 0.95; 0.20 0.20 0.20];
    for i=1:3
        figure('Name',[figTitle,'-',schemeNames{i}]);
        Y = Ych_all{i};
        imagesc((0:size(Y,2)-1)*dt_min, 1:size(Y,1), Y);
        colormap(cmap);
        caxis([0 1]);
        xlabel('时间 / min'); ylabel('车辆编号');
        title([figTitle,' - ',schemeNames{i}]);
        grid on;
    end
end
