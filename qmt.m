%% ========================================================================
%  Quantum Tunneling & Tunnel Diode — Extended MATLAB Analysis
%  ESE 4310 Final Project
%  ========================================================================
%  This script extends the baseline Python simulation with three deeper
%  investigations that go beyond the original code:
%
%   Part 1 — Analytical Validation
%            Compare the transfer-matrix solver against the exact
%            closed-form transmission formula for a rectangular barrier.
%            Proves numerical accuracy before building on the method.
%
%   Part 2 — Realistic GaAs/AlGaAs Resonant Tunneling Diode
%            Use real effective masses, a physically tilted potential
%            under applied bias, and Fermi–Dirac statistics to compute
%            temperature-dependent I-V curves with NDR.
%
%   Part 3 — Superlattice Miniband Formation
%            Extend from 2 barriers to N = 3, 5, 10 and show how
%            isolated resonant peaks broaden into continuous minibands.
%
%  All figures are saved to ./matlab_outputs/
%  ========================================================================

clear; clc; close all;

%% ======================== Physical Constants ============================
m0   = 9.1093837015e-31;   % free electron mass (kg)
hbar = 1.054571817e-34;    % reduced Planck constant (J·s)
h    = 6.62607015e-34;     % Planck constant (J·s)
q    = 1.602176634e-19;    % elementary charge (C)
kB   = 1.380649e-23;       % Boltzmann constant (J/K)
nm   = 1e-9;               % metres per nanometre

%% ======================== Output Directory ==============================
out_dir = 'matlab_outputs';
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

%% ========================================================================
%       PART 1 — ANALYTICAL VALIDATION OF THE TRANSFER-MATRIX METHOD
%  ========================================================================
%  The exact transmission coefficient for a single rectangular barrier
%  of height V0 and width a is (Griffiths, Eq. 2.169):
%
%    E < V0:  T = [1 + V0^2 sinh^2(kappa*a) / (4*E*(V0-E))]^(-1)
%             where kappa = sqrt(2*m*(V0-E)) / hbar
%
%    E > V0:  T = [1 + V0^2 sin^2(k'*a) / (4*E*(E-V0))]^(-1)
%             where k' = sqrt(2*m*(E-V0)) / hbar
%
%  We compare this against our transfer-matrix numerical result.
%  ========================================================================
fprintf('\n=== PART 1: Analytical Validation ===\n');

V0_val   = 0.35;                       % barrier height (eV)
a_val    = 1.1;                        % barrier width (nm)
mass_val = m0;                         % free electron mass
x_val    = linspace(-3.0, 3.0, 500);   % spatial grid (nm)
V_val    = rectangular_barrier(x_val, V0_val, a_val);

energies_val = linspace(0.005, 0.8, 300);  % energy sweep (eV)

% ---- Numerical: transfer-matrix method ---------------------------------
T_numerical = zeros(size(energies_val));
for i = 1:numel(energies_val)
    [T_numerical(i), ~, ~] = solve_scattering(x_val, V_val, ...
        energies_val(i), mass_val, q, hbar, nm);
end

% ---- Analytical: exact closed-form formula ------------------------------
T_analytical = zeros(size(energies_val));
a_m = a_val * nm;   % barrier width in metres
for i = 1:numel(energies_val)
    E = energies_val(i);
    if abs(E - V0_val) < 1e-10
        % At E = V0, use the limit: T = 1/(1 + m*V0*a^2/(2*hbar^2))
        T_analytical(i) = 1 / (1 + mass_val * V0_val*q * a_m^2 ...
                           / (2 * hbar^2));
    elseif E < V0_val
        kappa = sqrt(2 * mass_val * (V0_val - E) * q) / hbar;
        T_analytical(i) = 1 / (1 + V0_val^2 * sinh(kappa * a_m)^2 ...
                           / (4 * E * (V0_val - E)));
    else
        kp = sqrt(2 * mass_val * (E - V0_val) * q) / hbar;
        T_analytical(i) = 1 / (1 + V0_val^2 * sin(kp * a_m)^2 ...
                           / (4 * E * (E - V0_val)));
    end
end

% ---- Plot: overlay comparison -------------------------------------------
figure('Position', [80 80 800 500]);
semilogy(energies_val, T_analytical, 'b-', 'LineWidth', 2.5, ...
    'DisplayName', 'Analytical (exact)');
hold on;
semilogy(energies_val, T_numerical, 'r--', 'LineWidth', 1.8, ...
    'DisplayName', 'Transfer matrix (numerical)');
xline(V0_val, ':', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.3, ...
    'DisplayName', sprintf('V_0 = %.2f eV', V0_val));
xlabel('Electron Energy (eV)');
ylabel('Transmission T(E)');
title('Validation: Analytical vs Transfer-Matrix Method');
legend('Location', 'southeast');
grid on; set(gca, 'GridAlpha', 0.25);
annotation('textbox', [0.15 0.78 0.3 0.08], 'String', ...
    sprintf('Barrier: V_0 = %.2f eV, a = %.1f nm\nMass: m_0 (free electron)', ...
    V0_val, a_val), 'EdgeColor', 'none', 'FontSize', 9);
exportgraphics(gcf, fullfile(out_dir, 'part1_validation_overlay.png'), ...
    'Resolution', 200);
close;

% ---- Plot: absolute and relative error ----------------------------------
abs_error = abs(T_numerical - T_analytical);
rel_error = abs_error ./ max(T_analytical, 1e-30) * 100;   % percent

figure('Position', [80 80 800 600]);
subplot(2,1,1);
semilogy(energies_val, abs_error, 'Color', [0.64 0.08 0.18], 'LineWidth', 1.6);
ylabel('Absolute Error |T_{num} - T_{exact}|');
title('Transfer-Matrix Error Analysis');
grid on; set(gca, 'GridAlpha', 0.25);
xline(V0_val, ':', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.2);

subplot(2,1,2);
semilogy(energies_val, rel_error, 'Color', [0.0 0.35 0.65], 'LineWidth', 1.6);
xlabel('Electron Energy (eV)');
ylabel('Relative Error (%)');
grid on; set(gca, 'GridAlpha', 0.25);
xline(V0_val, ':', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.2);
exportgraphics(gcf, fullfile(out_dir, 'part1_validation_error.png'), ...
    'Resolution', 200);
close;

fprintf('  Max absolute error: %.2e\n', max(abs_error));
fprintf('  Mean relative error: %.4f %%\n', mean(rel_error(T_analytical > 1e-12)));

%% ========================================================================
%       PART 2 — REALISTIC GaAs/AlGaAs RESONANT TUNNELING DIODE
%  ========================================================================
%  Key improvements over the baseline Python model:
%    • GaAs effective mass  m* = 0.067 m0  (well / contacts)
%    • Al_{0.3}Ga_{0.7}As  m* = 0.092 m0  (barrier)
%    • Barrier height ≈ 0.23 eV for 30% Al fraction
%    • Under applied bias V_a, the potential tilts linearly across
%      the entire structure (not just a rigid energy shift)
%    • Current is computed from the Tsu-Esaki formula with
%      Fermi-Dirac occupation at finite temperature
%  ========================================================================
fprintf('\n=== PART 2: GaAs/AlGaAs Resonant Tunneling Diode ===\n');

% ---- Material parameters ------------------------------------------------
m_well    = 0.067 * m0;   % GaAs effective mass
m_barrier = 0.092 * m0;   % Al_{0.3}Ga_{0.7}As effective mass
V0_rtd    = 0.23;         % conduction band offset (eV)
b_width   = 3.0;          % barrier width (nm) — typical RTD
w_width   = 5.0;          % well width (nm)
L_contact = 8.0;          % flat contact region each side (nm)

% Total structure length
L_total = 2*L_contact + 2*b_width + w_width;
N_pts   = 800;           % fine grid for T(E) and wavefunction plots
x_rtd   = linspace(-L_total/2, L_total/2, N_pts);
N_coarse = 300;          % coarser grid for the I-V loop (speed)
x_coarse = linspace(-L_total/2, L_total/2, N_coarse);

% ---- Build effective mass profile (position-dependent) ------------------
m_eff = m_well * ones(size(x_rtd));
hw = w_width / 2;
left_b  = (x_rtd >= -hw - b_width) & (x_rtd <= -hw);
right_b = (x_rtd >=  hw) & (x_rtd <=  hw + b_width);
m_eff(left_b)  = m_barrier;
m_eff(right_b) = m_barrier;

% Coarse mass profile (for I-V computation)
m_eff_c = m_well * ones(size(x_coarse));
left_bc  = (x_coarse >= -hw - b_width) & (x_coarse <= -hw);
right_bc = (x_coarse >=  hw) & (x_coarse <=  hw + b_width);
m_eff_c(left_bc)  = m_barrier;
m_eff_c(right_bc) = m_barrier;

% ---- 2a. Zero-bias transmission spectrum --------------------------------
fprintf('  Computing zero-bias T(E)...\n');
V_rtd_0 = double_barrier(x_rtd, V0_rtd, b_width, w_width);
E_rtd   = linspace(0.001, 0.35, 400);
T_rtd_0 = zeros(size(E_rtd));
for i = 1:numel(E_rtd)
    [T_rtd_0(i), ~, ~] = solve_scattering_meff(x_rtd, V_rtd_0, ...
        E_rtd(i), m_eff, q, hbar, nm);
end

% Find resonant peak positions
[pks, locs] = findpeaks(T_rtd_0, 'MinPeakHeight', 0.05);
E_res_peaks = E_rtd(locs);
fprintf('  Resonant levels found at: ');
fprintf('%.4f eV  ', E_res_peaks); fprintf('\n');

% Plot zero-bias T(E) with resonance markers
figure('Position', [80 80 800 500]);
semilogy(E_rtd*1000, T_rtd_0, 'Color', [0.0 0.44 0.75], 'LineWidth', 2.2);
hold on;
for pk = 1:numel(locs)
    semilogy(E_rtd(locs(pk))*1000, pks(pk), 'rv', 'MarkerSize', 10, ...
        'MarkerFaceColor', [0.9 0.2 0.2], 'HandleVisibility', 'off');
    text(E_rtd(locs(pk))*1000 + 5, pks(pk), ...
        sprintf('E_%d = %.1f meV', pk, E_res_peaks(pk)*1000), ...
        'FontSize', 9, 'Color', [0.7 0.1 0.1]);
end
xlabel('Electron Energy (meV)');
ylabel('Transmission T(E)');
title('GaAs/AlGaAs RTD — Zero-Bias Transmission');
grid on; set(gca, 'GridAlpha', 0.25);
annotation('textbox', [0.15 0.15 0.35 0.15], 'String', ...
    {sprintf('Barrier: Al_{0.3}Ga_{0.7}As, V_0 = %.0f meV', V0_rtd*1000), ...
     sprintf('Barrier width = %.1f nm', b_width), ...
     sprintf('Well width = %.1f nm', w_width), ...
     sprintf('m^*_{well} = 0.067 m_0,  m^*_{bar} = 0.092 m_0')}, ...
    'EdgeColor', [0.5 0.5 0.5], 'BackgroundColor', [1 1 1 0.85], ...
    'FontSize', 8.5);
exportgraphics(gcf, fullfile(out_dir, 'part2_rtd_transmission.png'), ...
    'Resolution', 200);
close;

% ---- 2b. Wavefunction at first resonance --------------------------------
fprintf('  Computing wavefunction at first resonance...\n');
if ~isempty(E_res_peaks)
    E_wave = E_res_peaks(1);
else
    E_wave = 0.05;  % fallback
end
[~, ~, psi_rtd] = solve_scattering_meff(x_rtd, V_rtd_0, E_wave, ...
    m_eff, q, hbar, nm);
prob_density = abs(psi_rtd).^2;
prob_density = prob_density / max(prob_density);

figure('Position', [80 80 800 500]);
yyaxis left;
plot(x_rtd, prob_density, 'Color', [0.0 0.44 0.75], 'LineWidth', 2.0);
ylabel('|\psi(x)|^2 (normalized)');
yyaxis right;
plot(x_rtd, V_rtd_0 * 1000, 'Color', [0.8 0.4 0.0], 'LineWidth', 1.8);
ylabel('Potential (meV)');
xlabel('Position (nm)');
title(sprintf('RTD Wavefunction at First Resonance E = %.1f meV', E_wave*1000));
grid on; set(gca, 'GridAlpha', 0.25);
exportgraphics(gcf, fullfile(out_dir, 'part2_rtd_wavefunction.png'), ...
    'Resolution', 200);
close;

% ---- 2c. Tilted potential under bias & I-V curve -----------------------
fprintf('  Computing I-V curves with tilted potential...\n');

% Tsu–Esaki current formula:
%   J = (q*m*kT)/(2*pi^2*hbar^3) * integral T(E,V) * 
%       [ln(1+exp((Ef-E)/kT)) - ln(1+exp((Ef-E-qV)/kT))] dE
%
% We evaluate this at multiple temperatures.

Ef       = 0.01;              % Fermi level above conduction band (eV)
temps    = [77, 150, 300];    % temperatures (K)
V_bias   = linspace(0, 0.50, 60);   % reduced from 180
J_all    = zeros(numel(temps), numel(V_bias));

for ti = 1:numel(temps)
    T_K = temps(ti);
    kT  = kB * T_K / q;   % in eV
    fprintf('    T = %d K: ', T_K);
    
    for vi = 1:numel(V_bias)
        Va = V_bias(vi);
        
        % Build tilted potential on COARSE grid for speed
        V_tilted = double_barrier(x_coarse, V0_rtd, b_width, w_width);
        V_tilted = V_tilted - Va * (x_coarse - x_coarse(1)) / (x_coarse(end) - x_coarse(1));
        
        % Compute T(E) at this bias (coarse grid + fewer energies)
        E_scan = linspace(0.001, 0.30, 80);
        T_bias = zeros(size(E_scan));
        for ei = 1:numel(E_scan)
            [T_bias(ei), ~, ~] = solve_scattering_meff(x_coarse, V_tilted, ...
                E_scan(ei), m_eff_c, q, hbar, nm);
        end
        
        % Tsu-Esaki integrand with Fermi-Dirac supply function
        f_left  = log(1 + exp((Ef - E_scan) / kT));
        f_right = log(1 + exp((Ef - E_scan - Va) / kT));
        supply  = f_left - f_right;
        
        integrand = T_bias .* supply;
        J_all(ti, vi) = q * m_well * kT / (2 * pi^2 * hbar^3) ...
                         * trapz(E_scan * q, integrand);
        
        % Progress indicator
        if mod(vi, 15) == 0
            fprintf('%d%% ', round(vi/numel(V_bias)*100));
        end
    end
    fprintf('done\n');
end

% Normalize for plotting (current density in A/m^2 — large numbers)
figure('Position', [80 80 800 550]);
colors_t = [0.0 0.30 0.70;  0.0 0.60 0.35;  0.85 0.15 0.15];
for ti = 1:numel(temps)
    plot(V_bias * 1000, J_all(ti,:), 'Color', colors_t(ti,:), ...
        'LineWidth', 2.2, 'DisplayName', sprintf('T = %d K', temps(ti)));
    hold on;
end
xlabel('Applied Bias (mV)');
ylabel('Current Density J (A/m^2)');
title('GaAs/AlGaAs RTD I-V Characteristics — Temperature Dependence');
legend('Location', 'northwest');
grid on; set(gca, 'GridAlpha', 0.25);
annotation('textbox', [0.55 0.70 0.35 0.12], 'String', ...
    {'Tilted potential under bias', ...
     'Tsu-Esaki formula', ...
     'Fermi-Dirac statistics'}, ...
    'EdgeColor', [0.5 0.5 0.5], 'BackgroundColor', [1 1 1 0.85], ...
    'FontSize', 9, 'FontAngle', 'italic');
exportgraphics(gcf, fullfile(out_dir, 'part2_rtd_iv_temperature.png'), ...
    'Resolution', 200);
close;

% ---- 2d. Potential profile visualisation at several biases --------------
figure('Position', [80 80 800 450]);
bias_show = [0, 0.10, 0.20, 0.35];
cmap = lines(numel(bias_show));
for bi = 1:numel(bias_show)
    Va = bias_show(bi);
    V_show = double_barrier(x_rtd, V0_rtd, b_width, w_width);
    V_show = V_show - Va * (x_rtd - x_rtd(1)) / (x_rtd(end) - x_rtd(1));
    plot(x_rtd, V_show * 1000, 'Color', cmap(bi,:), 'LineWidth', 1.8, ...
        'DisplayName', sprintf('V_a = %.0f mV', Va*1000));
    hold on;
end
xlabel('Position (nm)');
ylabel('Potential (meV)');
title('RTD Potential Profile Under Applied Bias');
legend('Location', 'northeast');
grid on; set(gca, 'GridAlpha', 0.25);
exportgraphics(gcf, fullfile(out_dir, 'part2_rtd_potential_profiles.png'), ...
    'Resolution', 200);
close;

%% ========================================================================
%       PART 3 — SUPERLATTICE MINIBAND FORMATION
%  ========================================================================
%  A superlattice is a periodic sequence of N barriers and (N-1) wells.
%  As N increases, the isolated resonant transmission peaks of a double
%  barrier split and broaden into continuous minibands — the quantum
%  origin of band structure in periodic potentials.
%
%  We use GaAs/AlGaAs parameters for physical realism.
%  ========================================================================
fprintf('\n=== PART 3: Superlattice Miniband Formation ===\n');

sl_barrier_h = 0.23;   % eV
sl_barrier_w = 2.0;    % nm
sl_well_w    = 4.0;    % nm
N_barriers   = [2, 3, 5, 10];

E_sl = linspace(0.001, 0.25, 800);  % energy grid for sharp peaks

figure('Position', [60 60 900 750]);
T_sl_cache = cell(1, numel(N_barriers));   % cache for overlay plot

for ni = 1:numel(N_barriers)
    Nb = N_barriers(ni);
    fprintf('  N = %d barriers ... ', Nb);
    
    % Build the N-barrier superlattice potential and mass profile
    [x_sl, V_sl, m_sl] = build_superlattice(Nb, sl_barrier_h, ...
        sl_barrier_w, sl_well_w, m_well, m_barrier, 6.0);
    
    % Compute T(E)
    T_sl = zeros(size(E_sl));
    for ei = 1:numel(E_sl)
        [T_sl(ei), ~, ~] = solve_scattering_meff(x_sl, V_sl, ...
            E_sl(ei), m_sl, q, hbar, nm);
    end
    T_sl_cache{ni} = T_sl;   % cache for overlay plot
    
    subplot(numel(N_barriers), 1, ni);
    plot(E_sl * 1000, T_sl, 'Color', [0.0 0.35 0.65], 'LineWidth', 1.4);
    ylabel('T(E)');
    title(sprintf('N = %d barriers  (%d wells)', Nb, Nb-1));
    ylim([-0.05 1.1]);
    grid on; set(gca, 'GridAlpha', 0.2);
    if ni == numel(N_barriers)
        xlabel('Electron Energy (meV)');
    end
    fprintf('done\n');
end
sgtitle('Superlattice: Resonant Peaks → Miniband Formation', ...
    'FontSize', 14, 'FontWeight', 'bold');
exportgraphics(gcf, fullfile(out_dir, 'part3_superlattice_miniband.png'), ...
    'Resolution', 200);
close;

% ---- 3b. Overlay plot: all N on same axes (log scale) -------------------
%  Uses cached T(E) data from above — no recomputation needed.
figure('Position', [80 80 850 500]);
colors_n = [0.85 0.33 0.10; 0.47 0.67 0.19; 0.0 0.45 0.74; 0.49 0.18 0.56];
for ni = 1:numel(N_barriers)
    semilogy(E_sl*1000, T_sl_cache{ni} + 1e-15, 'Color', colors_n(ni,:), ...
        'LineWidth', 1.5, 'DisplayName', sprintf('N = %d', N_barriers(ni)));
    hold on;
end
xlabel('Electron Energy (meV)');
ylabel('Transmission T(E)');
title('Miniband Formation — Log Scale Comparison');
legend('Location', 'southeast');
grid on; set(gca, 'GridAlpha', 0.25);
ylim([1e-12 2]);
exportgraphics(gcf, fullfile(out_dir, 'part3_superlattice_overlay.png'), ...
    'Resolution', 200);
close;

% ---- 3c. Superlattice potential profile (N=5 example) -------------------
[x_sl5, V_sl5, ~] = build_superlattice(5, sl_barrier_h, ...
    sl_barrier_w, sl_well_w, m_well, m_barrier, 6.0);
figure('Position', [80 80 800 350]);
area(x_sl5, V_sl5*1000, 'FaceColor', [0.85 0.92 1.0], ...
    'EdgeColor', [0.0 0.35 0.65], 'LineWidth', 1.5);
xlabel('Position (nm)');
ylabel('Potential (meV)');
title('Superlattice Potential Profile (N = 5 barriers)');
grid on; set(gca, 'GridAlpha', 0.2);
exportgraphics(gcf, fullfile(out_dir, 'part3_superlattice_potential.png'), ...
    'Resolution', 200);
close;

%% ========================================================================
%       PART 4 — PARAMETRIC EXPLORATION (2-D HEATMAPS & SWEEPS)
%  ========================================================================
%  These plots map out how transmission depends on TWO parameters at once,
%  giving a much richer picture than the 1-D sweeps in the Python code.
%  ========================================================================
fprintf('\n=== PART 4: Parametric Exploration ===\n');

% ---- 4a. 2-D heatmap: T(Energy, Barrier Width) — single barrier --------
%  Shows how the exponential tunnelling suppression grows with width,
%  and how it vanishes once E crosses above V0.
fprintf('  4a: T(E, barrier width) heatmap...\n');

E_hm       = linspace(0.01, 0.6, 150);
W_hm       = linspace(0.3, 3.0, 60);
V0_hm      = 0.35;    % eV
T_heatmap  = zeros(numel(W_hm), numel(E_hm));

for wi = 1:numel(W_hm)
    x_hm = linspace(-max(W_hm)/2 - 2, max(W_hm)/2 + 2, 300);
    V_hm = rectangular_barrier(x_hm, V0_hm, W_hm(wi));
    for ei = 1:numel(E_hm)
        [T_heatmap(wi, ei), ~, ~] = solve_scattering(x_hm, V_hm, ...
            E_hm(ei), m0, q, hbar, nm);
    end
    if mod(wi, 15) == 0, fprintf('    %d%% ', round(wi/numel(W_hm)*100)); end
end
fprintf('\n');

figure('Position', [60 60 850 550]);
imagesc(E_hm, W_hm, log10(T_heatmap + 1e-16));
set(gca, 'YDir', 'normal');
colormap(parula(256));
cb = colorbar; cb.Label.String = 'log_{10} T(E)';
hold on;
xline(V0_hm, 'w--', 'LineWidth', 1.8);
text(V0_hm + 0.01, W_hm(end) - 0.15, sprintf('V_0 = %.2f eV', V0_hm), ...
    'Color', 'w', 'FontSize', 10, 'FontWeight', 'bold');
xlabel('Electron Energy (eV)');
ylabel('Barrier Width (nm)');
title('Transmission Heatmap: Energy vs Barrier Width');
exportgraphics(gcf, fullfile(out_dir, 'part4_heatmap_energy_width.png'), ...
    'Resolution', 200);
close;

% ---- 4b. Barrier height sweep — single barrier -------------------------
%  Shows how the tunnelling "cliff" shifts rightward with increasing V0,
%  and highlights the classical-transmission onset at E = V0.
fprintf('  4b: Barrier height sweep...\n');

heights_sweep = [0.15, 0.25, 0.35, 0.50, 0.75];  % eV
E_hs = linspace(0.01, 1.0, 250);
x_hs = linspace(-2.5, 2.5, 400);
bw_hs = 1.0;   % nm, fixed width

figure('Position', [60 60 850 550]);
cmap_h = turbo(numel(heights_sweep) + 2);
cmap_h = cmap_h(2:end-1, :);  % trim extremes for readability
for hi = 1:numel(heights_sweep)
    V_hs = rectangular_barrier(x_hs, heights_sweep(hi), bw_hs);
    T_hs = zeros(size(E_hs));
    for ei = 1:numel(E_hs)
        [T_hs(ei), ~, ~] = solve_scattering(x_hs, V_hs, E_hs(ei), ...
            m0, q, hbar, nm);
    end
    semilogy(E_hs, T_hs, 'Color', cmap_h(hi,:), 'LineWidth', 2.0, ...
        'DisplayName', sprintf('V_0 = %.2f eV', heights_sweep(hi)));
    hold on;
    % Mark E = V0 (classical threshold)
    xline(heights_sweep(hi), ':', 'Color', cmap_h(hi,:), ...
        'LineWidth', 1.0, 'HandleVisibility', 'off');
end
xlabel('Electron Energy (eV)');
ylabel('Transmission T(E)');
title(sprintf('Barrier Height Sweep (width = %.1f nm, free electron mass)', bw_hs));
legend('Location', 'southeast');
grid on; set(gca, 'GridAlpha', 0.25);
ylim([1e-12 2]);
annotation('textbox', [0.15 0.78 0.25 0.06], 'String', ...
    'Dotted lines: E = V_0 (classical threshold)', ...
    'EdgeColor', 'none', 'FontSize', 9, 'FontAngle', 'italic');
exportgraphics(gcf, fullfile(out_dir, 'part4_barrier_height_sweep.png'), ...
    'Resolution', 200);
close;

% ---- 4c. 2-D heatmap: T(Energy, Barrier Height) — single barrier -------
%  Continuous version of 4b.  The diagonal E = V0 line marks the boundary
%  between tunnelling (below) and classical transmission (above).
fprintf('  4c: T(E, barrier height) heatmap...\n');

H_hm2 = linspace(0.05, 0.8, 60);
E_hm2 = linspace(0.01, 1.0, 150);
bw_hm2 = 1.0;  % nm
T_hm2  = zeros(numel(H_hm2), numel(E_hm2));

x_hm2 = linspace(-2.5, 2.5, 300);
for hi = 1:numel(H_hm2)
    V_temp = rectangular_barrier(x_hm2, H_hm2(hi), bw_hm2);
    for ei = 1:numel(E_hm2)
        [T_hm2(hi, ei), ~, ~] = solve_scattering(x_hm2, V_temp, ...
            E_hm2(ei), m0, q, hbar, nm);
    end
    if mod(hi, 15) == 0, fprintf('    %d%% ', round(hi/numel(H_hm2)*100)); end
end
fprintf('\n');

figure('Position', [60 60 850 550]);
imagesc(E_hm2, H_hm2, log10(T_hm2 + 1e-16));
set(gca, 'YDir', 'normal');
colormap(parula(256));
cb = colorbar; cb.Label.String = 'log_{10} T(E)';
hold on;
% Draw E = V0 diagonal (classical boundary)
plot([H_hm2(1) H_hm2(end)], [H_hm2(1) H_hm2(end)], 'w--', 'LineWidth', 2.0);
text(0.55, 0.45, 'E = V_0', 'Color', 'w', 'FontSize', 11, ...
    'FontWeight', 'bold', 'Rotation', 30);
xlabel('Electron Energy (eV)');
ylabel('Barrier Height V_0 (eV)');
title(sprintf('Transmission Heatmap: Energy vs Barrier Height (width = %.1f nm)', bw_hm2));
exportgraphics(gcf, fullfile(out_dir, 'part4_heatmap_energy_height.png'), ...
    'Resolution', 200);
close;

% ---- 4d. Well width sweep — double barrier (resonant peak shifts) -------
%  As the quantum well widens, the quasi-bound energy levels drop and
%  more resonances fit below the barrier.  This is a direct visualisation
%  of the particle-in-a-box eigenvalue scaling E_n ~ 1/L^2.
fprintf('  4d: Well width sweep (double barrier)...\n');

well_widths = [2.0, 3.5, 5.0, 7.0, 10.0];   % nm
db_height   = 0.23;   % eV
db_bwidth   = 2.0;    % nm
E_ww        = linspace(0.001, 0.25, 400);

figure('Position', [60 60 850 550]);
cmap_w = [0.00 0.35 0.70; 0.00 0.55 0.35; 0.80 0.47 0.00;
          0.70 0.13 0.13; 0.50 0.18 0.56];
for wi = 1:numel(well_widths)
    ww = well_widths(wi);
    pad = max(8, ww/2 + db_bwidth + 4);
    x_ww = linspace(-pad, pad, 500);
    V_ww = double_barrier(x_ww, db_height, db_bwidth, ww);
    
    % Use GaAs effective mass for physical realism
    m_ww = m_well * ones(size(x_ww));
    hw_ww = ww / 2;
    m_ww((x_ww >= -hw_ww - db_bwidth) & (x_ww <= -hw_ww)) = m_barrier;
    m_ww((x_ww >=  hw_ww) & (x_ww <=  hw_ww + db_bwidth)) = m_barrier;
    
    T_ww = zeros(size(E_ww));
    for ei = 1:numel(E_ww)
        [T_ww(ei), ~, ~] = solve_scattering_meff(x_ww, V_ww, ...
            E_ww(ei), m_ww, q, hbar, nm);
    end
    semilogy(E_ww*1000, T_ww + 1e-15, 'Color', cmap_w(wi,:), ...
        'LineWidth', 1.8, 'DisplayName', sprintf('L_w = %.1f nm', ww));
    hold on;
end
xlabel('Electron Energy (meV)');
ylabel('Transmission T(E)');
title('Double Barrier: Resonant Peak Shift with Well Width');
legend('Location', 'southeast');
grid on; set(gca, 'GridAlpha', 0.25);
ylim([1e-10 2]);
annotation('textbox', [0.14 0.78 0.35 0.08], 'String', ...
    {sprintf('Barrier: %.0f meV, width = %.1f nm', db_height*1000, db_bwidth), ...
     'GaAs/AlGaAs effective masses'}, ...
    'EdgeColor', [0.5 0.5 0.5], 'BackgroundColor', [1 1 1 0.85], 'FontSize', 9);
exportgraphics(gcf, fullfile(out_dir, 'part4_well_width_sweep.png'), ...
    'Resolution', 200);
close;

% ---- 4e. 2-D heatmap: T(Energy, Well Width) — double barrier -----------
%  The continuous version.  Bright bands trace the resonant levels as
%  they sweep downward in energy with increasing well width.
fprintf('  4e: T(E, well width) heatmap...\n');

WW_hm     = linspace(1.5, 12.0, 50);
E_wwhm    = linspace(0.001, 0.25, 200);
T_wwhm    = zeros(numel(WW_hm), numel(E_wwhm));

for wi = 1:numel(WW_hm)
    ww = WW_hm(wi);
    pad = max(8, ww/2 + db_bwidth + 4);
    x_tmp = linspace(-pad, pad, 500);
    V_tmp = double_barrier(x_tmp, db_height, db_bwidth, ww);
    m_tmp = m_well * ones(size(x_tmp));
    hw_tmp = ww / 2;
    m_tmp((x_tmp >= -hw_tmp - db_bwidth) & (x_tmp <= -hw_tmp)) = m_barrier;
    m_tmp((x_tmp >=  hw_tmp) & (x_tmp <=  hw_tmp + db_bwidth)) = m_barrier;
    
    for ei = 1:numel(E_wwhm)
        [T_wwhm(wi, ei), ~, ~] = solve_scattering_meff(x_tmp, V_tmp, ...
            E_wwhm(ei), m_tmp, q, hbar, nm);
    end
    if mod(wi, 10) == 0, fprintf('    %d%% ', round(wi/numel(WW_hm)*100)); end
end
fprintf('\n');

figure('Position', [60 60 850 550]);
imagesc(E_wwhm*1000, WW_hm, log10(T_wwhm + 1e-16));
set(gca, 'YDir', 'normal');
colormap(hot(256));
cb = colorbar; cb.Label.String = 'log_{10} T(E)';
xlabel('Electron Energy (meV)');
ylabel('Well Width (nm)');
title('Resonant Level Map: Energy vs Well Width (Double Barrier)');
annotation('textbox', [0.15 0.80 0.32 0.08], 'String', ...
    {sprintf('Barrier: %.0f meV, width = %.1f nm', db_height*1000, db_bwidth), ...
     'Bright curves = resonant levels ~ 1/L_w^2'}, ...
    'EdgeColor', 'none', 'FontSize', 9, 'Color', 'w');
exportgraphics(gcf, fullfile(out_dir, 'part4_heatmap_energy_wellwidth.png'), ...
    'Resolution', 200);
close;

fprintf('  Part 4 complete.\n');

%% ======================== Summary =======================================
fprintf('\n========================================================\n');
fprintf('  All plots saved to ./%s/\n', out_dir);
fprintf('  Part 1: 2 figures  (validation overlay + error)\n');
fprintf('  Part 2: 4 figures  (T(E), wavefunction, I-V, profiles)\n');
fprintf('  Part 3: 3 figures  (miniband subplots, overlay, potential)\n');
fprintf('  Part 4: 5 figures  (heatmaps + parametric sweeps)\n');
fprintf('  Total:  14 figures\n');
fprintf('========================================================\n');

%% ========================================================================
%                         LOCAL FUNCTIONS
%  ========================================================================

function V = rectangular_barrier(x_nm, height_ev, width_nm)
% Single rectangular barrier centred at x = 0.
    V = zeros(size(x_nm));
    V(abs(x_nm) <= width_nm/2) = height_ev;
end

function V = double_barrier(x_nm, barrier_h, barrier_w, well_w)
% Two barriers with a quantum well between them, centred at x = 0.
    V = zeros(size(x_nm));
    hw = well_w / 2;
    V((x_nm >= -hw - barrier_w) & (x_nm <= -hw)) = barrier_h;
    V((x_nm >=  hw) & (x_nm <=  hw + barrier_w)) = barrier_h;
end

function [x, V, m_profile] = build_superlattice(N_barriers, barrier_h, ...
    barrier_w, well_w, m_well, m_barrier, contact_pad)
% BUILD_SUPERLATTICE  Generate the potential and effective-mass profiles
%   for an N-barrier periodic structure with flat contact regions.
%
%   Structure layout (example N=3):
%     contact | bar | well | bar | well | bar | contact
%
%   N barriers enclose (N-1) wells.  Total active region length is
%   N*barrier_w + (N-1)*well_w.  contact_pad nm of flat GaAs on each side.

    active_length = N_barriers * barrier_w + (N_barriers - 1) * well_w;
    total_length  = active_length + 2 * contact_pad;
    N_pts = max(800, round(total_length / 0.02));  % ~0.02 nm resolution
    
    x = linspace(-total_length/2, total_length/2, N_pts);
    V = zeros(size(x));
    m_profile = m_well * ones(size(x));
    
    % Active region starts at -active_length/2
    x_start = -active_length / 2;
    
    for ib = 1:N_barriers
        % Barrier ib starts at:
        b_left  = x_start + (ib-1) * (barrier_w + well_w);
        b_right = b_left + barrier_w;
        mask = (x >= b_left) & (x <= b_right);
        V(mask) = barrier_h;
        m_profile(mask) = m_barrier;
    end
end

function k = wave_number_scalar(E_ev, V_ev, mass, q_e, hbar_val)
% Complex wave number for energy E in a region of potential V.
    k = sqrt(2 * mass * complex((E_ev - V_ev) * q_e, 0)) / hbar_val;
end

function [T, R, psi] = solve_scattering(x_nm, V_ev, E_ev, mass, q_e, hbar_val, nm_val)
% Transfer-matrix scattering solver (uniform mass).
    N = numel(x_nm);
    dx = (x_nm(2) - x_nm(1)) * nm_val;
    
    k = arrayfun(@(v) wave_number_scalar(E_ev, v, mass, q_e, hbar_val), V_ev);
    
    M = eye(2);
    for i = 1:N-1
        kL = k(i); kR = k(i+1);
        if abs(kR) < 1e-18, kR = 1e-18; end
        if abs(kL) < 1e-18, kL = 1e-18; end
        I_mat = 0.5 * [1+kL/kR, 1-kL/kR; 1-kL/kR, 1+kL/kR];
        P_mat = [exp(1j*kR*dx), 0; 0, exp(-1j*kR*dx)];
        M = P_mat * I_mat * M;
    end
    
    r = -M(2,1)/M(2,2);
    t =  M(1,1) + M(1,2)*r;
    T = abs(t)^2 * real(k(end)) / real(k(1));
    R = abs(r)^2;
    
    state = [1; r];
    psi = zeros(1, N);
    psi(1) = state(1) + state(2);
    for i = 1:N-1
        kL = k(i); kR = k(i+1);
        if abs(kR) < 1e-18, kR = 1e-18; end
        if abs(kL) < 1e-18, kL = 1e-18; end
        I_mat = 0.5 * [1+kL/kR, 1-kL/kR; 1-kL/kR, 1+kL/kR];
        P_mat = [exp(1j*kR*dx), 0; 0, exp(-1j*kR*dx)];
        state = P_mat * I_mat * state;
        psi(i+1) = state(1) + state(2);
    end
    
    T = min(max(real(T), 0), 1);
    R = max(real(R), 0);
end

function [T, R, psi] = solve_scattering_meff(x_nm, V_ev, E_ev, m_eff, q_e, hbar_val, nm_val)
% SOLVE_SCATTERING_MEFF  Transfer-matrix solver with position-dependent
%   effective mass.  Interface matching uses (1/m*)(dpsi/dx) continuity
%   (Ben Daniel–Duke boundary conditions).
%
%   The transfer matrix at each interface accounts for the mass mismatch:
%     ratio = (m_right * k_left) / (m_left * k_right)
%   instead of simply k_left / k_right.

    N  = numel(x_nm);
    dx = (x_nm(2) - x_nm(1)) * nm_val;
    
    % Wave numbers with local effective mass
    k = zeros(1, N);
    for i = 1:N
        k(i) = sqrt(2 * m_eff(i) * complex((E_ev - V_ev(i)) * q_e, 0)) / hbar_val;
    end
    
    % ---- Total transfer matrix ------------------------------------------
    M = eye(2);
    for i = 1:N-1
        kL = k(i);  kR = k(i+1);
        mL = m_eff(i); mR = m_eff(i+1);
        if abs(kR) < 1e-18, kR = 1e-18; end
        if abs(kL) < 1e-18, kL = 1e-18; end
        
        % Ben Daniel-Duke: ratio = (mR*kL)/(mL*kR)
        ratio = (mR * kL) / (mL * kR);
        I_mat = 0.5 * [1 + ratio, 1 - ratio; ...
                        1 - ratio, 1 + ratio];
        P_mat = [exp(1j*kR*dx), 0; 0, exp(-1j*kR*dx)];
        M = P_mat * I_mat * M;
    end
    
    r = -M(2,1) / M(2,2);
    t =  M(1,1) + M(1,2) * r;
    
    % Flux-corrected transmission
    T = abs(t)^2 * real(k(end)/m_eff(end)) / real(k(1)/m_eff(1));
    R = abs(r)^2;
    
    % ---- Wavefunction ----------------------------------------------------
    state = [1; r];
    psi = zeros(1, N);
    psi(1) = state(1) + state(2);
    for i = 1:N-1
        kL = k(i); kR = k(i+1);
        mL = m_eff(i); mR = m_eff(i+1);
        if abs(kR) < 1e-18, kR = 1e-18; end
        if abs(kL) < 1e-18, kL = 1e-18; end
        ratio = (mR * kL) / (mL * kR);
        I_mat = 0.5 * [1 + ratio, 1 - ratio; ...
                        1 - ratio, 1 + ratio];
        P_mat = [exp(1j*kR*dx), 0; 0, exp(-1j*kR*dx)];
        state = P_mat * I_mat * state;
        psi(i+1) = state(1) + state(2);
    end
    
    T = min(max(real(T), 0), 1);
    R = max(real(R), 0);
end