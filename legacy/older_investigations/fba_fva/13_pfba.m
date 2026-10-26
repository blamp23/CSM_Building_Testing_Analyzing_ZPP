% s01_pfba -- Parallel pFBA on the 15 v6 flavin-forced extractions.
%
% For each cell:
%   1) load v6 flavin-forced model
%   2) apply sampling-style bounds (same as s02_ll_fva)
%   3) pFBA: max biomass then min L1 flux at that biomass (biomass floor
%      0.9 * biomass_max_full)
%   4) write per-cell canonical flux vector CSV
%
% Features (mirror s02_ll_fva):
%   - Resumes: skips cells where pfba_{cell}.csv already exists.
%   - Live progress via DataQueue.
%   - Parallel over cells via parfor with ensure_solver on each worker.
%
% LP only (no MILP), so ~30-60 s per cell. On 12 workers, ~1-2 min wall.

if ~exist('N_WORKERS','var'), N_WORKERS = 12; end

here = fileparts(mfilename('fullpath'));
repo = fileparts(here);
addpath(fullfile(repo,'functions'));
res  = fullfile(here,'results');
if ~exist(res,'dir'), mkdir(res); end
diary(fullfile(res,'s01_pfba.txt')); diary on;
cleanupObj = onCleanup(@() diary('off'));

try, initCobraToolbox(false); catch, end
changeCobraSolver('gurobi','LP');

CELLS = {'BL_24','BL_48','BL_72','BL_96','BL_120', ...
         'D_24','D_48','D_72','D_96','D_120', ...
         'LD_24','LD_48','LD_72','LD_96','LD_120'};
OCR = containers.Map({24,48,72,96,120},{97,173,227,273,313});
BIOMASS_RXN = 'MAR00021'; O2_EX = 'MAR09048'; BOUND_EX = 10.0;
BIO_LB_FRAC = 0.9;

mb = readtable(fullfile(repo,'results','qc','phase2_exchange_bounds_long.csv'), ...
    'FileType','text','Delimiter',',','VariableNamingRule','preserve');

% --- resume: filter out cells that already have a CSV ---
todo_mask = false(numel(CELLS),1);
for kc = 1:numel(CELLS)
    csv_path = fullfile(res, sprintf('pfba_%s.csv', CELLS{kc}));
    if ~exist(csv_path, 'file')
        todo_mask(kc) = true;
    end
end
CELLS_TODO = CELLS(todo_mask);
n_cells = numel(CELLS_TODO);
fprintf('Resume: %d of %d cells already done, %d to run.\n', ...
    sum(~todo_mask), numel(CELLS), n_cells);
if n_cells == 0
    fprintf('All cells complete. Skipping to summary.\n');
else
    fprintf('To run: %s\n', strjoin(CELLS_TODO, ', '));
end

% --- pool setup ---
pool = gcp('nocreate');
if isempty(pool)
    pool = parpool('local', min(N_WORKERS, max(n_cells,1)));
end

helper_dir     = fullfile(repo,'functions');
helper_dir_esc = strrep(helper_dir, '''', '''''');
pctRunOnAll(sprintf('addpath(''%s'');', helper_dir_esc));
pctRunOnAll('ensure_solver();');

% --- Gurobi threads=1 per worker (avoid oversubscription) ---
try
    pctRunOnAll('changeCobraSolverParams(''LP'', ''threads'', 1);');
    fprintf('Gurobi threads=1 on all workers.\n');
catch ME
    warning('Could not set Gurobi thread param: %s', ME.message);
end

% --- progress DataQueue ---
if n_cells > 0
    progQ = parallel.pool.DataQueue;
    afterEach(progQ, @(msg) fprintf('[%s] %-7s %s\n', ...
        datestr(now,'HH:MM:SS'), msg{1}, msg{2}));
else
    progQ = [];
end

out = cell(n_cells, 1);
t_run = tic;

parfor kc = 1:n_cells
    ensure_solver();
    cname = CELLS_TODO{kc};
    parts = strsplit(cname,'_'); cond = parts{1}; hpf = str2double(parts{2});
    ocrv  = OCR(hpf);

    send(progQ, {cname, 'loading model...'});
    S_data = load(fullfile(repo,'reviewer_packet_v4_slim','v7_models', ...
        sprintf('trans_rfastcormics_%s.mat', cname)));
    m = S_data.parent;

    % --- v7 canonical bounds via shared apply_context_bounds (2026-09-30) ---
    % Replaces the old "reset all exchanges to ±BOUND_EX" block which
    % discarded the v7 boundary. apply_context_bounds preserves stored
    % bounds and only edits (1) fc on the 105 measured, (2) OCR cap.
    [m, ~, ~] = apply_context_bounds(m, cond, hpf, 'A');
    % Note: 'A' inside apply_context_bounds already sets biomass_lb =
    % 0.9 * biomass_max computed under (1)+(2).
    Snnz  = full(sum(m.S ~= 0, 1));
    is_ex = (Snnz == 1)';
    j_bio = find(strcmp(m.rxns, BIOMASS_RXN), 1);
    m = changeObjective(m, BIOMASS_RXN);
    sol_max = optimizeCbModel(m, 'max');
    bio_max = sol_max.f;

    % --- pFBA: two-stage LP, L1 min of fluxes at biomass_max ---
    t = tic;
    sol_pfba = optimizeCbModel(m, 'max', 'one');
    dt_s = toc(t);

    if sol_pfba.stat ~= 1
        send(progQ, {cname, sprintf('pFBA LP not optimal (stat=%d) -- SKIPPED', sol_pfba.stat)});
        out{kc} = struct('cell', cname, 'ok', false);
        continue;
    end

    v = sol_pfba.x;
    n_active = sum(abs(v) > 1e-8);
    L1_internal = sum(abs(v(~is_ex)));

    T_out = table(m.rxns, v, m.lb, m.ub, ...
        'VariableNames', {'rxn','v_pfba','lb','ub'});
    writetable(T_out, fullfile(res, sprintf('pfba_%s.csv', cname)));

    send(progQ, {cname, sprintf('DONE in %.1f s  (bio_max=%.4f, bio_pfba=%.4f, active=%d, L1_int=%.4g)', ...
        dt_s, bio_max, sol_pfba.f, n_active, L1_internal)});

    out{kc} = struct('cell', cname, 'cond', cond, 'hpf', hpf, ...
        'biomass_max', bio_max, 'biomass_pfba', sol_pfba.f, ...
        'n_active', n_active, 'L1_internal', L1_internal, 'ok', true);
end

fprintf('\nAll parfor iters returned in %.1f min wall.\n', toc(t_run)/60);

% --- summary across all cells (this session + prior) ---
summary_rows = cell(0, 6);
for kc = 1:numel(CELLS)
    cname = CELLS{kc};
    csv_path = fullfile(res, sprintf('pfba_%s.csv', cname));
    if exist(csv_path, 'file')
        Tc = readtable(csv_path);
        parts = strsplit(cname,'_'); cond = parts{1}; hpf = str2double(parts{2});
        n_active = sum(abs(Tc.v_pfba) > 1e-8);
        % re-flag is_ex from bounds columns (single-connectivity assumed)
        L1_all = sum(abs(Tc.v_pfba));
        summary_rows(end+1, :) = { cname, cond, hpf, n_active, L1_all, height(Tc) };
    end
end
if ~isempty(summary_rows)
    Tsum = cell2table(summary_rows, 'VariableNames', ...
        {'cell','condition','hpf','n_active','L1_all','n_rxns'});
    writetable(Tsum, fullfile(res,'pfba_summary.csv'));
    fprintf('wrote pfba_summary.csv\n');
end

% --- biomass_max sanity check vs v6 reference ---
ref_csv = fullfile(repo,'reviewer_packet_v4_slim','v7_models', ...
    's125_biomass_max_with_flavin_forced.csv');
if exist(ref_csv,'file') && n_cells > 0
    Tref = readtable(ref_csv,'FileType','text','Delimiter',',','VariableNamingRule','preserve');
    ref_map = containers.Map(Tref.cell, num2cell(Tref.biomass_max));
    fprintf('\n=== biomass_max reconciliation (this session vs v6 reference) ===\n');
    for kc = 1:n_cells
        r = out{kc};
        if ~r.ok, continue; end
        p = ref_map(r.cell);
        fprintf('  %-7s got %.4f  ref %.4f  d=%+.4f\n', ...
            r.cell, r.biomass_max, p, r.biomass_max - p);
    end
end
fprintf('DONE\n');
