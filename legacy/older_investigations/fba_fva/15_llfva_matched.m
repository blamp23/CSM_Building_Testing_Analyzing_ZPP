% s02b_ll_fva_matched_floor -- ll-FVA with a MATCHED biomass floor across
% conditions at each hpf, so between-condition Δrange comparisons aren't
% confounded by the Frame A polytope-size artefact.
%
% For each hpf, we compute biomass_max_full for BL, D, and LD, then
% set every cell's biomass floor to
%
%     0.9 * min(biomass_max_full over the three conditions at that hpf).
%
% Parallelism model: ONE cell at a time in the main session, with a
% parfor loop over the ~25 reaction chunks inside each cell. This makes
% each cell finish in ~n_chunks / n_workers x avg_chunk_time
% (roughly 10-20 min per cell instead of 3-6 hr), gives visible
% progress (one cell fully completes and writes its CSV before the
% next starts), and avoids the Gurobi thread over-subscription that
% cell-level parfor caused in the original s02 run.
%
% All other bounds are identical to s02_ll_fva.m (BOUND_EX * fc metabolomics
% envelope, OCR cap on O2, loopless MILP with 30 s Gurobi time limit,
% 250-rxn chunking with progress reporting, resume on per-cell CSVs).
%
% Output directory:  loopless_fva/results_matched/
% (kept separate from results/ so both regimes coexist).

if ~exist('N_WORKERS','var'), N_WORKERS = 10; end
if ~exist('ALLOWLOOPS_ARG','var'), ALLOWLOOPS_ARG = 0; end
if ~exist('CHUNK_SIZE','var'), CHUNK_SIZE = 250; end
if ~exist('MILP_TIME_LIMIT','var'), MILP_TIME_LIMIT = 30; end

here = fileparts(mfilename('fullpath'));
repo = fileparts(here);
addpath(fullfile(repo,'functions'));
res  = fullfile(here,'results_matched');
if ~exist(res,'dir'), mkdir(res); end
diary(fullfile(res,'s02b_ll_fva_matched.txt')); diary on;
cleanupObj = onCleanup(@() diary('off'));

try, initCobraToolbox(false); catch, end
changeCobraSolver('gurobi','LP');
changeCobraSolver('gurobi','MILP');

CELLS = {'BL_24','BL_48','BL_72','BL_96','BL_120', ...
         'D_24','D_48','D_72','D_96','D_120', ...
         'LD_24','LD_48','LD_72','LD_96','LD_120'};
HPFS  = [24, 48, 72, 96, 120];
OCR   = containers.Map({24,48,72,96,120},{97,173,227,273,313});
BIOMASS_RXN = 'MAR00021'; O2_EX = 'MAR09048'; BOUND_EX = 10.0;
MODEL_DIR = fullfile(repo,'reviewer_packet_v4_slim','v7_models');

mb = readtable(fullfile(repo,'results','qc','phase2_exchange_bounds_long.csv'), ...
    'FileType','text','Delimiter',',','VariableNamingRule','preserve');

% ------------------------------------------------------------------
% Step 1: compute biomass_max_full per cell (serial, LP only, ~1 min total)
% ------------------------------------------------------------------
fprintf('\n=== Step 1: compute biomass_max_full per cell ===\n');
bio_max_map = containers.Map();
for kc = 1:numel(CELLS)
    cname = CELLS{kc};
    parts = strsplit(cname,'_'); cond = parts{1}; hpf = str2double(parts{2});
    ocrv  = OCR(hpf);
    S_data = load(fullfile(MODEL_DIR, sprintf('trans_rfastcormics_%s.mat', cname)));
    m = S_data.parent;
    % --- v7 canonical bounds via shared apply_context_bounds (2026-09-30) ---
    [m, ~, ~] = apply_context_bounds(m, cond, hpf, 'none');
    m = changeObjective(m, BIOMASS_RXN);
    sol = optimizeCbModel(m, 'max');
    bio_max_map(cname) = sol.f;
    fprintf('  [%s] biomass_max_full = %.4f\n', cname, sol.f);
end

% Per-hpf matched floor = 0.9 * min(biomass_max_full over BL/D/LD)
matched_floor = containers.Map('KeyType','double','ValueType','double');
fprintf('\nMatched biomass floors per hpf (0.9 * min over conditions):\n');
for h = HPFS
    vals = [bio_max_map(sprintf('BL_%d',h)), ...
            bio_max_map(sprintf('D_%d',h)), ...
            bio_max_map(sprintf('LD_%d',h))];
    matched_floor(h) = 0.9 * min(vals);
    fprintf('  hpf=%3d:  BL=%.4f  D=%.4f  LD=%.4f  ->  floor=%.4f\n', ...
        h, vals(1), vals(2), vals(3), matched_floor(h));
end

% ------------------------------------------------------------------
% Step 2: resume-aware parfor over cells at the matched floor
% ------------------------------------------------------------------
todo_mask = false(numel(CELLS),1);
for kc = 1:numel(CELLS)
    csv_path = fullfile(res, sprintf('llfva_%s.csv', CELLS{kc}));
    if ~exist(csv_path,'file'), todo_mask(kc) = true; end
end
CELLS_TODO = CELLS(todo_mask);
n_cells = numel(CELLS_TODO);
fprintf('\nResume: %d of %d cells already done, %d to run.\n', ...
    sum(~todo_mask), numel(CELLS), n_cells);

pool = gcp('nocreate');
if isempty(pool)
    pool = parpool('local', min(N_WORKERS, max(n_cells,1)));
end

helper_dir     = fullfile(repo,'functions');
helper_dir_esc = strrep(helper_dir, '''', '''''');
pctRunOnAll(sprintf('addpath(''%s'');', helper_dir_esc));
llfva_dir_esc  = strrep(here, '''', '''''');
pctRunOnAll(sprintf('addpath(''%s'');', llfva_dir_esc));
pctRunOnAll('ensure_solver();');
try
    pctRunOnAll(sprintf('changeCobraSolverParams(''MILP'', ''timeLimit'', %g);', MILP_TIME_LIMIT));
    fprintf('Gurobi MILP timeLimit = %g s on all workers.\n', MILP_TIME_LIMIT);
catch ME
    warning('Could not set MILP timeLimit: %s', ME.message);
end

if n_cells > 0
    progQ = parallel.pool.DataQueue;
    afterEach(progQ, @(msg) fprintf('[%s] %-7s %s\n', ...
        datestr(now,'HH:MM:SS'), msg{1}, msg{2}));
else
    progQ = [];
end

out = cell(n_cells,1);
t_run = tic;

% -------- outer serial loop over cells, inner parfor over chunks --------
for kc = 1:n_cells
    cname = CELLS_TODO{kc};
    parts = strsplit(cname,'_'); cond = parts{1}; hpf = str2double(parts{2});
    ocrv  = OCR(hpf);
    floor_val = matched_floor(hpf);

    fprintf('\n[%d/%d]  [%s]  %s   loading model...\n', kc, n_cells, ...
        datestr(now,'HH:MM:SS'), cname);
    S_data = load(fullfile(MODEL_DIR, sprintf('trans_rfastcormics_%s.mat', cname)));
    m = S_data.parent;
    % --- v7 canonical bounds + matched floor (2026-09-30) ---
    % apply_context_bounds preserves stored v7 boundary; matched floor
    % overrides its Frame-A biomass_lb after the fact.
    [m, ~, ~] = apply_context_bounds(m, cond, hpf, 'matched', floor_val / 0.9);
    Snnz  = full(sum(m.S ~= 0, 1));
    is_ex = (Snnz == 1)';
    m = changeObjective(m, BIOMASS_RXN);
    sol_max = optimizeCbModel(m, 'max');
    bio_max = sol_max.f;
    j_bio = find(strcmp(m.rxns, BIOMASS_RXN), 1);
    m.lb(j_bio) = floor_val;  % MATCHED FLOOR (explicit, in case rounding)

    n_rxns_c = numel(m.rxns);
    n_chunks = ceil(n_rxns_c / CHUNK_SIZE);
    fprintf('[%s]  [%s] own bio_max=%.4f, matched floor=%.4f, FVA over %d rxns in %d chunks\n', ...
        datestr(now,'HH:MM:SS'), cname, bio_max, floor_val, n_rxns_c, n_chunks);

    % chunk ranges precomputed so parfor can slice cleanly
    los = zeros(n_chunks,1); his = zeros(n_chunks,1);
    for kchunk = 1:n_chunks
        los(kchunk) = (kchunk-1)*CHUNK_SIZE + 1;
        his(kchunk) = min(kchunk*CHUNK_SIZE, n_rxns_c);
    end
    chunk_rxns_list = arrayfun(@(i) m.rxns(los(i):his(i)), 1:n_chunks, ...
                               'UniformOutput', false);

    minChunks = cell(n_chunks,1); maxChunks = cell(n_chunks,1);
    tChunks   = zeros(n_chunks,1);
    t_cell = tic;

    % --- parfor over chunks; m and ALLOWLOOPS_ARG broadcast once ---
    parfor kchunk = 1:n_chunks
        ensure_solver();
        t_c = tic;
        [mn, mx] = fluxVariability_tl(m, 90, 'max', chunk_rxns_list{kchunk}, ...
                                       0, ALLOWLOOPS_ARG, 'FBA');
        minChunks{kchunk} = mn;
        maxChunks{kchunk} = mx;
        tChunks(kchunk)   = toc(t_c);
        send(progQ, {cname, sprintf('chunk %2d/%d done in %.1f min', ...
            kchunk, n_chunks, toc(t_c)/60)});
    end

    % assemble
    minFlux = nan(n_rxns_c,1); maxFlux = nan(n_rxns_c,1);
    for kchunk = 1:n_chunks
        minFlux(los(kchunk):his(kchunk)) = minChunks{kchunk};
        maxFlux(los(kchunk):his(kchunk)) = maxChunks{kchunk};
    end
    dt_min = toc(t_cell)/60;

    tol = 1e-6;
    is_blocked   = (abs(minFlux) < tol) & (abs(maxFlux) < tol);
    is_required  = (minFlux > tol)  | (maxFlux < -tol);
    is_flexible  = (maxFlux - minFlux) > 1.0;

    subs = strings(numel(m.subSystems), 1);
    for i = 1:numel(m.subSystems)
        v = m.subSystems{i};
        while iscell(v)
            if isempty(v), v = ''; break; end
            v = v{1};
        end
        if isempty(v), subs(i) = ""; else, subs(i) = string(v); end
    end

    T_out = table(m.rxns, m.rxnNames, subs, ...
                  minFlux, maxFlux, maxFlux - minFlux, ...
                  is_blocked, is_required, is_flexible, ...
                  'VariableNames', {'rxn','name','subsystem', ...
                     'min_flux','max_flux','range','blocked','required_at_matched_floor','range_gt_1'});
    writetable(T_out, fullfile(res, sprintf('llfva_%s.csv', cname)));

    fprintf('[%s]  [%s]  DONE in %.1f min (blocked=%d, required=%d, flex=%d)\n', ...
        datestr(now,'HH:MM:SS'), cname, dt_min, ...
        sum(is_blocked), sum(is_required), sum(is_flexible));

    out{kc} = struct('cell', cname, 'n_rxns', numel(m.rxns), ...
        'n_blocked', sum(is_blocked), 'n_required', sum(is_required), ...
        'n_flexible', sum(is_flexible), 'biomass_max', bio_max, ...
        'matched_floor', floor_val, 'dt_min', dt_min);
end

fprintf('\nAll cells done in %.1f min wall.\n', toc(t_run)/60);

% summary
rows = cell(numel(CELLS),7);
for kc = 1:numel(CELLS)
    cname = CELLS{kc};
    csv_path = fullfile(res, sprintf('llfva_%s.csv', cname));
    if exist(csv_path,'file')
        Tc = readtable(csv_path);
        rows(kc,:) = {cname, height(Tc), sum(Tc.blocked), ...
                      sum(Tc.required_at_matched_floor), sum(Tc.range_gt_1), NaN, NaN};
    else
        rows(kc,:) = {cname, NaN, NaN, NaN, NaN, NaN, NaN};
    end
end
for kc = 1:n_cells
    r = out{kc};
    idx = find(strcmp(CELLS, r.cell), 1);
    rows{idx, 6} = r.biomass_max;
    rows{idx, 7} = r.matched_floor;
end
T = cell2table(rows, 'VariableNames', ...
    {'cell','n_rxns','n_blocked','n_required','n_flexible','biomass_max','matched_floor'});
writetable(T, fullfile(res,'llfva_summary_matched.csv'));
fprintf('DONE\n');
