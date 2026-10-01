% s01_ll_fva -- Loopless FVA across 15 context-specific extractions.
%
% For each cell:
%   1) load the fastcore-extracted model
%   2) apply sampling-style bounds (BOUND_EX * fc on 105 measured
%      exchanges, OCR cap on O2, 90% biomass floor)
%   3) run loopless fluxVariability (MILP)
%   4) classify each reaction as blocked / essential / flexible
%   5) write per-cell CSV
%
% Features:
%   - Resumes: skips cells whose llfva_{cell}.csv already exists.
%   - Live progress: prints [HH:MM:SS] event lines from workers via
%     DataQueue, so you can see progress without waiting for a whole
%     parfor iter to finish.
%   - Gurobi single-thread per worker (avoids core oversubscription
%     when running 12 workers on a 12-core machine).
%
% Prereq: run s00_preflight.m first to confirm allowLoops sign
% convention and packet biomass_max reproduction.

if ~exist('N_WORKERS','var'), N_WORKERS = 12; end
if ~exist('ALLOWLOOPS_ARG','var'), ALLOWLOOPS_ARG = 0; end   % 0 = loopless
if ~exist('CHUNK_SIZE','var'), CHUNK_SIZE = 250; end         % rxns per FVA sub-call (progress granularity)

here = fileparts(mfilename('fullpath'));
repo = fileparts(here);
addpath(fullfile(repo,'functions'));
res  = fullfile(here,'results');
if ~exist(res,'dir'), mkdir(res); end
diary(fullfile(res,'s01_ll_fva.txt')); diary on;
cleanupObj = onCleanup(@() diary('off'));

try, initCobraToolbox(false); catch, end
changeCobraSolver('gurobi','LP');
changeCobraSolver('gurobi','MILP');

CELLS = {'BL_24','BL_48','BL_72','BL_96','BL_120', ...
         'D_24','D_48','D_72','D_96','D_120', ...
         'LD_24','LD_48','LD_72','LD_96','LD_120'};
OCR = containers.Map({24,48,72,96,120},{97,173,227,273,313});
BIOMASS_RXN = 'MAR00021'; O2_EX = 'MAR09048'; BOUND_EX = 10.0;

mb = readtable(fullfile(repo,'results','qc','phase2_exchange_bounds_long.csv'), ...
    'FileType','text','Delimiter',',','VariableNamingRule','preserve');

% --- resume: filter out cells that already have a CSV ---
todo_mask = false(numel(CELLS),1);
for kc = 1:numel(CELLS)
    csv_path = fullfile(res, sprintf('llfva_%s.csv', CELLS{kc}));
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
% Also expose loopless_fva/ on workers so fluxVariability_tl.m is found
% and takes precedence over COBRA's stock fluxVariability.
llfva_dir_esc  = strrep(here, '''', '''''');
pctRunOnAll(sprintf('addpath(''%s'');', llfva_dir_esc));
pctRunOnAll('ensure_solver();');

% --- Gurobi MILP time limit per solve ---
% Without this, some loopless MILPs will grind for hours trying to
% prove optimality on hard branch-and-bound trees. TimeLimit=30s means
% each MILP returns after 30s with its best-known bound; loopless
% constraints are still enforced -- reported range on the few
% pathological reactions is just slightly wider than proven optimum.
% This is the standard fix in loopless FVA implementations.
if ~exist('MILP_TIME_LIMIT','var'), MILP_TIME_LIMIT = 30; end
try
    pctRunOnAll(sprintf('changeCobraSolverParams(''MILP'', ''timeLimit'', %g);', MILP_TIME_LIMIT));
    fprintf('Gurobi MILP timeLimit = %g s on all workers.\n', MILP_TIME_LIMIT);
catch ME
    warning('Could not set MILP timeLimit: %s', ME.message);
end

% --- progress DataQueue (workers -> main session) ---
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
    ensure_solver();   % insurance if a worker gets recycled mid-sweep
    cname = CELLS_TODO{kc};
    parts = strsplit(cname,'_'); cond = parts{1}; hpf = str2double(parts{2});
    ocrv  = OCR(hpf);

    send(progQ, {cname, 'loading model...'});
    S_data = load(fullfile(repo,'reviewer_packet_v4_slim','v7_models', ...
        sprintf('trans_rfastcormics_%s.mat', cname)));
    m = S_data.parent;

    % --- v7 canonical bounds via shared apply_context_bounds (2026-09-30) ---
    [m, ~, ~] = apply_context_bounds(m, cond, hpf, 'A');
    Snnz  = full(sum(m.S ~= 0, 1));
    is_ex = (Snnz == 1)';
    m = changeObjective(m, BIOMASS_RXN);
    sol = optimizeCbModel(m, 'max');
    bio_max = sol.f;
    j_bio = find(strcmp(m.rxns, BIOMASS_RXN), 1);

    n_rxns_c = numel(m.rxns);
    n_chunks = ceil(n_rxns_c / CHUNK_SIZE);
    send(progQ, {cname, sprintf('bounds applied, biomass_max=%.4f, FVA over %d rxns in %d chunks of %d', ...
        bio_max, n_rxns_c, n_chunks, CHUNK_SIZE)});

    % --- loopless FVA (chunked for intra-cell progress reporting) ---
    % Chunking is math-neutral: each reaction's min/max is an independent
    % MILP whose constraints don't depend on which other reactions are
    % probed. Overhead is repeated initial-FBA + loopless-MILP setup
    % per chunk (few seconds each).
    minFlux = nan(n_rxns_c, 1);
    maxFlux = nan(n_rxns_c, 1);
    t = tic;
    for kchunk = 1:n_chunks
        lo = (kchunk-1)*CHUNK_SIZE + 1;
        hi = min(kchunk*CHUNK_SIZE, n_rxns_c);
        chunk_rxns = m.rxns(lo:hi);
        t_c = tic;
        [mn, mx] = fluxVariability_tl(m, 90, 'max', chunk_rxns, 0, ALLOWLOOPS_ARG, 'FBA');
        minFlux(lo:hi) = mn;
        maxFlux(lo:hi) = mx;
        pct = 100*hi/n_rxns_c;
        send(progQ, {cname, sprintf('%3.0f%%  (%d/%d rxns, chunk %d/%d in %.1f min, cell %.1f min)', ...
            pct, hi, n_rxns_c, kchunk, n_chunks, toc(t_c)/60, toc(t)/60)});
    end
    dt_min = toc(t)/60;

    tol = 1e-6;
    is_blocked   = (abs(minFlux) < tol) & (abs(maxFlux) < tol);
    is_essential = (minFlux > tol)  | (maxFlux < -tol);
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
                  is_blocked, is_essential, is_flexible, ...
                  'VariableNames', {'rxn','name','subsystem', ...
                     'min_flux','max_flux','range','blocked','essential','flexible'});
    writetable(T_out, fullfile(res, sprintf('llfva_%s.csv', cname)));

    send(progQ, {cname, sprintf('DONE in %.1f min  (blocked=%d, ess=%d, flex=%d)', ...
        dt_min, sum(is_blocked), sum(is_essential), sum(is_flexible))});

    out{kc} = struct('cell', cname, 'n_rxns', numel(m.rxns), ...
        'n_blocked', sum(is_blocked), 'n_essential', sum(is_essential), ...
        'n_flexible', sum(is_flexible), 'biomass_max', bio_max, ...
        'dt_min', dt_min);
end

fprintf('\nAll parfor iters returned in %.1f min wall.\n', toc(t_run)/60);

% --- rebuild full summary by reading every per-cell CSV back in ---
rows = cell(numel(CELLS), 6);
for kc = 1:numel(CELLS)
    cname = CELLS{kc};
    csv_path = fullfile(res, sprintf('llfva_%s.csv', cname));
    if exist(csv_path, 'file')
        Tc = readtable(csv_path);
        rows(kc,:) = {cname, height(Tc), sum(Tc.blocked), ...
                      sum(Tc.essential), sum(Tc.flexible), NaN};
    else
        rows(kc,:) = {cname, NaN, NaN, NaN, NaN, NaN};
    end
end
% fill biomass_max from this run's out{}
for kc = 1:n_cells
    r = out{kc};
    idx = find(strcmp(CELLS, r.cell), 1);
    rows{idx, 6} = r.biomass_max;
end
T = cell2table(rows, 'VariableNames', ...
    {'cell','n_rxns','n_blocked','n_essential','n_flexible','biomass_max'});
writetable(T, fullfile(res,'llfva_summary.csv'));

% --- biomass_max report (this session, sampling-bounds regime) ---
% Note: s125_biomass_max_with_flavin_forced.csv is at EXTRACTION bounds,
% which are looser than the sampling-bounds regime we're using here (fc
% scaling on 105 exchanges + OCR cap). So our biomass_max values are
% expected to be <= s125 reference. Comparison intentionally omitted.
if n_cells > 0
    fprintf('\n=== biomass_max at sampling bounds (this session) ===\n');
    for kc = 1:n_cells
        r = out{kc};
        fprintf('  %-7s biomass_max=%.4f  FVA in %.1f min\n', ...
            r.cell, r.biomass_max, r.dt_min);
    end
end
fprintf('DONE\n');
