% s140 -- rFASTCORMICS extraction on the v7 curated bases.
%
% Same core-set logic as s02_extract_rfastcormics.m and s125, but:
%   - Loads baked_ocr_{24,48,72,96,120}_v7.mat (v7 curated: Part 1 + Part 2
%     + ledger closes + pool coef edits + inositol/vit E/carotene/lipoprotein
%     medium_uptake + ysl_lipid_delivery)
%   - NO forced MAR IDs beyond the standard L1/L2/L3 (removed the s125
%     flavin-forcing block MAR06506/MAR06508/MAR09143 — riboflavin exchange
%     stays open at ±10 default via ledger keep_exchange, so FASTCORE picks
%     it up automatically)
%   - Core = expressed ∪ biomass ∪ open exchanges (ledger medium_uptake,
%     ysl_lipid_delivery, 105 measured, plus the standard L1 bulk medium)
%   - Writes to reviewer_packet_v4_slim/v7_models/
%
% Run:
%   matlab -batch "cd('reviewer_packet_v4_slim'); s140_extract_v7"

clear cell cname

if ~exist('EPSILON','var'),   EPSILON = 1e-4; end
if ~exist('N_WORKERS','var'), N_WORKERS = 12; end   % local Mac, per user 2026-09-29

here = fileparts(mfilename('fullpath'));
repo = fileparts(here);
addpath(fullfile(repo,'functions'));

pdir_out = fullfile(here,'v7_models');
if ~exist(pdir_out,'dir'), mkdir(pdir_out); end

diary(fullfile(pdir_out,'s140_extract.txt')); diary on
cleanupObj = onCleanup(@() diary('off'));

fprintf('=== s140 v7 extraction (no forced MAR IDs) ===\n');
fprintf('date: %s   epsilon=%g   N_WORKERS=%d\n\n', datestr(now,'yyyy-mm-dd HH:MM:SS'), EPSILON, N_WORKERS);

BIOMASS_RXN = 'MAR00021';

% Standard L1 (bulk medium + biomass reaction anchors). No flavin forcing.
L1 = {'MAR09077','MAR09081','MAR09082','MAR13072','MAR09150','MAR09074', ...
      'MAR09047','MAR09079','MAR09076','MAR09080','MAR13073','MAR09048', ...
      'MAR09058','MAR09078','MAR09072','MAR09073','MAR11420','MAR09096', ...
      'MAR09378','MAR00021','MAR10023','MAR10024','MAR06916','MAR03964'};
fprintf('L1 (standard, no flavin forcing): %d rxns\n', numel(L1));

T2 = readtable(fullfile(repo,'results','qc','phase2_essential_exchanges.csv'), ...
               'FileType','text','Delimiter',',','VariableNamingRule','preserve');
L2 = T2.rxn_id; if ~iscell(L2), L2 = cellstr(L2); end
T3 = readtable(fullfile(repo,'results','qc','phase2_exchange_bounds_long.csv'), ...
               'FileType','text','Delimiter',',','VariableNamingRule','preserve');
L3 = unique(T3.ex_rxn); if ~iscell(L3), L3 = cellstr(L3); end
fprintf('L1=%d  L2=%d  L3=%d\n', numel(L1), numel(L2), numel(L3));

% ---- v7 ledger keep-open exchanges (medium_uptake + ysl_lipid_delivery + keep_exchange)
LEDGER = fullfile(here, 's130_v7_ledger.csv');
Tl = readtable(LEDGER, 'FileType','text','Delimiter',',','VariableNamingRule','preserve');
keep_mask = strcmp(Tl.action,'keep_exchange') | strcmp(Tl.action,'medium_uptake') | ...
            strcmp(Tl.action,'ysl_lipid_delivery');
L4 = Tl.exchange_MAR_id(keep_mask);
L4 = L4(~cellfun('isempty', L4));
fprintf('L4 (v7 ledger keep-open, includes ysl_lipid_delivery): %d\n\n', numel(L4));

% ---- load per-hpf v7 parents (all 5) -----------------------------------
hpf_list = [24 48 72 96 120];
parents_by_hpf = struct(); anchors_by_hpf = struct();
for h = hpf_list
    p = fullfile(repo,'ocr_anchored_extraction','models', sprintf('baked_ocr_%d_v7.mat', h));
    BL = load(p); parent = ensureCobra(BL.m);
    sol_p = optimizeCbModel(changeObjective(parent, BIOMASS_RXN), 'max');
    parents_by_hpf.(sprintf('h%d', h)) = parent;
    anchors_by_hpf.(sprintf('h%d', h)) = 0.9 * sol_p.f;
    fprintf('v7 parent hpf=%d: bio_max=%.4g  anchor=%.4g\n', h, sol_p.f, 0.9*sol_p.f);
end
fprintf('\n');

% ---- gene labels + gene map (same as s02) ------------------------------
labels_csv = fullfile(repo,'rfastcoreomics','results','gene_labels_by_cell.csv');
LT = readtable(labels_csv, 'FileType','text','Delimiter',',','VariableNamingRule','preserve');
gene_ids  = string(LT{:,1});
cell_names = LT.Properties.VariableNames(2:end);
label_mat  = LT{:, 2:end};
map_csv = fullfile(repo,'results','qc','gene_map_ensdarg_to_model.csv');
map = readtable(map_csv, 'FileType','text','Delimiter',',','VariableNamingRule','preserve');
map_model = string(map.model_gene);
map_ens   = string(map.ensembl_gene_id);
ref_model = parents_by_hpf.h24;
n_gene = numel(ref_model.genes);
score_row_map = containers.Map(cellstr(gene_ids), num2cell(1:numel(gene_ids)));
gene_to_score_rows = cell(n_gene, 1);
for i = 1:size(map,1)
    mg = char(map_model(i)); eg = char(map_ens(i));
    idx = find(strcmp(ref_model.genes, mg), 1);
    if isempty(idx), continue; end
    if isKey(score_row_map, eg)
        gene_to_score_rows{idx}(end+1) = score_row_map(eg);
    end
end
fprintf('model genes with any measured paralog: %d / %d\n\n', ...
        sum(~cellfun('isempty',gene_to_score_rows)), n_gene);

% ---- job list (15 cells) -----------------------------------------------
conds = {'BL','D','LD'}; tps = [24 48 72 96 120];
jobs = cell(0,2);
for cc = 1:numel(conds), for tt = 1:numel(tps)
    jobs(end+1,:) = { conds{cc}, tps(tt) }; %#ok<AGROW>
end, end
n_jobs = size(jobs,1);
job_parent = cell(n_jobs, 1);
job_anchor = zeros(n_jobs, 1);
for ji = 1:n_jobs
    fld = sprintf('h%d', jobs{ji,2});
    job_parent{ji} = parents_by_hpf.(fld);
    job_anchor(ji) = anchors_by_hpf.(fld);
end

% ---- parpool (8 workers, Gurobi Threads=1) -----------------------------
pool = gcp('nocreate');
if isempty(pool)
    try
        cluster = parcluster('local');
        pool = parpool('local', min([N_WORKERS, n_jobs, cluster.NumWorkers]));
    catch, pool = []; end
end
if ~isempty(pool)
    fprintf('pool: %d workers\n\n', pool.NumWorkers);
    helper_dir = fullfile(repo,'functions');
    pctRunOnAll(sprintf('addpath(''%s'');', strrep(helper_dir,'''','''''')));
end

out = cell(n_jobs, 1);
wall_start = tic;

parfor ji = 1:n_jobs
    if exist('ensure_solver','file') == 2, ensure_solver(); end
    % Force Gurobi Threads=1 per worker (avoid oversubscription with 8 workers)
    try
        changeCobraSolverParams('LP','Threads',1);
        changeCobraSolverParams('MILP','Threads',1);
    catch, end

    cond = jobs{ji,1}; hpf = jobs{ji,2};
    tag = sprintf('%s_%d', cond, hpf);
    t_c = tic;
    fprintf('[%2d/%2d] %s ...\n', ji, n_jobs, tag);

    k_col = find(strcmp(cell_names, tag));
    ens_labels = double(label_mat(:, k_col));
    model_gene_label = zeros(n_gene, 1);
    for g = 1:n_gene
        rows = gene_to_score_rows{g};
        if isempty(rows), continue; end
        vals = ens_labels(rows);
        vals = vals(~isnan(vals));
        if ~isempty(vals), model_gene_label(g) = max(vals); end
    end

    m = job_parent{ji}; BIOMASS_ANCHOR = job_anchor(ji);
    n_r = numel(m.rxns);
    rxn_label = zeros(n_r, 1);
    for r = 1:n_r
        rule = m.rules{r};
        if isempty(rule) || all(isspace(rule)), continue; end
        rxn_label(r) = evalGPR_discrete(rule, model_gene_label);
    end
    for k = 1:numel(L1), j = findRxnIDs(m, L1{k}); if j>0, rxn_label(j)=1; end, end
    for k = 1:numel(L2), j = findRxnIDs(m, L2{k}); if j>0, rxn_label(j)=1; end, end
    for k = 1:numel(L3), j = findRxnIDs(m, L3{k}); if j>0, rxn_label(j)=1; end, end
    for k = 1:numel(L4), j = findRxnIDs(m, L4{k}); if j>0, rxn_label(j)=1; end, end
    expressed_idx = find(rxn_label == 1);

    stat = 'ok'; err_msg = '';
    n_active = NaN; bio_max = NaN; bio_feasible = false;
    try
        j_bio = findRxnIDs(m, BIOMASS_RXN);
        m_prelim = m; m_prelim.lb(j_bio) = BIOMASS_ANCHOR;
        [~, A_biomass] = run_single_core_fastcore(m_prelim, j_bio, EPSILON);

        consistent_idx = fastcc(m, EPSILON);
        if numel(consistent_idx) < numel(m.rxns)
            consistent = removeRxns(m, m.rxns(setdiff(1:numel(m.rxns), consistent_idx)));
        else, consistent = m; end
        core_ids_all = unique([L1(:); L2(:); L3(:); L4(:); m.rxns(expressed_idx); m.rxns(A_biomass)]);
        core_ids = core_ids_all(ismember(core_ids_all, consistent.rxns));
        map_c = containers.Map(consistent.rxns, num2cell(1:numel(consistent.rxns)));
        core_idx = cell2mat(values(map_c, core_ids));
        extracted = fastcore(consistent, core_idx, EPSILON);
        n_active = numel(extracted.rxns);
        j_bio_e = findRxnIDs(extracted, BIOMASS_RXN);
        if j_bio_e > 0
            sb = optimizeCbModel(changeObjective(extracted, BIOMASS_RXN), 'max');
            bio_max = sb.f;
            bio_feasible = (sb.stat == 1) && (sb.f > 1e-6);
        end
        parent = extracted; %#ok<NASGU>
        parent_meta = struct('cond',cond,'hpf',hpf,'method','trans_rfastcormics_v7', ...
            'epsilon',EPSILON,'n_active',n_active,'biomass_max',bio_max, ...
            'bio_feasible',bio_feasible, 'base_model','baked_ocr_*_v7.mat', ...
            'curation','v7_ledger_part1_part2');
        parsave(fullfile(pdir_out, sprintf('trans_rfastcormics_%s.mat', tag)), ...
                parent, parent_meta);
    catch ME
        stat = 'error'; err_msg = ME.message;
        fprintf('  [%s] ERROR: %s\n', tag, err_msg);
    end
    out{ji} = struct('tag',tag,'cond',cond,'hpf',hpf,'stat',stat, ...
        'n_active',n_active,'biomass_max',bio_max,'bio_feasible',bio_feasible, ...
        'err_msg',err_msg,'wall_sec',toc(t_c));
    fprintf('  [%s] n_active=%d  biomass_max=%.4g  wall=%.1fs\n', ...
            tag, n_active, bio_max, toc(t_c));
end

fprintf('\nTOTAL WALL: %.1f s\n\n', toc(wall_start));

% ---- summary table -----------------------------------------------------
rows = cell(n_jobs, 5);
for ji = 1:n_jobs
    r = out{ji};
    rows(ji,:) = { r.tag, r.cond, r.hpf, r.n_active, r.biomass_max };
end
T = cell2table(rows, 'VariableNames', {'cell','cond','hpf','n_active','biomass_max'});
writetable(T, fullfile(pdir_out, 's140_v7_biomass_max.csv'));
fprintf('wrote %s\n', fullfile(pdir_out,'s140_v7_biomass_max.csv'));

% comparison against v6 packet biomass_max values
packet = struct('BL_24',9.5238,'BL_48',9.3945,'BL_72',7.7832,'BL_96',8.6762,'BL_120',8.5619, ...
                'D_24',9.5238,'D_48',9.4008,'D_72',7.0797,'D_96',8.7835,'D_120',8.8956, ...
                'LD_24',9.5238,'LD_48',9.4952,'LD_72',8.5950,'LD_96',9.8116,'LD_120',8.7725);
fprintf('\n==== biomass_max: v6 packet vs v7 ====\n');
fprintf('%-8s %10s %10s %10s\n','cell','v6','v7','delta');
for ji = 1:n_jobs
    r = out{ji};
    p = packet.(r.tag);
    delta = r.biomass_max - p;
    fprintf('%-8s %10.4f %10.4f %+10.4f\n', r.tag, p, r.biomass_max, delta);
end

clear cleanupObj

% ---- helper functions (identical to s125) ------------------------------
function parsave(fp, parent, parent_meta)
    save(fp, 'parent','parent_meta','-v7.3');
end
function [extracted, A_biomass_idx] = run_single_core_fastcore(model, core_idx, epsilon)
    consistent_idx = fastcc(model, epsilon);
    if numel(consistent_idx) < numel(model.rxns)
        consistent = removeRxns(model, model.rxns(setdiff(1:numel(model.rxns), consistent_idx)));
    else, consistent = model; end
    map_c = containers.Map(consistent.rxns, num2cell(1:numel(consistent.rxns)));
    biomass_rxn_id = model.rxns{core_idx};
    core_in_consistent = map_c(biomass_rxn_id);
    extracted = fastcore(consistent, core_in_consistent, epsilon);
    A_biomass_idx = zeros(numel(extracted.rxns), 1);
    for k = 1:numel(extracted.rxns)
        j = find(strcmp(model.rxns, extracted.rxns{k}), 1);
        if ~isempty(j), A_biomass_idx(k) = j; end
    end
    A_biomass_idx = A_biomass_idx(A_biomass_idx > 0);
end
function v = evalGPR_discrete(rule, gene_label)
    rule = strtrim(rule);
    if isempty(rule), v = 0; return; end
    [v, ~] = parse_expr(rule, gene_label);
end
function [v, rest] = parse_expr(s, ge)
    [v, s] = parse_term(s, ge); s = strtrim(s);
    while ~isempty(s) && s(1) == '|'
        s = strtrim(s(2:end)); [v2, s] = parse_term(s, ge); v = max(v, v2); s = strtrim(s);
    end
    rest = s;
end
function [v, rest] = parse_term(s, ge)
    [v, s] = parse_factor(s, ge); s = strtrim(s);
    while ~isempty(s) && s(1) == '&'
        s = strtrim(s(2:end)); [v2, s] = parse_factor(s, ge); v = min(v, v2); s = strtrim(s);
    end
    rest = s;
end
function [v, rest] = parse_factor(s, ge)
    s = strtrim(s);
    if isempty(s), v = 0; rest = ''; return; end
    if s(1) == '('
        [v, s] = parse_expr(s(2:end), ge); s = strtrim(s);
        if ~isempty(s) && s(1) == ')', s = s(2:end); end
        rest = s; return;
    end
    tok = regexp(s, '^x\((\d+)\)', 'tokens', 'once');
    if isempty(tok), v = 0; rest = s(2:end); return; end
    idx = str2double(tok{1});
    if idx >= 1 && idx <= numel(ge)
        val = ge(idx); if isnan(val), v = 0; else, v = val; end
    else, v = 0; end
    rest = regexprep(s, '^x\(\d+\)', '', 'once');
end
