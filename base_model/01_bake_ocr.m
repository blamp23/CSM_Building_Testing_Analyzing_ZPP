% s01 -- bake per-hpf parent models with OCR-realistic O2 bounds.
%
% Alternative to bounds10_extraction/s01_bake_base.m.
% Same as that script EXCEPT: the O2 exchange (MAR09048) is bounded by
% measured OCR per timepoint rather than the uniform +-10.
%
% Produces FIVE parent models (one per hpf), each with:
%   - +-10 on every exchange EXCEPT O2 (MAR09048)
%   - O2 upper bound = -OCR_umol_O2_gDW_hr (uptake is negative in Wang convention)
%   - O2 lower bound = 0 (secretion blocked)
%   - o2_leak_fixes applied
%   - cofactor pool modification already in the base model
%
% OCR values (from Data/Metabolomics/physiology_by_timepoint.csv):
%     24 hpf:  96.79 umol/gDW/hr
%     48 hpf: 173.35
%     72 hpf: 226.63
%     96 hpf: 273.07
%    120 hpf: 313.00
%
% Save: ocr_anchored_extraction/models/baked_ocr_{24,48,72,96,120}.mat
%
% Run from repo root:
%   addpath('functions'); initCobraToolbox(false);
%   changeCobraSolver('gurobi','LP');
%   run('ocr_anchored_extraction/s01_bake_base_ocr.m')

here = fileparts(mfilename('fullpath'));
repo = fileparts(here);
addpath(fullfile(repo,'functions'));

diary_dir = fullfile(here,'results');
if ~exist(diary_dir,'dir'), mkdir(diary_dir); end
diary(fullfile(diary_dir,'s01_bake_base_ocr.txt')); diary on
cleanupObj = onCleanup(@() diary('off'));

fprintf('=== s01 bake base OCR-anchored ===\n');
fprintf('date: %s\n\n', datestr(now,'yyyy-mm-dd HH:MM:SS'));

% ---- load base ----------------------------------------------------------
base_path = fullfile(repo,'models','zebrafishGEM_v2_modcofpool.mat');
assert(exist(base_path,'file')==2, 'missing base model: %s', base_path);
L = load(base_path);
m0 = find_model_struct(L);
fprintf('loaded base: %d rxns, %d mets, %d genes\n', ...
        numel(m0.rxns), numel(m0.mets), numel(m0.genes));

% ---- apply O2 / H2O2 leak fixes ----------------------------------------
fix_csv = fullfile(repo,'gimme_transcriptomics_pipeline','results', ...
                   'pipeline','o2_leak_fixes.csv');
assert(exist(fix_csv,'file')==2, 'missing %s', fix_csv);
FIXES = readtable(fix_csv, 'FileType','text', 'Delimiter',',', ...
                  'VariableNamingRule','preserve');
fprintf('loaded %d O2/H2O2 leak fixes from %s\n', height(FIXES), fix_csv);

n_applied = 0; n_missing = 0;
for fi = 1:height(FIXES)
    rid = char(FIXES.rxn_id{fi});
    act = char(FIXES.action{fi});
    j   = findRxnIDs(m0, rid);
    if j == 0, n_missing = n_missing + 1; continue; end
    switch act
        case 'block_fwd'
            m0.ub(j) = min(m0.ub(j), 0);
        case 'block_rev'
            m0.lb(j) = max(m0.lb(j), 0);
        case {'block_both','close_both'}
            m0.lb(j) = 0; m0.ub(j) = 0;
        otherwise
            warning('unknown action "%s" for %s -- skipping', act, rid);
            continue;
    end
    n_applied = n_applied + 1;
end
fprintf('  applied %d fixes  (missing rxn: %d)\n', n_applied, n_missing);

% ---- uniform +-10 on every exchange ------------------------------------
is_ex = false(numel(m0.rxns),1);
for r = 1:numel(m0.rxns)
    col = m0.S(:, r);
    nz = find(col ~= 0);
    if numel(nz) == 1
        is_ex(r) = true;
    end
end
ex_idx = find(is_ex);
fprintf('exchange rxns detected: %d\n', numel(ex_idx));

BOUND = 10;
m0.lb(ex_idx) = -BOUND;
m0.ub(ex_idx) =  BOUND;
fprintf('set lb=-%g, ub=+%g on all %d exchanges (before OCR override)\n', ...
        BOUND, BOUND, numel(ex_idx));

% ---- load OCR bounds ---------------------------------------------------
phys_csv = fullfile(repo,'Data','Metabolomics','physiology_by_timepoint.csv');
assert(exist(phys_csv,'file')==2, 'missing %s', phys_csv);
PHY = readtable(phys_csv, 'FileType','text','Delimiter',',',...
                'VariableNamingRule','preserve');
fprintf('loaded physiology table (%d rows)\n', height(PHY));

% Wang O2 exchange = MAR09048
O2_RXN = 'MAR09048';
j_o2 = findRxnIDs(m0, O2_RXN);
assert(j_o2 > 0, 'O2 exchange %s not in model', O2_RXN);
fprintf('O2 exchange rxn %s -> col %d\n', O2_RXN, j_o2);

% ---- build per-hpf parents ---------------------------------------------
hpf_list = [24 48 72 96 120];
summary  = struct('hpf',{}, 'ocr',{}, 'biomass_max',{}, 'stat',{}, 'path',{});

for h = 1:numel(hpf_list)
    hpf = hpf_list(h);
    ocr = PHY.ocr_umol_O2_gDW_hr(PHY.hpf == hpf);
    assert(~isempty(ocr), 'no OCR for hpf=%d', hpf);
    ocr = ocr(1);

    m = m0;
    % Wang / Human-GEM convention: uptake = negative flux for boundary rxn
    % So OCR upper cap on UPTAKE = lb of exchange (more negative).
    % Keep secretion at 0 (ub = 0) so O2 can only be consumed.
    m.lb(j_o2) = -ocr;
    m.ub(j_o2) =  0;

    j_bio = findRxnIDs(m, 'MAR00021');
    assert(j_bio > 0, 'MAR00021 (biomass) not in model');
    m_test = changeObjective(m, 'MAR00021');
    sol = optimizeCbModel(m_test, 'max');
    fprintf('  hpf=%3d  OCR=%7.2f  biomass_max=%8.4f  stat=%d\n', ...
            hpf, ocr, sol.f, sol.stat);

    baked_meta = struct();
    baked_meta.source           = 'zebrafishGEM_v2_modcofpool.mat';
    baked_meta.o2_h2o2_fix_csv  = fix_csv;
    baked_meta.n_o2_h2o2_fixes  = n_applied;
    baked_meta.exchange_bound   = BOUND;
    baked_meta.n_exchanges      = numel(ex_idx);
    baked_meta.hpf              = hpf;
    baked_meta.ocr_umol_O2_gDW_hr = ocr;
    baked_meta.o2_rxn           = O2_RXN;
    baked_meta.o2_lb            = -ocr;
    baked_meta.o2_ub            =  0;
    baked_meta.n_rxns           = numel(m.rxns);
    baked_meta.n_mets           = numel(m.mets);
    baked_meta.biomass_max_fba  = sol.f;
    baked_meta.date             = datestr(now,'yyyy-mm-dd HH:MM:SS');

    out_path = fullfile(here,'models', sprintf('baked_ocr_%d.mat', hpf));
    save(out_path, 'm', 'baked_meta', '-v7.3');
    fprintf('    wrote %s\n', out_path);

    summary(h).hpf = hpf;
    summary(h).ocr = ocr;
    summary(h).biomass_max = sol.f;
    summary(h).stat = sol.stat;
    summary(h).path = out_path;
end

% ---- summary CSV -------------------------------------------------------
T = table([summary.hpf]', [summary.ocr]', [summary.biomass_max]', ...
          [summary.stat]', 'VariableNames', ...
          {'hpf','ocr_umol_O2_gDW_hr','parent_biomass_max','solver_stat'});
writetable(T, fullfile(here,'results','baked_ocr_summary.csv'));
fprintf('\nwrote baked_ocr_summary.csv\n');
disp(T);

clear cleanupObj

function m = find_model_struct(L)
    fn = fieldnames(L);
    for i = 1:numel(fn)
        v = L.(fn{i});
        if isstruct(v) && isfield(v,'rxns'), m = v; return; end
    end
    error('no struct with .rxns field found in loaded mat');
end
