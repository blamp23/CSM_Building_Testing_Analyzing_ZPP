% s135 -- v7 base build with the FULL boundary rule.
% Supersedes s131. Applies, in order:
%   (a) MAR00022 pool coef edits from s130_pool_table.csv
%       (4 reduced-redox substrates -> coef 0)
%   (b) close_exchange rows from s130_v7_ledger.csv  (lb = ub = 0)
%   (c) Part 1 (secretion_uncapped): for every exchange NOT in the 105 AND
%       NOT on the close list AND NOT on the medium list -> ub = 1000.
%   (d) Part 2 (medium_uptake + unmeasured-organic-close):
%       - Medium list (hardcoded 20 MAR IDs of bulk inorganics/water/gases)
%         gets lb = -1000, ub = 1000.
%       - Everything else NOT in the 105 AND NOT on the close/keep list
%         gets lb = 0 (Part 2: close unmeasured organic uptake).
%       - Class-1 keep_exchange rows are left at their current bounds
%         (fc applied at extraction time).
%   (e) O2 (MAR09048) is left untouched (OCR bake).
%
% Saves as baked_ocr_72_v7.mat (overwrites the s131 output).

REPO      = '/Users/lamp_b/Library/CloudStorage/OneDrive-TexasA&MUniversity/Hala, David''s files - Benji_COBRA/Tanguay_Data/Discrete_Models';
LEDGER    = fullfile(REPO, 'reviewer_packet_v4_slim/s130_v7_ledger.csv');
POOL      = fullfile(REPO, 'reviewer_packet_v4_slim/s130_pool_table.csv');
BOUNDS    = fullfile(REPO, 'results/qc/phase2_exchange_bounds_long.csv');

% Loop over all 5 hpf baselines (each has its own OCR cap baked in).
% Produces baked_ocr_{24,48,72,96,120}_v7.mat.
hpf_list = [24 48 72 96 120];
for hpf = hpf_list
    IN_MAT  = fullfile(REPO, sprintf('ocr_anchored_extraction/models/baked_ocr_%d.mat', hpf));
    OUT_MAT = fullfile(REPO, sprintf('ocr_anchored_extraction/models/baked_ocr_%d_v7.mat', hpf));
    fprintf('\n======= hpf %d =======\n', hpf);
    fprintf('loading %s\n', IN_MAT);
    L = load(IN_MAT); m = L.m;
    fprintf('base: %d rxns, %d mets\n', numel(m.rxns), numel(m.mets));

% ---- 105 measured set ---------------------------------------------------
mb = readtable(BOUNDS);
measured = unique(mb.ex_rxn);
fprintf('metabolomics 105: %d exchange MAR IDs\n', numel(measured));

% ---- Ledger classification ---------------------------------------------
ledger = readtable(LEDGER);
close_list = ledger.exchange_MAR_id(strcmp(ledger.action, 'close_exchange'));
close_list = close_list(~cellfun('isempty', close_list));
% keep-open actions: keep_exchange (class 1 default), medium_uptake (class 1
% vitamin/carotenoid additions from s137), ysl_lipid_delivery (class 1
% yolk syncytial lipoprotein additions from s137).
keep_mask  = strcmp(ledger.action, 'keep_exchange') | ...
             strcmp(ledger.action, 'medium_uptake') | ...
             strcmp(ledger.action, 'ysl_lipid_delivery');
keep_list  = ledger.exchange_MAR_id(keep_mask);
keep_list  = keep_list(~cellfun('isempty', keep_list));
fprintf('ledger: %d close, %d keep-open (keep_exchange + medium_uptake + ysl_lipid_delivery)\n', ...
        numel(close_list), numel(keep_list));

% ---- Medium list (bulk inorganics + water + gases + protected O2) -------
% Whitelisted at lb = -1000, ub = 1000. O2 handled separately (OCR bake).
medium_list = { ...
    'MAR09047', ... % H2O
    'MAR09058', ... % CO2
    'MAR09072', ... % Pi
    'MAR09073', ... % NH3
    'MAR09074', ... % sulfate
    'MAR09076', ... % Fe2+
    'MAR09077', ... % Na+
    'MAR09078', ... % HCO3-
    'MAR09079', ... % H+
    'MAR09080', ... % Zn2+
    'MAR09081', ... % K+
    'MAR09082', ... % Ca2+
    'MAR09148', ... % iodide
    'MAR09150', ... % Cl-
    'MAR13066', ... % carbonate
    'MAR13072', ... % Mg2+
    'MAR13073', ... % Cu2+
};
O2_EX = 'MAR09048';

% ---- (a) pool edits -----------------------------------------------------
% (a1) MAR00022 substrate zeros for reduced-redox states
poolTbl = readtable(POOL);
j_pool = find(strcmp(m.rxns, 'MAR00022'));
fprintf('\n(a1) editing MAR00022 pool coefficients ...\n');
zeroed = 0;
for i = 1:height(poolTbl)
    mm  = poolTbl.MAM_id{i};
    nc  = poolTbl.new_coef(i);
    idx = find(strcmp(m.mets, mm));
    if isempty(idx), continue; end
    old = full(m.S(idx, j_pool));
    if abs(nc) < 1e-12 && abs(old) > 0
        m.S(idx, j_pool) = 0;
        fprintf('    zeroed %s (%s): %+g -> 0\n', mm, poolTbl.metabolite{i}, old);
        zeroed = zeroed + 1;
    end
end

% (a2) MAR00021 pool pseudo-metabolite coef -> -EPSILON (amendment 2026-09-30).
% Wang's -1 is a placeholder; ε=1e-3 lets the pool track trace consumption
% without inflating it into a growth ceiling. Coefficient parameterized via
% the EPSILON variable so the {1e-2, 1e-3, 1e-4} scan is a two-line change.
if ~exist('EPSILON','var'), EPSILON = 1e-3; end
j_bio = find(strcmp(m.rxns, 'MAR00021'));
i_pool_met = find(strcmp(m.mets, 'MAM01602c'));
if ~isempty(j_bio) && ~isempty(i_pool_met)
    old_coef = full(m.S(i_pool_met, j_bio));
    m.S(i_pool_met, j_bio) = -EPSILON;
    fprintf('\n(a2) MAR00021 pool coef (MAM01602c): %+g -> %+g (EPSILON=%g)\n', old_coef, -EPSILON, EPSILON);
end

% ---- (b) close_exchange list -------------------------------------------
fprintf('\n(b) closing %d ledger exchanges ...\n', numel(close_list));
closed = 0;
for i = 1:numel(close_list)
    rid = close_list{i};
    j = find(strcmp(m.rxns, rid));
    if isempty(j), warning('  %s not found', rid); continue; end
    m.lb(j) = 0; m.ub(j) = 0;
    closed = closed + 1;
end
fprintf('    closed %d exchanges\n', closed);

% ---- (c) + (d) enumerate ALL exchanges and classify --------------------
% Determine is_ex via S connectivity (single-nonzero column)
n_rxn = numel(m.rxns);
n_mets_per_rxn = zeros(1, n_rxn);
for j = 1:n_rxn
    n_mets_per_rxn(j) = nnz(m.S(:, j));
end
is_ex = n_mets_per_rxn == 1;
ex_js = find(is_ex);
fprintf('\n(c+d) applying Part 1 + Part 2 to %d exchange reactions ...\n', numel(ex_js));

n_secretion_uncapped = 0;
n_medium = 0;
n_uptake_closed_part2 = 0;
n_kept_measured = 0;
n_kept_keeplist = 0;
n_kept_o2 = 0;
n_kept_close = 0;

for k = 1:numel(ex_js)
    j = ex_js(k);
    rid = m.rxns{j};
    if strcmp(rid, O2_EX)
        n_kept_o2 = n_kept_o2 + 1; continue;
    end
    if any(strcmp(rid, close_list))
        n_kept_close = n_kept_close + 1; continue;
    end
    if any(strcmp(rid, measured))
        n_kept_measured = n_kept_measured + 1; continue;
    end
    if any(strcmp(rid, medium_list))
        m.lb(j) = -1000; m.ub(j) = 1000;
        n_medium = n_medium + 1; continue;
    end
    if any(strcmp(rid, keep_list))
        % class-1 vitamin/AA at default; leave as-is (fc applied later)
        n_kept_keeplist = n_kept_keeplist + 1; continue;
    end
    % All remaining unmeasured non-close non-keep: Part 1 uncap secretion,
    % Part 2 close uptake. v7 base biomass_max = 9.839 vs v6 10.0; the 1.6%
    % gap is diffuse import of unmeasured organics that the ±10 default
    % permitted; an embryo doesn't have those uptakes, so the lower number
    % is the more honest one. Accepted 2026-09-29.
    m.ub(j) = 1000;
    m.lb(j) = 0;
    n_secretion_uncapped = n_secretion_uncapped + 1;
    n_uptake_closed_part2 = n_uptake_closed_part2 + 1;
end

fprintf('\n  Part 1 + Part 2 summary:\n');
fprintf('    kept (O2, OCR bake)                : %d\n', n_kept_o2);
fprintf('    kept (ledger close_exchange)       : %d\n', n_kept_close);
fprintf('    kept (in 105 measured set)         : %d\n', n_kept_measured);
fprintf('    kept (ledger keep_exchange class 1): %d\n', n_kept_keeplist);
fprintf('    set to medium (lb=-1000, ub=1000)  : %d\n', n_medium);
fprintf('    Part 1+2 (lb=0, ub=1000)           : %d\n', n_secretion_uncapped);
fprintf('    (total exchanges)                  : %d\n', numel(ex_js));

% ---- save --------------------------------------------------------------
save(OUT_MAT, 'm', '-v7.3');
fprintf('\nwrote %s\n', OUT_MAT);
end   % for hpf
