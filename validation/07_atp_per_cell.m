% s141 -- ATP tests 1 and 3 on each of the 15 v7-extracted models.
%
% Gate:
%   Test 1 (ATP from nothing, O2 closed) = 0 for every model.
%   Test 3 (aerobic yield, glc = 1) <= 31.5 for every model.
%
% Uses the same MATLAB+Gurobi ATP recipe you ran manually on the base
% (MAR03964, MAR09034, MAR09048).

REPO = '/Users/lamp_b/Library/CloudStorage/OneDrive-TexasA&MUniversity/Hala, David''s files - Benji_COBRA/Tanguay_Data/Discrete_Models';
MODELS_DIR = fullfile(REPO,'reviewer_packet_v4_slim','v7_models');
BOUNDS_CSV = fullfile(REPO,'results','qc','phase2_exchange_bounds_long.csv');

ATP = 'MAR03964'; GLC = 'MAR09034'; O2 = 'MAR09048';
MEDIUM = {'MAR09047','MAR09058','MAR09072','MAR09073','MAR09074','MAR09076', ...
          'MAR09077','MAR09078','MAR09079','MAR09080','MAR09081','MAR09082', ...
          'MAR09148','MAR09150','MAR13066','MAR13072','MAR13073'};

conds = {'BL','D','LD'}; tps = [24 48 72 96 120];
rows = {};
for cc = 1:numel(conds), for tt = 1:numel(tps)
    tag = sprintf('%s_%d', conds{cc}, tps(tt));
    fp = fullfile(MODELS_DIR, sprintf('trans_rfastcormics_%s.mat', tag));
    if ~exist(fp,'file'), fprintf('MISSING %s\n', tag); continue; end
    S = load(fp); m = S.parent;
    j_atp = findRxnIDs(m, ATP);
    j_glc = findRxnIDs(m, GLC);
    j_o2  = findRxnIDs(m, O2);
    if any([j_atp j_glc j_o2] == 0)
        fprintf('%s: missing key reaction (atp=%d glc=%d o2=%d)\n', tag, j_atp, j_glc, j_o2);
        continue
    end

    % Test 1: ATP from nothing, O2 closed
    m1 = m;
    n_r = numel(m1.rxns);
    is_ex = (sum(m1.S ~= 0, 1) == 1)';
    m1.lb(is_ex) = 0; m1.ub(is_ex) = 1000;
    for k = 1:numel(MEDIUM)
        j = findRxnIDs(m1, MEDIUM{k}); if j>0, m1.lb(j) = -1000; end
    end
    m1.lb(j_o2) = 0; m1.ub(j_o2) = 0;
    m1 = changeObjective(m1, ATP);
    s1 = optimizeCbModel(m1, 'max'); v1 = s1.f;

    % Test 3: aerobic yield, glc = 1
    m3 = m;
    m3.lb(is_ex) = 0; m3.ub(is_ex) = 1000;
    for k = 1:numel(MEDIUM)
        j = findRxnIDs(m3, MEDIUM{k}); if j>0, m3.lb(j) = -1000; end
    end
    m3.lb(j_o2) = -1000; m3.ub(j_o2) = 0;
    m3.lb(j_glc) = -1; m3.ub(j_glc) = 0;
    m3 = changeObjective(m3, ATP);
    s3 = optimizeCbModel(m3, 'max'); v3 = s3.f;

    gate1 = abs(v1) < 1e-6;
    gate3 = v3 <= 31.5 + 1e-4;
    rows(end+1,:) = { tag, v1, gate1, v3, gate3 }; %#ok<AGROW>
    fprintf('%s: test1=%.4f (%s)  test3=%.4f (%s)\n', tag, v1, ...
        ternary(gate1,'PASS','FAIL'), v3, ternary(gate3,'PASS','FAIL'));
end, end

T = cell2table(rows, 'VariableNames', {'cell','atp_from_nothing','test1_pass','aerobic_glc1','test3_pass'});
writetable(T, fullfile(REPO,'reviewer_packet_v4_slim','s141_v7_atp_per_cell.csv'));
fprintf('\nwrote s141_v7_atp_per_cell.csv\n');

fprintf('\n==== SUMMARY ====\n');
fprintf('cells passing test1 (=0):        %d/%d\n', sum([rows{:,3}]), size(rows,1));
fprintf('cells passing test3 (<=31.5):    %d/%d\n', sum([rows{:,5}]), size(rows,1));

function s = ternary(b, a, b_str)
    if b, s = a; else, s = b_str; end
end
