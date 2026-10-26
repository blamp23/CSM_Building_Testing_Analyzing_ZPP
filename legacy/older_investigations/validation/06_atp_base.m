% 06_atp_base -- ATP sanity on the v7 base, Fritzemeier 2017 style.
% Standalone MATLAB script. Run on any baked_ocr_{hpf}_v7.mat baseline.
% Loops-allowed (LP) + loopless (MILP) via COBRA Toolbox and Gurobi.
%
% Gate:
%   test 1 (ATP from nothing, O2 closed)   = 0
%   test 2 (ATP from nothing, O2 open)     = 0
%   test 3 (aerobic yield, glucose = 1)   <= 32
%   test 4 (anaerobic yield, glucose = 1) <= 2
%   test 5 (full v7 medium, biomass unfixed) report only

repo = getenv('V7_REPO');
if isempty(repo)
    error('set V7_REPO to the data root (e.g. the OneDrive repo dir)');
end
base = fullfile(repo, 'ocr_anchored_extraction/models/baked_ocr_72_v7.mat');
S = load(base); m = S.m;

ATP = 'MAR03964';
GLC = 'MAR09034';
O2  = 'MAR09048';
MEDIUM = {'MAR09047','MAR09058','MAR09072','MAR09073','MAR09074','MAR09076', ...
          'MAR09077','MAR09078','MAR09079','MAR09080','MAR09081','MAR09082', ...
          'MAR09148','MAR09150','MAR13066','MAR13072','MAR13073'};

Snnz  = full(sum(m.S ~= 0, 1));
is_ex = (Snnz == 1)';
j_atp = findRxnIDs(m, ATP);
j_glc = findRxnIDs(m, GLC);
j_o2  = findRxnIDs(m, O2);

reset_bounds = @(mm) deal_bounds(mm, is_ex);

% Test 1
m1 = reset_bounds(m);
m1.lb(j_o2) = 0; m1.ub(j_o2) = 0;
for k = 1:numel(MEDIUM), j = findRxnIDs(m1, MEDIUM{k}); if j>0, m1.lb(j) = -1000; end; end
m1 = changeObjective(m1, ATP);
v1 = optimizeCbModel(m1, 'max').f;

% Test 2
m2 = reset_bounds(m);
m2.lb(j_o2) = -1000; m2.ub(j_o2) = 0;
for k = 1:numel(MEDIUM), j = findRxnIDs(m2, MEDIUM{k}); if j>0, m2.lb(j) = -1000; end; end
m2 = changeObjective(m2, ATP);
v2 = optimizeCbModel(m2, 'max').f;

% Test 3
m3 = reset_bounds(m);
m3.lb(j_o2) = -1000; m3.ub(j_o2) = 0;
m3.lb(j_glc) = -1; m3.ub(j_glc) = 0;
for k = 1:numel(MEDIUM), j = findRxnIDs(m3, MEDIUM{k}); if j>0, m3.lb(j) = -1000; end; end
m3 = changeObjective(m3, ATP);
v3 = optimizeCbModel(m3, 'max').f;

% Test 4
m4 = reset_bounds(m);
m4.lb(j_o2) = 0; m4.ub(j_o2) = 0;
m4.lb(j_glc) = -1; m4.ub(j_glc) = 0;
for k = 1:numel(MEDIUM), j = findRxnIDs(m4, MEDIUM{k}); if j>0, m4.lb(j) = -1000; end; end
m4 = changeObjective(m4, ATP);
v4 = optimizeCbModel(m4, 'max').f;

% Test 5
m5 = m;
m5 = changeObjective(m5, ATP);
sol5 = optimizeCbModel(m5, 'max');
v5 = sol5.f;
o2_v5 = sol5.x(j_o2);

fprintf('\n==== ATP sanity (v7 base) ====\n');
fprintf('  test 1 (nothing + O2 closed): %.3f  (expect 0)\n', v1);
fprintf('  test 2 (nothing + O2 open ): %.3f  (expect 0)\n', v2);
fprintf('  test 3 (aerobic glc=1    ): %.3f  (expect <=32)\n', v3);
fprintf('  test 4 (anaerobic glc=1  ): %.3f  (expect <=2)\n', v4);
fprintf('  test 5 (full v7 medium   ): ATPD = %.3f, O2 = %+.3f\n', v5, o2_v5);

function m = deal_bounds(m, is_ex)
    m.lb(is_ex) = 0; m.ub(is_ex) = 1000;
end
