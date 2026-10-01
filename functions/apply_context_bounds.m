function [m, fc_applied, n_p12] = apply_context_bounds(m, cond, hpf, frame, matched_biomax)
% apply_context_bounds -- canonical per-cell bound application.
%
% Starts from the model's stored bounds (v7 boundary: Part 1 secretion
% uncapped, Part 2 unmeasured uptake closed, ledger closes, medium open).
% NEVER resets exchanges to a uniform default. Applies:
%
%   (1) BOUND_EX * fc on the 105 measured exchanges
%   (2) OCR cap on MAR09048
%   (3) Frame A: biomass_lb = 0.9 * biomass_max computed here;
%       'matched' with matched_biomax = per-hpf minimum biomax
%
% Asserts:
%   - every ledger close_exchange row has lb == ub == 0
%   - every unmeasured non-medium exchange has lb == 0 && ub == 1000
%   - every medium row has lb < 0
%   - exactly (# measured present in model) exchanges carry fc scaling

if nargin < 4, frame = 'none'; end
if nargin < 5, matched_biomax = NaN; end

REPO = '/Users/lamp_b/Library/CloudStorage/OneDrive-TexasA&MUniversity/Hala, David''s files - Benji_COBRA/Tanguay_Data/Discrete_Models';
BOUNDS_CSV = fullfile(REPO,'results','qc','phase2_exchange_bounds_long.csv');
LEDGER_CSV = fullfile(REPO,'reviewer_packet_v4_slim','s130_v7_ledger.csv');
BOUND_EX   = 10.0;   % relative-envelope convention; fc-multiplied per cell
O2_EX      = 'MAR09048';
BIOMASS    = 'MAR00021';
OCR = containers.Map({24,48,72,96,120}, {97,173,227,273,313});
MEDIUM = {'MAR09047','MAR09058','MAR09072','MAR09073','MAR09074','MAR09076', ...
          'MAR09077','MAR09078','MAR09079','MAR09080','MAR09081','MAR09082', ...
          'MAR09148','MAR09150','MAR13066','MAR13072','MAR13073'};

mb = readtable(BOUNDS_CSV,'FileType','text','Delimiter',',','VariableNamingRule','preserve');
if exist(LEDGER_CSV,'file')
    led = readtable(LEDGER_CSV,'FileType','text','Delimiter',',','VariableNamingRule','preserve');
    close_list = led.exchange_MAR_id(strcmp(led.action,'close_exchange'));
    close_list = close_list(~cellfun('isempty',close_list));
    keep_mask = strcmp(led.action,'keep_exchange') | ...
                strcmp(led.action,'medium_uptake') | ...
                strcmp(led.action,'ysl_lipid_delivery');
    keep_list = led.exchange_MAR_id(keep_mask);
    keep_list = keep_list(~cellfun('isempty',keep_list));
    medium_uptake_list = led.exchange_MAR_id(strcmp(led.action,'medium_uptake'));
    medium_uptake_list = medium_uptake_list(~cellfun('isempty',medium_uptake_list));
else
    close_list = {}; keep_list = {}; medium_uptake_list = {};
end

% detect exchange reactions (single-nonzero column)
n_rxn = numel(m.rxns);
is_ex = false(n_rxn,1);
for j = 1:n_rxn
    if nnz(m.S(:,j)) == 1, is_ex(j) = true; end
end

% -------- (1) fc scaling on 105 measured --------------------------------
sel = strcmp(mb.condition, cond) & (mb.hpf == hpf);
ex_list  = mb.ex_rxn(sel);
fc_list  = mb.fc(sel);
measured_all = unique(mb.ex_rxn);
fc_applied = 0;
for k = 1:numel(ex_list)
    j = findRxnIDs(m, ex_list{k});
    if j > 0
        fc = fc_list(k);
        if ~isfinite(fc) || fc <= 0, fc = 1.0; end
        m.lb(j) = -BOUND_EX * fc;
        m.ub(j) = +BOUND_EX * fc;
        fc_applied = fc_applied + 1;
    end
end

% -------- (2) OCR cap ---------------------------------------------------
j_o2 = findRxnIDs(m, O2_EX);
if j_o2 > 0
    m.lb(j_o2) = -OCR(hpf);
    m.ub(j_o2) = 0;
end

% -------- (2b) Free supply for unmeasured class-1 medium_uptake ---------
% 2026-10-01: medium_uptake rows not in the 105 get lb=-1000.
for k = 1:numel(medium_uptake_list)
    rid = medium_uptake_list{k};
    if any(strcmp(measured_all, rid)), continue; end
    j = findRxnIDs(m, rid);
    if j > 0, m.lb(j) = -1000; end
end

% -------- assertions ----------------------------------------------------
% a) ledger closes
for k = 1:numel(close_list)
    j = findRxnIDs(m, close_list{k});
    if j > 0
        assert(m.lb(j) == 0 && m.ub(j) == 0, ...
               'assertion: ledger close %s not closed [lb=%g, ub=%g]', ...
               close_list{k}, m.lb(j), m.ub(j));
    end
end
% b) Part 1+2 rows (unmeasured, not on any ledger row)
n_p12 = 0;
for jj = find(is_ex)'
    rid = m.rxns{jj};
    if strcmp(rid, O2_EX), continue; end
    if any(strcmp(close_list, rid)), continue; end
    if any(strcmp(keep_list, rid)), continue; end
    if any(strcmp(MEDIUM, rid)), continue; end
    if any(strcmp(measured_all, rid)), continue; end
    assert(m.lb(jj) == 0 && m.ub(jj) == 1000, ...
           'assertion: Part 1+2 row %s bounds violated [lb=%g, ub=%g]', ...
           rid, m.lb(jj), m.ub(jj));
    n_p12 = n_p12 + 1;
end
% c) medium rows
for k = 1:numel(MEDIUM)
    j = findRxnIDs(m, MEDIUM{k});
    if j > 0
        assert(m.lb(j) < 0, 'assertion: medium %s lb not < 0 [lb=%g]', MEDIUM{k}, m.lb(j));
    end
end
% d) exactly (# measured present) got fc
n_present = 0;
for k = 1:numel(measured_all)
    if findRxnIDs(m, measured_all{k}) > 0, n_present = n_present + 1; end
end
assert(fc_applied == n_present, ...
       'assertion: fc applied to %d, expected %d present in model', fc_applied, n_present);

% -------- (3) biomass floor --------------------------------------------
j_bio = findRxnIDs(m, BIOMASS);
if strcmpi(frame, 'A')
    m_probe = m;
    sol = optimizeCbModel(changeObjective(m_probe, BIOMASS), 'max');
    if sol.stat == 1 && sol.f > 1e-9
        m.lb(j_bio) = 0.9 * sol.f;
    end
elseif strcmpi(frame, 'matched')
    if isfinite(matched_biomax) && matched_biomax > 0
        m.lb(j_bio) = 0.9 * matched_biomax;
    end
end
end
