function model = ensureCobra(model)
% ensureCobra  Add COBRA-required fields missing from a RAVEN-style struct.
%
% Lifted from Data/Biomass_Search/interrogate_biomass.m. This is the ONLY place
% the shim should live going forward. Anything in the pipeline that touches the
% model must run it through this first.

n = numel(model.rxns);
m = numel(model.mets);

if ~isfield(model,'rules') || isempty(model.rules)
    if isfield(model,'grRules')
        try
            model = generateRules(model);
            fprintf('ensureCobra: generated model.rules from grRules\n');
        catch
            model.rules = repmat({''}, n, 1);
            fprintf('ensureCobra: generateRules failed, rules left empty\n');
        end
    else
        model.rules = repmat({''}, n, 1);
    end
end
if ~isfield(model,'osenseStr'), model.osenseStr = 'max'; end
if ~isfield(model,'csense') || numel(model.csense) ~= m
    model.csense = repmat('E', m, 1);
end
if ~isfield(model,'b') || numel(model.b) ~= m
    model.b = zeros(m,1);
end
if isfield(model,'metCharges')
    model.metCharges = double(model.metCharges);
end
if ~isfield(model,'metCharge') && isfield(model,'metCharges')
    model.metCharge = model.metCharges;
end
if ~isfield(model,'rxnGeneMat') && isfield(model,'grRules')
    model = buildRxnGeneMat(model);
end

% verifyModel is chatty about the subSystems format on RAVEN-style structs;
% swallow its output but keep the pass/fail bool.
try
    [~, res] = evalc('verifyModel(model, ''simpleCheck'', true)');
    fprintf('ensureCobra: verifyModel simpleCheck = %d\n', res);
catch
    fprintf('ensureCobra: verifyModel skipped (raised)\n');
end
end
