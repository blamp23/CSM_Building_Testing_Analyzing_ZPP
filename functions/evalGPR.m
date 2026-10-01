% Benji Lamp

function v = evalGPR(rule, gene_expr)
% evalGPR — evaluate a COBRA GPR rule with continuous expression values.
%   rule       string like 'x(12) | (x(7) & x(9))' — model.rules{r} format.
%              AND (&) evaluates as min, OR (|) as max (E-Flux convention).
%   gene_expr  numeric vector indexed same as model.genes.
%   v          scalar. Returns 1 for an empty/whitespace rule (reactions
%              with no GPR are unconstrained by expression).

    rule = strtrim(rule);
    if isempty(rule)
        v = 1;
        return
    end
    [v, ~] = parse_expr(rule, gene_expr);
end

function [v, rest] = parse_expr(s, ge)
%  OR has the lowest precedence.
    [v, s] = parse_term(s, ge);
    s = strtrim(s);
    while ~isempty(s) && s(1) == '|'
        s = strtrim(s(2:end));
        [v2, s] = parse_term(s, ge);
        v = max(v, v2);
        s = strtrim(s);
    end
    rest = s;
end

function [v, rest] = parse_term(s, ge)
%  AND binds tighter than OR.
    [v, s] = parse_factor(s, ge);
    s = strtrim(s);
    while ~isempty(s) && s(1) == '&'
        s = strtrim(s(2:end));
        [v2, s] = parse_factor(s, ge);
        v = min(v, v2);
        s = strtrim(s);
    end
    rest = s;
end

function [v, rest] = parse_factor(s, ge)
    s = strtrim(s);
    if isempty(s)
        v = 1; rest = ''; return
    end
    if s(1) == '('
        [v, s] = parse_expr(s(2:end), ge);
        s = strtrim(s);
        if ~isempty(s) && s(1) == ')'
            s = s(2:end);
        end
        rest = s;
        return
    end
    tok = regexp(s, '^x\((\d+)\)', 'tokens', 'once');
    if isempty(tok)
        % Unknown token — skip one char and treat as identity so a malformed
        % rule doesn't halt the whole sweep. Matches Stage 2 tolerance.
        v = 1;
        rest = s(2:end);
        return
    end
    idx = str2double(tok{1});
    if idx >= 1 && idx <= numel(ge)
        v = ge(idx);
    else
        v = 1;
    end
    rest = regexprep(s, '^x\(\d+\)', '', 'once');
end
