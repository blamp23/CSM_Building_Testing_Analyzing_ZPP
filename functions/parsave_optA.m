function parsave_optA(fp, S_active, active_rxns, meta)
% Helper because MATLAB's `save` can't be called from inside a parfor loop.
save(fp, 'S_active', 'active_rxns', 'meta', '-v7.3');
end
