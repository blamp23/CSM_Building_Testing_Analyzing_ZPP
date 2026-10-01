function ensure_solver()
%ENSURE_SOLVER  Set COBRA LP solver globals on a worker without triggering
%initCobraToolbox -- which fires off git submodule updates that collide
%across parallel workers.
%
%Uses `persistent` so only the first parfor iter on each worker actually
%touches the globals; subsequent iters return immediately.
%
%Pattern lifted verbatim from
%  Tanguay_Data/GIMME/Final_Pipeline/7_Gene_Robustness_Analysis/ensure_solver.m
%where it is the known-good approach for this codebase.
persistent initialised
if isempty(initialised) || ~initialised
    global CBTLPSOLVER CBT_LP_SOLVER CBT_MILP_SOLVER
    CBTLPSOLVER     = 'gurobi';
    CBT_LP_SOLVER   = 'gurobi';
    CBT_MILP_SOLVER = 'gurobi';
    initialised = true;
end
end
