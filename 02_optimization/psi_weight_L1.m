function [X, info] = psi_weight_L1(V, gamma, options)
%PSI_WEIGHT_L1 Proximal operator matched to weighted_L1.m.
% [X,info] = psi_weight_L1(V,gamma,options) approximately minimizes
%   0.5*norm(X-V,'fro')^2 + gamma*Phi(X), where
%   S = std(reshape(X,Nz,Nx,[]),0,3);
%   Phi(X) = sum(abs(S(:))) + sum(abs(X(:)))/size(X,2).

if nargin < 3 || isempty(options), options = struct(); end
opt = parse_options(options);
validateattributes(V, {'double','single'}, {'2d','nonempty','finite','nonsparse'}, mfilename, 'V');
validateattributes(gamma, {'numeric'}, {'scalar','real','finite','nonnegative'}, mfilename, 'gamma');
if size(V,1) ~= opt.Nz*opt.Nx
    error('psi_weight_L1:DimensionMismatch', ...
        'V must have Nz*Nx=%d rows; got %d.',opt.Nz*opt.Nx,size(V,1));
end
timer = tic;
V = double(V);
gamma = double(gamma);
T = size(V,2);
scale = max_abs_chunked(V,opt.ChunkPixels);
info = struct('converged',false,'iterations',0,'gamma',gamma, ...
    'Nz',opt.Nz,'Nx',opt.Nx,'frames',T,'inputScale',scale, ...
    'primalObjective',NaN,'dualObjective',NaN,'gap',NaN, ...
    'relativeGap',NaN,'solutionErrorBound',NaN, ...
    'elapsedSeconds',NaN,'status','not_started', ...
    'method','dual coordinate minimization for entrywise temporal-std L1');

if gamma == 0 || scale == 0
    X = V;
    info = exact_info(info,0,'exact_identity_or_zero',timer);
    return;
end

% Phi is positively homogeneous: prox_gamma(s*V) = s*prox_(gamma/s)(V).
% Internal scaling improves numerical range without changing the model.
V = V/scale;
alpha = (gamma/scale)/T;
if alpha >= 1
    % Entrywise L1 alone is sufficient to make zero optimal; adding the
    % nonnegative temporal penalty, which is zero at zero, preserves it.
    X = zeros(size(V),'like',V);
    objective = (0.5*squared_norm_chunked(V,opt.ChunkPixels)*scale)*scale;
    info = exact_info(info,objective,'exact_zero_solution',timer);
    return;
end
if T > 1
    beta = (gamma/scale)/sqrt(T-1);
else
    beta = 0; % MATLAB std along a singleton time dimension is zero.
end
absoluteTolerance = (opt.AbsGapTol/scale)/scale;
Q = zeros(size(V),'like',V);
last = struct();

for iteration = 1:opt.MaxIter
    % Dual feasible sets:
    %   |P_pt| <= alpha;
    %   mean(Q,2)=0 and norm(Q_p,:,2) <= beta for every spatial pixel p.
    % Coordinate projections minimize 0.5*norm(V-P-Q,'fro')^2.
    P = project_entry_disks(V-Q,alpha,opt.ChunkPixels);
    Q = project_temporal_dual(V-P,beta,opt);

    if iteration == 1 || mod(iteration,opt.CheckEvery)==0 || iteration==opt.MaxIter
        X = V-P-Q;
        last = gap_statistics(X,P,Q,alpha,beta,opt);
        if opt.Verbose
            fprintf('psi_weight_L1: iter=%d, relative gap=%.3e\n',iteration,last.relativeGap);
        end
        if last.gap <= absoluteTolerance + opt.RelGapTol*last.objective
            info.converged = true;
            break;
        end
    end
end

X = X*scale;
info.iterations = iteration;
info.primalObjective = (last.objective*scale)*scale;
info.gap = (last.gap*scale)*scale;
info.dualObjective = info.primalObjective-info.gap;
info.relativeGap = last.relativeGap;
info.solutionErrorBound = sqrt(2*last.gap)*scale;
info.elapsedSeconds = toc(timer);
if info.converged
    info.status = 'gap_tolerance_met';
else
    info.status = 'maximum_iterations_reached';
    message = sprintf(['Prox did not meet the requested primal-dual gap tolerance. ' ...
        'Iterations=%d, relative gap=%.3e, absolute gap=%.3e. ' ...
        'Increase MaxIter or investigate scaling; do not silently accept this in TwIST.'], ...
        iteration,info.relativeGap,info.gap);
    if opt.FailOnNonconvergence
        error('psi_weight_L1:NoConvergence','%s',message);
    else
        warning('psi_weight_L1:NoConvergence','%s',message);
    end
end
end

function opt = parse_options(user)
if ~isstruct(user) || ~isscalar(user)
    error('psi_weight_L1:InvalidOptions','options must be a scalar struct.');
end
opt = struct('Nz',285,'Nx',128,'MaxIter',1000,'CheckEvery',5, ...
    'RelGapTol',1e-8,'AbsGapTol',1e-12,'ProjectionRelTol',1e-13, ...
    'ProjectionMaxIter',100,'ChunkPixels',1024, ...
    'FailOnNonconvergence',true,'Verbose',false);
names = fieldnames(user);
for k=1:numel(names)
    if ~isfield(opt,names{k})
        error('psi_weight_L1:UnknownOption','Unknown option: %s.',names{k});
    end
    opt.(names{k}) = user.(names{k});
end
integers = {'Nz','Nx','MaxIter','CheckEvery','ProjectionMaxIter','ChunkPixels'};
for k=1:numel(integers)
    validateattributes(opt.(integers{k}),{'numeric'}, ...
        {'scalar','real','finite','integer','positive'},mfilename,integers{k});
end
validateattributes(opt.RelGapTol,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(opt.AbsGapTol,{'numeric'},{'scalar','real','finite','nonnegative'});
validateattributes(opt.ProjectionRelTol,{'numeric'},{'scalar','real','finite','positive','<',1});
validateattributes(opt.FailOnNonconvergence,{'logical'},{'scalar'});
validateattributes(opt.Verbose,{'logical'},{'scalar'});
end

function P = project_entry_disks(P,radius,blockSize)
for first=1:blockSize:size(P,1)
    ids=first:min(first+blockSize-1,size(P,1));
    block=P(ids,:);
    magnitude=abs(block);
    active=magnitude>radius;
    block(active)=block(active).*(radius./magnitude(active));
    P(ids,:)=block;
end
end

function Q = project_temporal_dual(B,radius,opt)
% Projection onto independent mean-zero temporal L2 balls for all pixels.
Q=zeros(size(B),'like',B);
if radius == 0, return; end
for first=1:opt.ChunkPixels:size(B,1)
    ids=first:min(first+opt.ChunkPixels-1,size(B,1));
    C=B(ids,:)-mean(B(ids,:),2);
    rowNorm=sqrt(sum(abs(C).^2,2));
    factor=ones(numel(ids),1);
    active=rowNorm>radius;
    factor(active)=radius./rowNorm(active);
    Q(ids,:)=C.*factor;
end
end

function stat = gap_statistics(X,P,Q,alpha,beta,opt)
% V-X=P+Q by construction. Fenchel gaps give the primal-dual gap
% without subtracting two potentially large, almost equal objectives.
amplitude=0; innerP=0; innerQ=0; fidelity=0;
r=zeros(size(X,1),1);
for first=1:opt.ChunkPixels:size(X,1)
    ids=first:min(first+opt.ChunkPixels-1,size(X,1));
    xb=X(ids,:); pb=P(ids,:); qb=Q(ids,:);
    centered=xb-mean(xb,2);
    r(ids)=sqrt(sum(abs(centered).^2,2));
    amplitude=amplitude+sum(abs(xb(:)));
    innerP=innerP+real(sum(conj(pb(:)).*xb(:)));
    innerQ=innerQ+real(sum(conj(qb(:)).*centered(:)));
    residual=pb+qb;
    fidelity=fidelity+0.5*sum(abs(residual(:)).^2);
end
f=alpha*amplitude;
g=beta*sum(r);
rawGap=(f-innerP)+(g-innerQ);
roundoff=128*eps*max([abs(f),abs(g),abs(innerP),abs(innerQ),realmin])*max(1,log2(numel(X)+1));
if rawGap < -roundoff
    error('psi_weight_L1:InvalidDualGap','Negative dual gap beyond floating-point tolerance.');
end
stat.gap=max(rawGap,0);
stat.objective=fidelity+f+g;
stat.relativeGap=stat.gap/max(stat.objective,realmin);
end

function value = max_abs_chunked(X,blockSize)
value=0;
for first=1:blockSize:size(X,1)
    block=X(first:min(first+blockSize-1,size(X,1)),:);
    value=max(value,max(abs(block(:))));
end
end

function value = squared_norm_chunked(X,blockSize)
value=0;
for first=1:blockSize:size(X,1)
    block=X(first:min(first+blockSize-1,size(X,1)),:);
    value=value+sum(abs(block(:)).^2);
end
end

function info = exact_info(info,objective,status,timer)
info.converged=true;
info.primalObjective=objective;
info.dualObjective=objective;
info.gap=0;
info.relativeGap=0;
info.solutionErrorBound=0;
info.status=status;
info.elapsedSeconds=toc(timer);
end
