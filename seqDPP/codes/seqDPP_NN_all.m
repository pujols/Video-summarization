function seqDPP_NN_all(dataset, rng_seed)

seg_size = 10;
Inf_type = 'exact';
C = inf;
%% Splitting the training and tsting videos
if (strcmp(dataset, 'OVP'))
    rng(rng_seed);
    inds_order = randperm(50);
    inds_tr = sort(inds_order(1 : 40));
    inds_te = sort(inds_order(41 : 50));
elseif (strcmp(dataset, 'YouTube'))
    rng(rng_seed);
    valid_inds = [11:21, 23:50]; % Ignore the cartoon videos
    inds_order = randperm(39);
    inds_tr = valid_inds(sort(inds_order(1 : 31)));
    inds_te = valid_inds(sort(inds_order(32 : 39)));
end

[videos, int_W, int_V, int_alpha] = initialize_videos(seg_size, dataset);

% training
display('Start training!')
[W, V, alpha, ~, ~, ~] = traindpp_NN_MLE(videos(inds_tr), C, int_W, int_V, int_alpha);
% summarization
display('Start generating summaries!')
videos_te = testdpp_inside(videos(inds_te), W, V, alpha, [], Inf_type);

% VSUMM evaluation
display('Start evaluation!')
approach_name = ['NN_' dataset];
res = seqDPP_evaluate(videos_te, inds_te, 1, approach_name, dataset);
display('% F-score, Recall, Precision: ');
disp(res);
end


%%------------------ functions ------------------
function [new_videos, int_W, int_V, int_alpha] = initialize_videos(seg_size, dataset)
disp('Initializing videos');
new_videos = struct;

if (strcmp(dataset, 'OVP'))
    load ../video_summarization/Oracle_groundset/Oracle_OVP.mat Oracle_record
    load ../data/Init/OVP_NN.mat W V alpha_
    
    for t = 1:size(Oracle_record, 1)
        new_videos(t).ids = 1:length(Oracle_record{t, 3});
        [new_videos(t).Ys, new_videos(t).grounds] = genSeg(length(Oracle_record{t, 3}), map_ids(Oracle_record{t, 3}, Oracle_record{t, 2}), seg_size);
        load(['../video_summarization/feature/OVP_v' num2str(Oracle_record{t, 1}) '/saliency.mat'], 'saliency_score');
        load(['../video_summarization/feature/OVP_v' num2str(Oracle_record{t, 1}) '/context2.mat'], 'context_score');
        load(['../video_summarization/feature/OVP_v' num2str(Oracle_record{t, 1}) '/fishers_PCA90.mat'], 'fishers');
        new_videos(t).fts = [saliency_score, context_score, fishers];
        if(length(Oracle_record{t, 3}) ~= size(new_videos(t).fts, 1))
            display(['wrong_video: ' num2str(t)]);
        end
        clear saliency_score context_score fishers
    end
elseif (strcmp(dataset, 'YouTube'))
    load ../video_summarization/Oracle_groundset/Oracle_Youtube.mat Oracle_record
    load ../data/Init/YouTube_NN.mat W V alpha_
    
    for t = 1:size(Oracle_record, 1)
        new_videos(t).ids = 1:length(Oracle_record{t, 3});
        [new_videos(t).Ys, new_videos(t).grounds] = genSeg(length(Oracle_record{t, 3}), map_ids(Oracle_record{t, 3}, Oracle_record{t, 2}), seg_size);
        load(['../video_summarization/feature/Youtube_v' num2str(Oracle_record{t, 1}) '/saliency.mat'], 'saliency_score');
        load(['../video_summarization/feature/Youtube_v' num2str(Oracle_record{t, 1}) '/context2.mat'], 'context_score');
        load(['../video_summarization/feature/Youtube_v' num2str(Oracle_record{t, 1}) '/fishers_PCA90.mat'], 'fishers');
        new_videos(t).fts = [saliency_score, context_score, fishers];
        if(length(Oracle_record{t, 3}) ~= size(new_videos(t).fts, 1))
            display(['wrong_video: ' num2str(t)]);
        end
        clear saliency_score context_score fishers
    end
end
clear Oracle_record;

for t = 1:length(new_videos)
    new_videos(t).YpredSeq = [];
    new_videos(t).Ypred = [];
end

% parameters
int_alpha = alpha_;
int_W = W;
int_V = V;
end


function [videos, Ls] = testdpp_inside(videos, W, V, alpha, Ls, Inf_type)

if ~exist('Ls', 'var')
    Ls = cell(1, length(videos));
elseif (isempty(Ls))
    Ls = cell(1, length(videos));
end
parfor i = 1 : length(videos) %parfor
    %fprintf('testdpp %d\n', i)
    [videos(i).YpredSeq, videos(i).Ypred, Ls{i}] = predictY_inside(videos(i), W, V, Ls{i}, alpha, Inf_type);
end
end


function [Y, Y_record, L] = predictY_inside(video, W, V, L, alpha, Inf_type)
% for "each" video
if isempty(L)
    X = video.fts;
    sX = sigmoid(W * X');
    YY = 2 * sX - 1;
    ZZ = V * YY;
    L = ZZ' * ZZ + alpha * eye(size(X, 1));
end

% correct L from numerical errors
L = (L + L') / 2;
% [V, Lam] = eig(full(L));
% Lam(Lam <0) = 0;
% L = V * Lam * V';    % project L into the PSD cone
% L = (L + L') / 2;

Gs = cellfun(@(Y) map_ids(video.ids, Y), video.grounds, 'UniformOutput', false);
Gs = cat(2, {[]}, Gs); % add a dummy video segment V0
Y = cell(1, length(Gs)); % For recording predicted results
Y_record = []; % For recording predicted results

for t = 2 : length(Gs)
    V = [Y{t-1}(:); Gs{t}(:)];     % this is Y_t-1 U V_{t}
    [~, Y_loc] = intersect(V, Y{t-1});
    [~, Gs_loc] = intersect(V, Gs{t});
    Y_loc = Y_loc(:); Gs_loc = Gs_loc(:);
    L_window = L(V, V);            % this is L^{t,t-1}
    Iv = diag([zeros(length(Y{t-1}), 1); ones(length(Gs{t}), 1)]);
    
    if (strcmp(Inf_type, 'exact'))
        whichs = {};
        for k = 1 : length(min(Gs{t}, 10))
            C = num2cell(nchoosek((1:length(Gs{t})), k), 2);
            whichs = [whichs; C];
        end

        % Starting with empty set
        best_comb = 0;
        best_J = det(L_window(Y_loc, Y_loc)) / det(L_window + Iv); 

        % Brute-force for each subset of Gs{t}(:)
        for at_comb = 1:length(whichs)
            J = det(L_window([Y_loc; Gs_loc(whichs{at_comb})], [Y_loc; Gs_loc(whichs{at_comb})])) / det(L_window + Iv);
            if (J > best_J)
                best_J = J;
                best_comb = at_comb;
            end
        end

        if (best_comb ~= 0)
            Y{t} = V(Gs_loc(whichs{best_comb}));
            Y_record = [Y_record, V(Gs_loc(whichs{best_comb}))'];
        end
        
    else        
        L_cond = inv(L_window + Iv);
        L_cond = inv(L_cond(Gs_loc, Gs_loc)) - eye(length(Gs_loc));
        select_loc_Gs = greedy_sym(L_cond);
        Y{t} = V(Gs_loc(select_loc_Gs));
        Y_record = [Y_record, V(Gs_loc(select_loc_Gs))'];
    end
end
Y_record = sort(unique(Y_record));
end


%% ================== train MLE =================
function [W, V, alpha, fval, exitflag, output] = traindpp_NN_MLE(videos, C, W0, V0, alpha_0)
addpath ./minFunc/
addpath ./minFunc/autoDif/
addpath ./minFunc/minFunc/compiled/
addpath ./minFunc/minFunc/mex/
addpath ./minFunc/minFunc/

% INPUT:
%   videos(k).fts: matrix, #frames-by-dim
%   videos(k).grounds: cell, 1-by-T
%   videos(k).Ys: cell, 1-by-T

cX = cell(size(videos));    % features
cG = cX;                    % Ground sets
cY = cX;                    % Labeled subsets
for c = 1 : length(cX)
    cX{c} = videos(c).fts;
    cG{c} = cellfun(@(Y)map_ids(videos(c).ids, Y), videos(c).grounds, 'UniformOutput', false);
    cY{c} = cellfun(@(Y)map_ids(videos(c).ids, Y), videos(c).Ys, 'UniformOutput', false); 
end

n = size(videos(1).fts, 2);
m = numel(W0) / n;
p = numel(V0) / m;
theta_reg = {zeros(m, n), zeros(p, m)};

% minimize the hinge loss
options.Display = 'off';
options.Method = 'lbfgs';
options.optTol = 1e-10;
options.progTol = 1e-10;
options.MaxIter = 500;
options.MaxFunEvals = 500;
% options.DerivativeCheck = 'on';

funObj = @(arg)compute_fg(arg, cX, cG, cY, C, theta_reg, n, m, p);
[theta, fval, exitflag, output] = minFunc(funObj, [W0(:); V0(:); alpha_0], options);

% recover W, V, and alpha
W = theta(1:m*n);           W = reshape(W, m, n);
V = theta(m*n+1:m*n+p*m);   V = reshape(V, p, m);
alpha = max(theta(end), 1e-6);
end


function [f, g] = compute_fg(theta, cX, cG, cY, C, theta_reg, n, m, p)

W = theta(1:m*n);           W = reshape(W, m, n);
V = theta(m*n+1:m*n+p*m);   V = reshape(V, p, m);
alpha = max(theta(end), 1e-6);
M = V' * V;

f = 0;
gW = zeros(size(W));
gV = zeros(size(V));
gA = 0;
parfor k = 1 : length(cX) %parfor
    %fprintf('process video %d\n',k)
    [fk, gWk, gVk, gAk] = compute_fg_one_data(W, V, alpha, M, cX{k}, cG{k}, cY{k}, n, m, p);
    f = f - fk;
    gW = gW - gWk;
    gV = gV - gVk;
    gA = gA - gAk;
end

% g = [gW(:); gV(:); gA]; 
g = [gW(:); zeros(size(gV(:))); gA];
if ~isinf(C)
    diff_W = (W(:) - theta_reg{1}(:));
    diff_V = (V(:) - theta_reg{2}(:));
    f = C * f + 0.5 * (diff_W' * diff_W) + 0.5 * (diff_V' * diff_V);
    g = C * g + [diff_W; diff_V, 0];
end
end


function [f, gW, gV, galpha] = compute_fg_one_data(W, V, alpha, M, X, Gs, Ys, n, m, p)

% input: 
%   theta: the parameter to learn
%   X: #frames-by-dim, features
%   Gs: ground sets, cell, 1-by-T
%   Ys: labeled subsets, cell, 1-by-T
%   n: # features
%   m: W is m-by-n
%   p: V is p-by-m

Gs = cat(2, {[]}, Gs); % add a dummy video segment V0
Ys = cat(2, {[]}, Ys); % add a dummy subset Y0  

sX = sigmoid(W * X');    % m-by-#data
YY = 2 * sX - 1;
ZZ = V * YY;             % p-by-#data
S = 2 * sX .* (1-sX);

f = 0;
gW = zeros(m, n);
gV = zeros(p, m);
galpha = 0;
for t = 2 : length(Gs)
    Y = [Ys{t-1}(:); Ys{t}(:)];     % this is Y_t-1 U Y_{t}
    G = [Ys{t-1}(:); Gs{t}(:)];     % this is Y_t-1 U V_{t}
    [~,Y] = intersect(G, Y);
    L = ZZ(:, G)' * ZZ(:, G) + alpha * eye(length(G));      % this is L^{t,t-1}
    
    % some extra work to fix the numerical issues
    L = (L + L') / 2;
%     [U, Gamma] = svd(L);
%     Gamma(Gamma < 0) = 0;
%     L = U * Gamma * U';
%     L = (L + L') / 2;
    
    Iv = diag([zeros(length(Ys{t-1}), 1); ones(length(Gs{t}), 1)]);
    
    % compute function value
    J = log(det(L(Y,Y))) - log(det(L + Iv));
    
    % compute gradients
    [gWt, gVt, gAt] = compute_g(X(G,:), YY(:,G), W, V, S(:,G), M, L, Y, Iv);
    
    % overall
    f = f + J;
    gW = gW + gWt;
    gV = gV + gVt;
    galpha = galpha + gAt;
end
end


function [gW, gV, gA] = compute_g(X, YY, W, V, S, M, L, Y, Iv)
% gradient

% partial J / parigal L
Ainv = inv(L + Iv);
if isempty(Y)
    LYinv = 0;
else
    LYinv = zeros(size(L));
    LYinv(Y, Y) = inv(L(Y,Y));
end
J = LYinv - Ainv;

% (partial J / partial V)
gV = 2 * V * (YY * J * YY');

% (partial J / partial W)
gW = zeros(size(W));
for i = 1 : size(W, 1)
    D = diag(S(i,:));
    gW(i,:) = 2 * M(i,:) * YY * J * D * X;
end

% (partial J / partial alpha)
gA = trace(J);
end


%% ================== small functions =======================
function Y_mapped = map_ids(ids, Y)
[~, Y_mapped] = intersect(ids(:), Y(:));
if length(Y_mapped) ~= length(Y)
    error('error in map_ids(): Y is out of ids');
end
end
  

function sigZ = sigmoid(Z)
sigZ = 1./(1+exp(-Z));
end


function [Ys, grounds] = genSeg(nFrm, gt, seg_size)
nSeg = ceil(nFrm / seg_size);
Ys = cell(1, nSeg);
grounds = cell(1, nSeg);
for n = 1:nSeg
    grounds{n} = (n-1)*seg_size+1 : min(n*seg_size, nFrm);
    Ys{n} = gt((gt >= grounds{n}(1)) & (gt <= grounds{n}(end)));
end
end


function S = greedy_sym(L)
% Initialize.
N = size(L, 1);
S = [];

% Iterate through items.
for i = 1:N
    lift0 = obj(L, [S, i]) - obj(L, S);
    lift1 = obj(L, [S, (i+1):N]) - obj(L, [S, i:N]);
    
    if lift0 > lift1
        S = [S, i];
    end
end
end