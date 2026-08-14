
clear; close all; clc;

%% ---------------------------------------------------------------------
%  PART 1: Dynamical IPR vs n  
%  -----------------------------------------------------------------

nvals = 4:2:70;

bar_clique_th  = zeros(size(nvals));  bar_clique_num  = zeros(size(nvals));
bar_bridge_th  = zeros(size(nvals));  bar_bridge_num  = zeros(size(nvals));

v1_hub_th      = zeros(size(nvals));  v1_hub_num      = zeros(size(nvals));
v1_clique_th   = zeros(size(nvals));  v1_clique_num   = zeros(size(nvals));

v2_hub_th      = zeros(size(nvals));  v2_hub_num      = zeros(size(nvals));
v2_bridge_th   = zeros(size(nvals));  v2_bridge_num   = zeros(size(nvals));
v2_clique_th   = zeros(size(nvals));  v2_clique_num   = zeros(size(nvals));

for idx = 1:numel(nvals)
    n = nvals(idx);

    % ---------------- Barbell B_n ----------------
    A = build_barbell(n);
    Pi = limiting_distribution(normalize_adj(A));
    bar_clique_num(idx) = dynamical_ipr(Pi, 1);     % vertex A1 (clique interior)
    bar_bridge_num(idx) = dynamical_ipr(Pi, n);     % vertex Br_A (bridge)
    bar_clique_th(idx)  = (1 - 4/n + 1/n^2);
    bar_bridge_th(idx)  = (1/2 - 2/n + 1/n^2);

    % ---------------- Star of Cliques, Variant 1 (full connection) ----------------
    A1 = build_star_v1(n);
    Pi1 = limiting_distribution(normalize_adj(A1));
    v1_hub_num(idx)    = dynamical_ipr(Pi1, 1);     % hub |0>
    v1_clique_num(idx) = dynamical_ipr(Pi1, 2);     % vertex |1,1>
    v1_hub_th(idx)     = (n^4 + 2*n^2 + 5) / (n+1)^4;   % exact, Eq. (40)
    v1_clique_th(idx)  = (1 - 4/n+ 1/n^2);

    % ---------------- Star of Cliques, Variant 2 (single connection) ----------------
    A2 = build_star_v2(n);
    Pi2 = limiting_distribution(normalize_adj(A2));
    v2_hub_num(idx)    = dynamical_ipr(Pi2, 1);     % hub |0>
    v2_bridge_num(idx) = dynamical_ipr(Pi2, 2);     % vertex |1,1> = b_1 (bridge)
    v2_clique_num(idx) = dynamical_ipr(Pi2, 3);     % vertex |1,2> (clique interior)
    v2_hub_th(idx)     = 1/4;                              % exact, Eq. (55)
    v2_bridge_th(idx)  = 1-8/n+1/n^2;                  % exact, Table I
    v2_clique_th(idx)  = 1-4/n+1/n^2;                  % exact, Table I
end

figure('Name','Numerical versus analytical dynamical IPR','Color','w');

subplot(1,3,1);
plot(nvals, bar_clique_th, '-', 'LineWidth', 1.5, 'Color', [0 0 0]); hold on;
plot(nvals, bar_clique_num, 'o', 'MarkerSize', 4, 'Color', [0 0 0]);
plot(nvals, bar_bridge_th, '-', 'LineWidth', 1.5, 'Color', [0 0.45 0.74]);
plot(nvals, bar_bridge_num, 's', 'MarkerSize', 4, 'Color', [0 0.45 0.74]);
xlabel('n'); ylabel('Dynamical IPR'); title('(a) Barbell');
legend('Clique (th)','Clique (num)','Bridge (th)','Bridge (num)', ...
    'Location','southeast');
ylim([0 1]); grid on;

subplot(1,3,2);
plot(nvals, v1_hub_th, '-', 'LineWidth', 1.5, 'Color', [0.85 0 0]); hold on;
plot(nvals, v1_hub_num, 'o', 'MarkerSize', 4, 'Color', [0.85 0 0]);
plot(nvals, v1_clique_th, '-', 'LineWidth', 1.5, 'Color', [0 0 0]);
plot(nvals, v1_clique_num, 's', 'MarkerSize', 4, 'Color', [0 0 0]);
xlabel('n'); title('(b) Star of Cliques V1');
legend('Centre (th)','Centre (num)','Clique (th)','Clique (num)', ...
    'Location','southeast');
ylim([0 1]); grid on;

subplot(1,3,3);
plot(nvals, v2_hub_th, '-', 'LineWidth', 1.5, 'Color', [0 0 0]); hold on;
plot(nvals, v2_hub_num, 'o', 'MarkerSize', 4, 'Color', [0 0 0]);
plot(nvals, v2_bridge_th, '-', 'LineWidth', 1.5, 'Color', [0.85 0 0]);
plot(nvals, v2_bridge_num, 's', 'MarkerSize', 4, 'Color', [0.85 0 0]);
plot(nvals, v2_clique_th, '-', 'LineWidth', 1.5, 'Color', [0 0.45 0.74]);
plot(nvals, v2_clique_num, '^', 'MarkerSize', 4, 'Color', [0 0.45 0.74]);
xlabel('n'); title('(c) Star of Cliques V2');
legend('Centre (th)','Centre (num)','Bridge (th)','Bridge (num)', ...
    'Clique (th)','Clique (num)','Location','east');
ylim([-1 1]); grid on;

%% ---------------------------------------------------------------------
%  PART 2: Limiting distribution heatmaps pi_ij  ( n = 6)
%  -----------------------------------------------------------------

n_fig = 6;   % Barbell -> 12 vertices; Star graphs -> 37 vertices

figure('Name','Limiting Distribution pi_ij','Color','w');

subplot(1,3,1);
A = build_barbell(n_fig);
Pi = limiting_distribution(normalize_adj(A));
imagesc(Pi); axis square; colorbar; colormap(parula);
xlabel('Initial vertex j'); ylabel('Vertex i'); title('Barbell');

subplot(1,3,2);
A1 = build_star_v1(n_fig);
Pi1 = limiting_distribution(normalize_adj(A1));
imagesc(Pi1); axis square; colorbar; colormap(parula);
xlabel('Initial vertex j'); ylabel('Vertex i'); title('Star of Cliques V1');

subplot(1,3,3);
A2 = build_star_v2(n_fig);
Pi2 = limiting_distribution(normalize_adj(A2));
imagesc(Pi2); axis square; colorbar; colormap(parula);
xlabel('Initial vertex j'); ylabel('Vertex i'); title('Star of Cliques V2');

%% =======================================================================
%  Local functions
%  =======================================================================

function A = build_barbell(n)
% Barbell B_n: two K_n cliques joined by ONE bridge edge.
% Vertices 1..n = clique A, vertices n+1..2n = clique B.
% Bridge edge connects vertex n (Br_A) to vertex n+1 (Br_B).
    N = 2*n;
    A = zeros(N);
    A(1:n,1:n)         = ones(n) - eye(n);
    A(n+1:2*n,n+1:2*n) = ones(n) - eye(n);
    A(n, n+1) = 1; A(n+1, n) = 1;
end

function A = build_star_v1(n)
% Star of Cliques, Variant 1 (full connection):
% hub (vertex 1) connected to EVERY vertex of n identical K_n cliques.
% Total vertices: 1 + n^2.
    N = 1 + n*n;
    A = zeros(N);
    hub = 1;
    for j = 1:n
        idx = 1 + (j-1)*n + (1:n);
        A(idx, idx) = ones(n) - eye(n);   % clique j (complete K_n)
        A(hub, idx) = 1; A(idx, hub) = 1; % hub attached to all clique verts
    end
end

function A = build_star_v2(n)
% Star of Cliques, Variant 2 (single connection):
% hub (vertex 1) connected to only ONE vertex |j,1> per clique; that
% vertex is still part of clique j's complete K_n. Total: 1 + n^2.
    N = 1 + n*n;
    A = zeros(N);
    hub = 1;
    for j = 1:n
        idx = 1 + (j-1)*n + (1:n);
        A(idx, idx) = ones(n) - eye(n);   % clique j (complete K_n, all n verts)
        bridge = idx(1);                  % |j,1> = the clique's bridge vertex
        A(hub, bridge) = 1; A(bridge, hub) = 1;
    end
end

function M = normalize_adj(A)
% Normalized adjacency M~ = D^{-1/2} A D^{-1/2} 
    d = sum(A, 2);
    Dinv = diag(1 ./ sqrt(d));
    M = Dinv * A * Dinv;
    M = (M + M.') / 2;   % symmetrize away roundoff
end

function Pi = limiting_distribution(M)
% Exact Cesaro-limit distribution 
%   pi_ij = sum over DISTINCT eigenvalues E of | (P_E)_ij |^2
% where P_E is the orthogonal projector onto the eigenspace of E.
    N = size(M,1);
    [V, D] = eig(M);
    E = diag(D);

    tol = 1e-8 * max(1, max(abs(E)));
    [Es, order] = sort(E);
    Vs = V(:, order);

    Pi = zeros(N);
    m = 1;
    while m <= N
        grp = m;
        while grp < N && abs(Es(grp+1) - Es(m)) < tol
            grp = grp + 1;
        end
        Vg = Vs(:, m:grp);          % eigenvectors spanning this eigenspace
        P = Vg * Vg.';              % projector onto eigenspace
        Pi = Pi + P.^2;             % elementwise square, then accumulate
        m = grp + 1;
    end
end

function ipr = dynamical_ipr(Pi, j)
% Dynamical IPR for a walk started at vertex j: IPR_j = sum_i pi_ij^2
    ipr = sum(Pi(:, j).^2);
end

function y = clamp01(x)
    y = min(max(x, 0), 1);
end
