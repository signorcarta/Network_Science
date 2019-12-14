close all
clear 
clc

G = importdata('CA-GrQc.txt', '\t', 4);

% adjacency matrix
N = max(max(G.data));
A = sparse(G.data(:,2),G.data(:,1),ones(size(G.data,1),1),N,N);

clear G;

%% Clustering coefficient _________________________________________________
cluster_coeffs = zeros(1, N);
d = full(sum(A,1)); % computing degrees 
d = d(d>0); % discard zero degrees
for i=find(d>1)
    neigh_indexes = find(A(i, :)); %indexed of nodes connected to i
    L = A(neigh_indexes, neigh_indexes);
    neigh_connections_i = sum(L(:)); % # of links between neighs
    max_possible_neighs_connections_i = (d(i) * (d(i)-1))/2;
    cluster_coeffs(i) = neigh_connections_i / max_possible_neighs_connections_i;
end

%average clustering coefficient of degrees
average_C_coeff = mean(cluster_coeffs);
disp([ 'average clustering coefficient: ---> ' num2str(average_C_coeff) ' <---']);

% Show result _____________________________________________________________
% figure(3);
% subplot(2,1,1)
% loglog(cluster_coeffs, 'x')
% grid
% xlabel('d')
% ylabel('Cluster coefficients')
% title('LogLog plot')
% subplot(2,1,2)
% semilogy(cluster_coeffs, 'x')
% grid
% xlabel('d')
% ylabel('Cluster coefficients')
% title('Semilog plot')
