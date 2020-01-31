close all
clear 
clc

G = importdata('CA-GrQc.txt', '\t', 4);

% adjacency matrix
N = max(max(G.data));
A = sparse(G.data(:,2),G.data(:,1),ones(size(G.data,1),1),N,N);

clear G;

%% Distance distribution and diameter _____________________________________
X = graph(A);
dist = distances(X);
dist = dist(dist<inf);
uniques = unique(dist); %show possible distances among nodes
pd = histc(dist, uniques);
pd = pd/sum(pd);
dist_vector = dist(:); %matrix to vector
diameter = max(dist_vector); %max distances between any two pairs

disp(['diameter: ---> ' num2str(diameter) ' <---'])

% Show result _____________________________________________________________
% figure(2);
% bar(pd)
% xlim([0,17])
% grid
% xlabel('distance')
% ylabel('PDF')
% title('Distance distribution')