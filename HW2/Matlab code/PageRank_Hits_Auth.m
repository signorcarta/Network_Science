close all
clear all
clc

G = importdata('sport_edges.csv', ',',1);
T = readtable('sport_ids.csv');
T{:,'id'}=T{:,'id'}+1;

% adjacency matrix
G.data(:,:) = G.data(:,:)+1;
N = max(max(G.data));
A = sparse(G.data(:,2),G.data(:,1),ones(size(G.data,1),1),N,N);
A = 1*(A+A'>0); % undirected network
N = length(A); %# nodes

% GrU = graph(A,'upper');
% GrL = graph(A,'lower');
% Gr = digraph(A);
d = full(sum(A,2)); %degree of each node
%% Preprocessing __________________________________________________________
%remove self loops ________________________________________________________
for i = 1:N 
    A(i,i) = 0;
end

% Remove nodes which are NOT connected ____________________________________
pos = find(sum(A)~=0);
A = A(pos,pos);

% Remove dead ends (until none avalable) __________________________________
exit = false;
while (~exit)
    pos = find(sum(A)~=0);
    A = A(pos,pos);
    N = size(A,1);
    exit = isempty(find(sum(A)==0, 1));
end

% Find the largest connected component ____________________________________
e1 = [1;zeros(N-1,1)];
exit = false;
while(~exit)
    e1_old = e1;
    e1 = 1*(A*e1>0);
    exit = (sum(e1-e1_old)==0);
end
pos = find(e1);
A = A(pos,pos);
N = size(A,1);

%% Applying Page-Rank _____________________________________________________
c = 0.85; % damping factor
q = ones(N,1)/N; % normalized teleportation vector
M = A*sparse(diag(1./sum(A))); % normalized M
 
disp('Linear system solution. Page Rank Computing')
r = sparse((eye(N)-c*M)/(1-c))\q;
r = r/sum(r);
rank = table(r);


%% Applying HITS __________________________________________________________
epsilon = 1e-8;
M = A*A';
[pp,ee] = eigs(M,2);
p = -pp(:,1)/norm(pp(:,1));

N = size(M,1);
p0 = ones(N,1)/sqrt(N);
s = [];
t_conv = 1000;
for k = 1:30
    p00 = p0;
    p0 = A*(A'*p0);
    p0 = p0/norm(p0);
    s(k) = norm(p0-p00)/sqrt(N);
    if(norm(p0-p00)<epsilon && t_conv > 30) % computing # of iterations to get
                                            % convergence
        t_conv = k;
    end
end

% Table showing top-10 ranked according to Page-Rank ______________________
% rankT = table(r);
% T = [T rankT];% add the simRank score as the fourth column
% T = sortrows(T,'r','descend'); % sort the dataset 
% T(1:10,:)

% Table showing top-10 ranked according to Hits____________________________
% rankT = table(p);
% T = [T rankT];% add the simRank score as the fourth column
% T = sortrows(T,'p','descend'); % sort the dataset 
% T(1:10,:)

%% Comparison _____________________________________________________________ 

% figure(1)
% plot([p/sum(p),-r/sum(r)])
% grid
% legend('HITS','PageRank')
% title('PageRank vs HITS(Authorities) ')
% 
% figure(2)
% loglog(p/sum(p),r/sum(r),'x')
% grid
% xlabel('HITS score')
% ylabel('PageRank score')
% title('PageRank vs HITS(Authorities)')