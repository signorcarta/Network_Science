close all
clear all
clc

G = importdata('sport_edges.csv', ',',1);
T = readtable('sport_ids.csv');
T{:,'id'}=T{:,'id'}+1; % all indexes must be >=1
T = sortrows(T,'new_id','ascend');
T{5258:end,'new_id'}=T{5258:end,'new_id'}-1; % fixing indexes 5258:end

% adjacency matrix
G.data(:,:) = G.data(:,:)+1;
N = max(max(G.data));
A = sparse(G.data(:,2),G.data(:,1),ones(size(G.data,1),1),N,N);
N = length(A); % # of nodes

%% Preprocessing __________________________________________________________
A = 1*(A+A'>0); % undirected network

% Remove self loops _______________________________________________________
for i=1:N
    A(i,i)=0;
end

% Remove nodes which are NOT connected ____________________________________

pos = find(sum(A)~=0);
A = A(pos,pos);

% Remove dead ends  _______________________________________________________
exit = false;

 while (~exit)
    pos = find(sum(A)~=0);
    A = A(pos,pos);
    exit = isempty(find(sum(A)==0, 1));
end

N = size(A,1); % update N

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
N = size(A,1); % update N

% Some useful things ______________________________________________________
Au = A;
d = full(sum(Au)); % degrees vector
D = sum(d); % degrees sum
I = spdiags(ones(N,1),0,N,N); % identity matrix
Di = spdiags(1./sqrt(d'),0,N,N); % diagonal degrees square-rooted
L = I - Di*Au*Di; % normalized Laplacian
M = Au*Di*Di; % normalized adjacancy matrix

%% Spectral Clustering technique __________________________________________

% extract eigenvectors
[V,DD] = eigs(L,2,'SA');
Vv = Di*V; % normalize eigenvectors
v1 = Vv(:,2)/norm(Vv(:,2)); % Fiedler's vector
% sweep wrt the ordering identified by v1
[v1s,pos] = sort(v1,'descend'); % reorder the adjacency matrix
Au1 = Au(pos,pos);

% evaluate the conductance measure
a = sum(triu(Au1));
b = sum(tril(Au1));
assoc = cumsum(a+b);
assoc = min(assoc,D-assoc);
cut = cumsum(b-a);
conduct = cut./assoc;
conduct = conduct(1:end-1);
% identify the minimum to be used as threshold
[~,mpos] = min(conduct);
threshold = mean(v1s(mpos:mpos+1));

disp('Spectral Clustering technique')
disp(['---> Minimum conductance: ' num2str(conduct(mpos))])
disp(['---> Community #1 size: ' num2str(mpos)])
disp(['---> Community #2 size: ' num2str(N-mpos)])
disp([' '])
disp(['Displaying the smallest community:'])
T(1:mpos,:)

G = graph(A);
p = plot(G);
highlight(p,1:mpos,'NodeColor','green')

%% PageRank-nibble technique ______________________________________________

if mpos<N-mpos  % select seed node from the smaller group
    i = pos(1); % select the more relevant 
else
    i = pos(end);
end

q = zeros(N,1);
q(i) = 1; % teleport vector
c = 0.85;
r = (I-c*M)\((1-c)*q); % ranking vector
ep = 1e-3; % precision

% run PageRank-nibble
u = zeros(N,1); % starting point
v = q; % starting point
th = full(ep*d/D)'; % thresholds
count = 0; % exit counter
complexity = 0; % complexity value (# of operations)
ii = i; % starting index used for Push operation

while (count<N)
    if v(ii)>th(ii) % push if above threshold
        tmp = v(ii);
        u(ii) = u(ii)+(1-c)*tmp;
        v(ii) = 0;
        v = v + c*M(:,ii)*tmp;    
        complexity = complexity + d(ii); % update complexity
        count = 0; % reset the exit counter
    else % go to next entry if below threshold
        count = count + 1; % increase the exit counter
        ii = mod(ii,N)+1; % update the index used for Push
    end
end

% sweep wrt the ordering identified by v1
% reorder the adjacency matrix
[u1s,pos2] = sort(u,'descend');
Nmax = find(u1s>0,1,'last'); % discard nodes with 0 values (never used in Push)
Au1 = Au(pos2,pos2(1:Nmax));

% evaluate the conductance measure
a = sum(triu(Au1));
b = sum(tril(Au1));
assoc = cumsum(a+b);
assoc = min(assoc,D-assoc);
cut = cumsum(b-a);
conduct = cut./assoc;
conduct = conduct(1:Nmax-1); 
% identify the minimum to be used as threshold
[~,mpos2] = min(conduct);
threshold2 = mean(u1s(mpos2:mpos2+1));

disp('PageRank-nibble technique')
disp(['---> Minimum conductance: ' num2str(conduct(mpos2))])
disp(['---> Community size #1: ' num2str(mpos2)])
disp(['---> Community size #2: ' num2str(N-mpos2)])
disp([' '])
disp(['Displaying the smallest community:'])
T(1:mpos2,:)

G = graph(A);
p = plot(G);
highlight(p,1:mpos,'NodeColor','red')

% show network with partition
% figure(1)
% plot(u,v1,'k.')
% hold on
% plot(threshold2*[1,1],ylim,'g-')
% plot(xlim,threshold*[1,1],'r-')
% hold off
% grid
% ylabel('Fiedler''s eigenvector value')
% xlabel('PageRank value')
% title('Communities')