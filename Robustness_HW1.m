close all
clear all
clc

G = importdata('CA-GrQc.txt', '\t', 4);
N = max(max(G.data));
A = sparse(G.data(:,2), G.data(:,1), ones(size(G.data,1),1), N, N);
clear G;

% Note: in this section some built-in Matlab function i.e. conccomp() and 
%      degree() will be used.

%% Robustness to random failures___________________________________________
X = graph(A);
labels = conncomp(X); %finding out the indexes of the connected components
giantcomp_size = max(hist(labels)); % computing giant connected component                                    
p_inf_0 = giantcomp_size/N;
p_inf = zeros(1,N); % vector of size N for the p_inf

for i=1:N
    delete_node = randi(numnodes(X)); %index of the node that will be deleted
    X = rmnode(X,delete_node);  %the new graph without the deleted node
    labels = conncomp(X);
    giantcomp_size = max(hist(labels));
    p_inf(i)=giantcomp_size/N;
end

%% Robustness to random attacks____________________________________________
X = graph(A);
labels = conncomp(X); %finding the indexes of the connected components
giantcomp_size = max(hist(labels));


p_inf_0_att = giantcomp_size/N;
p_inf_att = zeros(1,N);   

for t = 1:N
    degrees = degree(X);
    [~,deleted_node] = max(degrees); %now the target is the hub
    X = rmnode(X,deleted_node);
    label = conncomp(X);
    giantcomp_size = max(hist(label));
    p_inf_att(t) = giantcomp_size/N;
end

%% Show results ____________________________________________________________
f = 1:100:N;
p_inf = p_inf(f);
p_inf_att = p_inf_att(f);
figure(7);
plot(f,p_inf_att/p_inf_0_att,'X-');
hold on
plot(f,p_inf/p_inf_0,'O-');
legend('Attacks', 'Random failures', 'Location','northeast')
hold off
grid
xlabel('f');
ylabel('p_{inf}(f)/p_{inf}(0)');
title('Robustness');