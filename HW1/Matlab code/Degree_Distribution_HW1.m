close all
clear all
clc

G = importdata('CA-GrQc.txt', '\t', 4);
N = max(max(G.data));
A = sparse(G.data(:,2), G.data(:,1), ones(size(G.data,1),1), N, N);

clear G;

%% Degree distribution ____________________________________________________
% distribution
d = full(sum(A,1)); % computing degrees 
d = d(d>0); % discard zero degrees
k = unique(d); % degree samples
pk = histc(d,k)'; % counts occurrencies
pk = pk/sum(pk); % normalize to 1

% K_min, K_max, average
k_min = min(d);
k_max = max(d);
disp(['(min K, max K): ---> (' num2str(k_min) ', ' num2str(k_max) ') <---'])
avg_degree = mean(d);
disp(['average degree: ---> ' num2str(avg_degree) ' <---'])

% cumulative distribution
Pk = cumsum(pk,'reverse');

% log binning
klog = 10.^(0:0.1:ceil(log10(max(k)))); % identifies the intervals
pklog = histc(d,klog)'; % counts occurrencies
pklog = pklog/sum(pklog); % normalize to 1

%% Moments ________________________________________________________________
variance = var(d);
disp(['variance: ---> ' num2str(variance) ' <---'])
skewness = skewness(d);
disp(['skewness: ---> ' num2str(skewness) ' <---'])
kurtosis = kurtosis(d);
disp(['kurtosis: ---> ' num2str(kurtosis) ' <---'])

%% ML Estimation of Gamma _________________________________________________
kmin = 8; %% empirically derived from the graph

d2 = d(d>=kmin); % range restriction
ga = 1+1/mean(log(d2/kmin)); % exponent estimation
disp(['gamma ML estimated: ---> ' num2str(ga) ' <---'])


c = (ga-1)*kmin^(ga-1); % constant c
C = sum(pk(k>=kmin)); % constant C
                      % correction factor taking into account for  
                      % the fractional impact of the degrees k>=kmin

% fitting in the PDF signal
s1 = C*c*k.^(-ga); 
% fitting in the CCDF signal
s3 = C*c/(ga-1)*k.^(1-ga);
% fitting in the log version of the PDF signal
s2 = C*c/(ga-1)*klog.^(1-ga) *(1-(klog(2)/klog(1))^(1-ga));

% Show Results ____________________________________________________________
% figure(1)
% subplot(2,2,1)
% plot(k,pk,'.')
% grid
% xlabel('k')
% ylabel('PDF')
% title('linear PDF')
% subplot(2,2,2)
% loglog(k,pk,'.')
% hold on
% loglog(k,s1);
% hold off
% axis([xlim min(pk/2) 2*max(pk)])
% grid
% xlabel('k')
% ylabel('PDF')
% title('logarithmic PDF')
% subplot(2,2,3)
% loglog(klog,pklog,'.')
% hold on
% loglog(klog,s2);
% hold off
% axis([xlim min(pklog/2) 2*max(pklog)])
% grid
% xlabel('k')
% ylabel('PDF')
% title('logarithmic PDF (log bins)')
% subplot(2,2,4)
% loglog(k,Pk,'.')
% hold on
% loglog(k,s3);
% hold off
% axis([xlim min(Pk/2) 2])
% grid
% xlabel('k')
% ylabel('CCDF')
% title('logarithmic CCDF')