%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NETWORK ANALYSIS
% adapted from:
% Vasco M Carvalho
% From Micro to Macro via Production Networks
% Journal of Economic Perspectives, 2014
% Section 3 Computations and Figures 3 and 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOAD DATA 

%clear everything
clear
close all
clc

% 2014 World Input Output Data (WIOD)
load WIOD.mat

%% CREATE UNWEIGHTED VERSION OF ADJACENT MATRIX
Gamma(isnan(Gamma))=0; % coefficients

%Discretize IO into binary matrix (Create Unweighted Directed Graph)
%Pick threshold rule: >1% of total intermediate input bill
t=0.01;

for m=1:min(size(Gamma))
    for i=1:max(size(Gamma))
        for j=1:max(size(Gamma))
            if Gamma(i,j,m)>t
                io_bin(i,j,m)=1; 
            else
                io_bin(i,j,m)=0;
            end
        end
    end
end


%% NETWORK DENSITY

density=zeros(1,M);
for i=1:length(M_countries) 
    bin_matrix=io_bin(:,:,i);
    density(:,i)=nnz(bin_matrix)/(max(size(bin_matrix))^2);
end

%% DISTANCES AND DIAMETERS

% uses function "charpath": see description there
% uses function "distance_bin": see description there
for i=1:length(M_countries) 
    bin_matrix=io_bin(:,:,i);
    distance_matrix_bin=distance_bin(bin_matrix);
    [lambda(:,i),efficiency(:,i),ecc(:,i),radius(:,i),diameter(:,i)] = charpath(distance_matrix_bin);
end

distance=lambda;
clear lambda;




%% WEIGHTED DISTANCE

% uses function "distance_wei_floyd": see description there
% uses function "charpath": see description there
SPL = zeros(N,N,M); % aka distance matrix
hops = zeros(N,N,M); 
Pmat = zeros(N,N,M);
for m=1:M
    [SPL(:,:,m),hops(:,:,m),Pmat(:,:,m)] = distance_wei_floyd(Gamma(:,:,m),'log');
    [distance_w(:,m),efficiency_w(:,m),ecc_w(:,m),radius_w(:,m),diameter_w(:,m)] = charpath(SPL(:,:,m));
end



%% WEIGHTED INDEGREE 

idw = gammai;

%% WEIGHTED OUTDEGREE 

% gammao = sum(Gamma,2) for a country represents total intermediate share
% of an industry's buyers (downstream measure)
gammao = sum(Gamma,2);
gammao = (permute(gammao, [2 1 3])); 
gammao = reshape(gammao,[56,43]);

odw = gammao;


%% DOWNSTREAM CENTRALITY COMPUTATIONS

centrality_data = mu;

GammaX=Gamma;
GammaX(GammaX==0)=0.0000000001; % to avoid singularity


% % Katz-Bonacich % % 

% ``Because the length of walks in a graph is unbounded, Katz-Bonacich centrality requires
% a discount factor -- a factor ? between 0 and 1 -- to compute the discounted sum of walks
% emanating from the node.''
% c^KB = (I - Gamma)^(-1) Gamma*1
% Gamma*1 = sum of rows -> outdegree -> first term of the infinite sum


% no discounting and no rescaling, gammas are very small, series converge
cent_KB_down = zeros(N,M); 
for m=1:M
    cent_KB_down(:,m) = (eye(56)-GammaX(:,:,m))\sum(GammaX(:,:,m),2); 
end



% Carvalho
% Transforming share matrix so that all IO columns sum to 1 (following Acemoglu et
% al 2012, ECMA)
% CLARIFICATION: columns in the IO data are rows in the paper's W
% Acemoglu nao da ponto sem no -> matrix is transpose and indegrees are all equal 1
X=GammaX;
XX=zeros(N,N,M);
for m=1:M
    for i=1:N
        if sum(X(:,i,m))>0
            XX(:,i,m)=X(:,i,m)./sum(X(:,i,m));
        elseif sum(X(:,i,m))==0
            XX(:,i,m)=X(:,i,m);
        end
    end
end
W=XX;
clear X XX;


% Computing Centrality Score as Carvalho
% alpha was the share of intermediate inputs in production = 0.5
% n = number of sectors (56 industries)
% 
v = zeros(N,N,M); 
cent_Carv_down = zeros(N,M); 

for m=1:M
    v(:,:,m)=(0.5/56)*inv(eye(56)-0.5*W(:,:,m));
    cent_Carv_down(:,m)=sum(v(:,:,m)');
end



% eigenvector centrality
V = zeros(N,N,M); 
D = zeros(N,N,M); 
Ones = ones(N,M);
EigValues = zeros(N,M);
EigValues_max_down = zeros(1,M);
cent_eig_down= zeros(N,M); 
for m=1:M
    [V(:,:,m),D(:,:,m)] = eig(GammaX(:,:,m));
    EigValues(:,m)=D(:,:,m)*Ones(:,m);
    EigValues_max_down(1,m)=max(abs(EigValues(:,m)));
    EigVal_i=find(max(abs(EigValues(:,m))));
    cent_eig_down(:,m) = V(:,EigVal_i,m);
end

X=cent_eig_down;
XX= X < 0; % it's important that values are all neg or pos within vector
for m=1:M
        if X(:,m)<=0
            X(:,m)=X(:,m).*(-1);
        else
            X(:,m)=X(:,m);
        end
end
cent_eig_down=X;
clear V D Ones EigValues X XX;

EigValues_max_down=EigValues_max_down';


% for curiosity, check relationship Carvalho's alpha = lambda = 0.5
gammai0=gammai;
gammai0(isnan(gammai0))=0;
gammai0_avg = zeros(1,M);
for m=1:M
       gammai0_avg(:,m) = transpose(sum(gammai0(:,m), 1) / N);
end

gammai0_avg=gammai0_avg';
curious=EigValues_max_down-gammai0_avg % nothing systematic, but very close for US 



%% UPSTREAM CENTRALITY COMPUTATIONS

% transpose Gamma
for m=1:M
    GammaXT(:,:,m)=GammaX(:,:,m)';
end


% % Katz-Bonacich % % 
% no discounting and no rescaling, gammas are very small, series converge
cent_KB_up = zeros(N,M); 
for m=1:M
    cent_KB_up(:,m) = (eye(56)-GammaXT(:,:,m))\sum(GammaXT(:,:,m),2); 
end


% eigenvector centrality
V = zeros(N,N,M); 
D = zeros(N,N,M); 
Ones = ones(N,M);
EigValues = zeros(N,M);
EigValues_max_up = zeros(1,M);
cent_eig_up= zeros(N,M); 
for m=1:M
    [V(:,:,m),D(:,:,m)] = eig(GammaXT(:,:,m));
    EigValues(:,m)=D(:,:,m)*Ones(:,m);
    EigValues_max_up(1,m)=max(abs(EigValues(:,m)));
    EigVal_i=find(max(abs(EigValues(:,m))));
    cent_eig_up(:,m) = V(:,EigVal_i,m);
end

X=cent_eig_up;
XX= X < 0; % it's important that values are all neg or pos within vector
for m=1:M
        if X(:,m)<=0
            X(:,m)=X(:,m).*(-1);
        else
            X(:,m)=X(:,m);
        end
end
cent_eig_up=X;
clear V D Ones EigValues X XX;

EigValues_max_up=EigValues_max_up';


% Carvalho
X=zeros(N,N,M); 
% transpose each country's Gamma
for i=1:length(M_countries)
    X(:,:,i)=Gamma(:,:,i)';
end

XX=zeros(N,N,M);
for m=1:M
    for i=1:N
        if sum(X(:,i,m))>0
            XX(:,i,m)=X(:,i,m)./sum(X(:,i,m));
        elseif sum(X(:,i,m))==0
            XX(:,i,m)=X(:,i,m);
        end
    end
end

W1=XX;
clear X XX;

v1 = zeros(N,N,M); 
cent_Carv_up = zeros(N,M); 

for m=1:M
    v1(:,:,m)=(0.5/56)*inv(eye(56)-0.5*W1(:,:,m));
    cent_Carv_up(:,m)=sum(v1(:,:,m)');
end


%% CENTRALITY OF IMPORTED GOODS

Sigma(isnan(Sigma))=0; % coefficients


% using KB downstram
cent_Sigma_down = zeros(N,M); 
for m=1:M
    cent_Sigma_down(:,m)=Sigma(:,:,m)*cent_KB_down(:,m);
end



% using KB upstram
cent_Sigma_up = zeros(N,M); 
for m=1:M
    cent_Sigma_up(:,m)=Sigma(:,:,m)*cent_KB_up(:,m);
end




%% Export calculated data

% data per country
density=density';
diameter=diameter';
diameter_w=diameter_w';
distance=distance';
distance_w=distance_w';


     NET = table(M_countries,density,distance,distance_w,diameter,diameter_w);
        filepath = 'C:\Users\Veri\OneDrive - University of Essex\RESEARCH\PUBLISHING\Chapter3\Codes\Stata\Input'; 
        filename = 'WIOD_Network.xls';
        writetable(NET,fullfile(filepath,filename),'Sheet',1,'Range','A1','WriteVariableNames',true)


%% Extra        
% industry data per country
% industry = (1:N)';
% 
% for i=1:length(M_countries) 
%         IND = table(industry,ecc(:,i),gammao(:,i), cent_Carv_down(:,i), cent_Carv_up(:,i), cent_eig_down(:,i), cent_eig_up(:,i), cent_KB_down(:,i), cent_KB_up(:,i), cent_Sigma_down(:,i),cent_Sigma_up(:,i));
%         filepath = 'C:\Users\Veri\OneDrive - University of Essex\RESEARCH\PUBLISHING\Chapter3\Codes\Stata\Input'; 
%         filename = sprintf('Industries_Network_%s_2014.xls',M_countries(i));
%         writetable(IND,fullfile(filepath,filename),'Sheet',1,'Range','A1','WriteVariableNames',true)
% end
% 
% 
