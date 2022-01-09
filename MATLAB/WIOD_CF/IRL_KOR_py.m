%% WIOD, 2014, 48 industries, 2 countries
% assist in selecting industries to be merged based on negative py 
%%  
%clear everything
clear
close all
clc

%control parameters
N = 56; %no industries
M = 2; %no countries

%% OBS - DATA OBSERVED and ADJUSTED

%%% WIOD files for IRL and KOR have to be moved manually 
%%% from folder WIOD to folder WIOD_CF
%%% WIOD_IRL_2014.xls and WIOD_KOR_2014.xls

% Factors quantities and payment to factos

%capital
k_obs = zeros (N,M);
cap_obs = zeros (N,M);
%total labour employed
e_obs = zeros (N,M);
lab_obs = zeros (N,M);

M_countries =[
"IRL"
"KOR"
];

for i=1:length(M_countries) 
    filepath = 'C:\Users\Veri\OneDrive - University of Essex\RESEARCH\PUBLISHING\Chapter3\Codes\MATLAB\WIOD'; 
    filename = sprintf('WIOD_%s_2014.xls',M_countries(i));
k_obs(:,i) = readmatrix(fullfile(filepath, filename),'Sheet','SEA','Range','I2:I57');    
cap_obs(:,i) = readmatrix(fullfile(filepath, filename),'Sheet','SEA','Range','H2:H57');    
e_obs(:,i) = readmatrix(fullfile(filepath, filename),'Sheet','SEA','Range','F2:F57');    
lab_obs(:,i) = readmatrix(fullfile(filepath, filename),'Sheet','SEA','Range','G2:G57');    
end

% replace zeros for very small numbers
k_obs0 = k_obs;
k_obs0(k_obs0==0) = 1/1000000000;
k_obs = k_obs0;

cap_obs0 = cap_obs;
cap_obs0(cap_obs0==0) = 1/1000000000;
cap_obs = cap_obs0;

e_obs0 = e_obs;
e_obs0(e_obs0==0) = 1/1000000000;
e_obs = e_obs0;

lab_obs0 = lab_obs;
lab_obs0(lab_obs0==0) = 1/1000000000;
lab_obs = lab_obs0;

clear k_obs0 cap_obs0 e_obs0 lab_obs0

% AGGREGATES
K_obs = sum(k_obs);
PAY_K = sum(cap_obs);
E_obs = sum(e_obs);
PAY_E = sum(lab_obs);

% factor prices (average payment to factors) 
wK_obs = PAY_K./K_obs;
wE_obs = PAY_E./E_obs;

% payments to factors in each industry
pay_ki = wK_obs.*k_obs;
pay_ei = wE_obs.*e_obs;

% value added in each industry
VAi_obs = pay_ki + pay_ei;

% adjusted aggregate VA (millions of US$, current prices)
VA_obs = sum(VAi_obs);




% intermediate use - DOMESTIC
%d matrix, where d(j,i) is usage by industry i of
%industry j's good (row is origin, column is destination)

pd_obs = zeros (N,N,M);

for i=1:length(M_countries) 
    filepath = 'C:\Users\Veri\OneDrive - University of Essex\RESEARCH\PUBLISHING\Chapter3\Codes\MATLAB\WIOD'; 
    filename = sprintf('WIOD_%s_2014.xls',M_countries(i));
pd_obs(:,:,i) = readmatrix(fullfile(filepath, filename),'Sheet','II','Range','B2:BE57');       
end


% replace zeros for very small numbers
pd_obs0 = pd_obs;
pd_obs0(pd_obs0==0) = 1/1000000000000;
pd_obs = pd_obs0;
clear pd_obs0


ii_cons = sum(pd_obs,1); % Error using ' Transpose on ND array is not defined. Use PERMUTE instead.
ii_cons = (permute(ii_cons, [2 1 3])); 
ii_cons = reshape(ii_cons,[N,M]); 

ii_prod = sum(pd_obs,2);
ii_prod = reshape(ii_prod,[N,M]);

% intermediate use - FOREIGN
%f matrix, where f(j,i) is usage by industry i of
%industry j's good (row is origin, column is destination)

pf_obs = zeros (N,N,M);

for i=1:length(M_countries) 
    filepath = 'C:\Users\Veri\OneDrive - University of Essex\RESEARCH\PUBLISHING\Chapter3\Codes\MATLAB\WIOD'; 
    filename = sprintf('WIOD_%s_2014.xls',M_countries(i));
pf_obs(:,:,i) = readmatrix(fullfile(filepath, filename),'Sheet','M','Range','B2:BE57');       
end


% replace zeros for very small numbers e-12
pf_obs0 = pf_obs;
pf_obs0(pf_obs0==0) = 1/1000000000000;
pf_obs = pf_obs0;
clear pf_obs0


m_cons = sum(pf_obs,1); % Error using ' Transpose on ND array is not defined. Use PERMUTE instead.
m_cons = (permute(m_cons, [2 1 3])); % but it's 29x1x7 whereas 29x7 would do
m_cons = reshape(m_cons,[N,M]); % but then it doesn't do the calculations afterwards


% output (millions of US$, current prices)
pq_obs = VAi_obs + ii_cons + m_cons;

% final consumption (millions of US$, current prices)
py_obs = pq_obs - ii_prod;



%% Export py

industry = (1:N)';
py = table(industry,py_obs); % columns: industry, IRL, KOR
    filename = sprintf('IRL_KOR_py.xls');
    writetable(py,filename,'Sheet',1,'Range','A1','WriteVariableNames',true)
        