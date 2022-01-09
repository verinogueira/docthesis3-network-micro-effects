%% WIOD, 2014, 56 industries, 43 countries
% constructs the dataset for the analyis of the WIOD's sample of countries
%%  
%clear everything
clear
close all
clc

%control parameters
N = 56; %no industries
M = 43; %no countries

%% OBS - DATA OBSERVED and ADJUSTED


% Factors quantities and payment to factos

% gross output
go_data = zeros (N,M);

%capital
k_data = zeros (N,M);
cap_data = zeros (N,M);
%total labour employed
e_data = zeros (N,M);
lab_data = zeros (N,M);

M_countries =[
"AUS"
"AUT"
"BEL"
"BGR"
"BRA"
"CAN"
"CHE"
"CHN"
"CYP"
"CZE"
"DEU"
"DNK"
"ESP"
"EST"
"FIN"
"FRA"
"GBR"
"GRC"
"HRV"
"HUN"
"IDN"
"IND"
"IRL"
"ITA"
"JPN"
"KOR"
"LTU"
"LUX"
"LVA"
"MEX"
"MLT"
"NLD"
"NOR"
"POL"
"PRT"
"ROU"
"RUS"
"SVK"
"SVN"
"SWE"
"TUR"
"TWN"
"USA"
];

for i=1:length(M_countries) 
        filename = sprintf('WIOD_%s_2014.xls',M_countries(i));
go_data(:,i) = readmatrix(filename,'Sheet','SEA','Range','C2:C57');    
k_data(:,i) = readmatrix(filename,'Sheet','SEA','Range','I2:I57');    
cap_data(:,i) = readmatrix(filename,'Sheet','SEA','Range','H2:H57');    
e_data(:,i) = readmatrix(filename,'Sheet','SEA','Range','F2:F57');    
lab_data(:,i) = readmatrix(filename,'Sheet','SEA','Range','G2:G57');    
end


% AGGREGATES
K_data = sum(k_data)';
CAP_data = sum(cap_data)';
E_data = sum(e_data)';
LAB_data = sum(lab_data)';
GO_data = sum(go_data)';


% value added in each industry
va_data = cap_data + lab_data;

% adjusted aggregate VA (millions of US$, current prices)
VA_data = sum(va_data)';



% intermediate use - DOMESTIC
%d matrix, where d(j,i) is usage by industry i of
%industry j's good (row is origin, column is destination)

pd_data = zeros (N,N,M);

for i=1:length(M_countries) 
        filename = sprintf('WIOD_%s_2014.xls',M_countries(i));
pd_data(:,:,i) = readmatrix(filename,'Sheet','II','Range','B2:BE57');       
end

ii_cons = sum(pd_data,1); % Error using ' Transpose on ND array is not defined. Use PERMUTE instead.
ii_cons = (permute(ii_cons, [2 1 3])); % but it's 29x1x7 whereas 29x7 would do
ii_cons = reshape(ii_cons,[56,43]); % but then it doesn't do the calculations afterwards

II_cons = sum(ii_cons)';


% intermediate use - FOREIGN
%f matrix, where f(j,i) is usage by industry i of
%industry j's good (row is origin, column is destination)

pf_data = zeros (N,N,M);

for i=1:length(M_countries) 
        filename = sprintf('WIOD_%s_2014.xls',M_countries(i));
pf_data(:,:,i) = readmatrix(filename,'Sheet','M','Range','B2:BE57');       
end

m_cons = sum(pf_data,1); % Error using ' Transpose on ND array is not defined. Use PERMUTE instead.
m_cons = (permute(m_cons, [2 1 3])); % but it's 29x1x7 whereas 29x7 would do
m_cons = reshape(m_cons,[56,43]); % but then it doesn't do the calculations afterwards

M_cons = sum(m_cons)';

% total intermediate consumption
iim_cons = ii_cons + m_cons;
IIM_cons = sum(iim_cons)';

% ii SEA different sum WIOT
error_data = go_data - va_data - ii_cons - m_cons;

% production internediate sold domestically
ii_prod = sum(pd_data,2);
ii_prod = reshape(ii_prod,[56,43]);

II_prod = sum(ii_prod)';


% !!!! I was wrong, we don't need that GDP = VA doesn't include trade
% total production intermediate inputs (domestic and exported)
iix_prod_data = zeros (N,M);

for i=1:length(M_countries) 
        filename = sprintf('WIOD_%s_2014.xls',M_countries(i));
iix_prod_data(:,i) = readmatrix(filename,'Sheet','II_PROD','Range','D2:D57'); 
end

IIX_prod = sum(iix_prod_data)';


% MODEL: trade balanced doesnt hold
% GDP value-added (income) - vertical
VA_model = GO_data - II_cons - M_cons;
% GDP final demand (expenditure): C + X - horizontal
PY_model = GO_data - II_prod;

% peanuts: vertical difference
error_VA = VA_model - VA_data;


% balanced trade assumption
X_obs = IIX_prod - II_prod; % underestimated, not counting C + I (horizontal)
TRADE_obs = X_obs - M_cons;

% industries' leftover: domestic consumption and exports
py_obs = go_data - ii_prod;

% aggregate final consumption
PY_obs = sum(py_obs)'; % identical to PY_model by construction




% PARAMETERS 

% beta: final consumption shares on Y (sum to one)
beta = py_obs ./ PY_obs';

% mu: output shares on Y (sum greater than one)
mu = go_data ./ PY_obs';


%Gamma matrix, where gamma(j,i) is usage by industry i of
%industry j's good (row is origin, column is destination)
% Gamma = pd_obs./pq_obs';  % ERROR: Array dimensions must match for binary array op.
go_data_M = reshape(go_data,[56,1,43]);
go_data_M = (permute(go_data_M, [2 1 3]));
Gamma = pd_data./go_data_M;


%Gammai = sum(Gamma,1) is total intermediate share of sector i
gammai = sum(Gamma,1);
gammai = (permute(gammai, [2 1 3])); 
gammai = reshape(gammai,[56,43]);


%Sigma matrix, where sigma(j,i) is usage by industry i of
%industry j's good (row is origin, column is destination)

Sigma = pf_data./go_data_M;

sigmai = sum(Sigma,1);
sigmai = (permute(sigmai, [2 1 3])); 
sigmai = reshape(sigmai,[56,43]);

Sigma(isnan(Sigma))=0;



% aggregate factors shares (UK data on compensation)
alpha = CAP_data./VA_data;


% industry-specific factor shares 
alphai = cap_data./va_data; 


%% Export calculated data

industry = (1:N)';
save('WIOD.mat')


% all countries agg data 
ALL = table(M_countries,alpha,PY_obs,VA_data,GO_data,II_cons, II_prod, IIM_cons, IIX_prod,M_cons, X_obs, TRADE_obs);
    filepath = 'C:\Users\Veri\OneDrive - University of Essex\RESEARCH\PUBLISHING\Chapter3\Codes\Stata\Input'; 
    filename = 'WIOD_2014_ALL.xls';
    writetable(ALL,fullfile(filepath,filename),'Sheet',1,'Range','A1','WriteVariableNames',true)
 
