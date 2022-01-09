%% WIOD, 2014, 48 industries, 2 countries
% constructs the baseline dataset
%%  
%clear everything
clear
close all
clc

%control parameters
N = 48; %no industries
M = 2; %no countries

%% OBS - DATA OBSERVED and ADJUSTED


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
        filename = sprintf('WIOD_%s_2014_48.xls',M_countries(i));
k_obs(:,i) = readmatrix(filename,'Sheet','SEA','Range','I2:I49');    
cap_obs(:,i) = readmatrix(filename,'Sheet','SEA','Range','H2:H49');    
e_obs(:,i) = readmatrix(filename,'Sheet','SEA','Range','F2:F49');    
lab_obs(:,i) = readmatrix(filename,'Sheet','SEA','Range','G2:G49');    
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
        filename = sprintf('WIOD_%s_2014_48.xls',M_countries(i));
pd_obs(:,:,i) = readmatrix(filename,'Sheet','II','Range','B2:AW49');       
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
        filename = sprintf('WIOD_%s_2014_48.xls',M_countries(i));
pf_obs(:,:,i) = readmatrix(filename,'Sheet','M','Range','B2:AW49');       
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

% aggregate GDP
PY_obs = sum(py_obs);



% PARAMETERS 

% beta: final consumption shares on Y (sum to one)
beta = py_obs ./ PY_obs;
%sum(beta,1)

% mu: output shares on Y (sum greater than one)
mu = pq_obs ./ PY_obs;
%sum(mu,1)


%Gamma matrix, where gamma(j,i) is usage by industry i of
%industry j's good (row is origin, column is destination)
% Gamma = pd_obs./pq_obs';  % ERROR: Array dimensions must match for binary array op.
pq_obs_T = reshape(pq_obs,[N,1,M]);
pq_obs_T = (permute(pq_obs_T, [2 1 3]));
Gamma = pd_obs./pq_obs_T;

%Gammai = sum(Gamma,1) is total intermediate share of sector i
gammai = sum(Gamma,1);
gammai = (permute(gammai, [2 1 3])); 
gammai = reshape(gammai,[N,M]);


%Sigma matrix, where sigma(j,i) is usage by industry i of
%industry j's good (row is origin, column is destination)
Sigma = pf_obs./pq_obs_T;

%Sigmai = sum(Sigma,1) is total import share of sector i
sigmai = sum(Sigma,1);
sigmai = (permute(sigmai, [2 1 3])); 
sigmai = reshape(sigmai,[N,M]);



% test mu = [I - Gamma]^(-1)*beta
aux = (eye(N) - Gamma);
inverse = zeros(N,N,M);
for t = 1:M
auxname = aux(:,:,t);
inverse(:,:,t)=inv(auxname);
end
%size(inverse)
%size(beta)

mu_test = zeros (N,M);
for t = 1:M
mu_test(:,t)=inverse(:,:,t)*beta(:,t);
end

if abs(sum(mu - mu_test)) >= 1e-15
    disp('check parameters')
end




% aggregate factors shares 
alpha = wK_obs.*K_obs./VA_obs;

if max(abs((1-alpha)-wE_obs.*E_obs./VA_obs)) >= 1e-15
    disp('check agg factors shares')
end


% industry-specific factor shares 
alphai = (wK_obs.*k_obs)./(pq_obs.*(1-gammai-sigmai)); 
if max(abs((1-alphai)-wE_obs.*e_obs./VAi_obs)) >= 1e-13
    disp('check ind factors shares')
end


% industry-specific factor shares 
if max(abs(pq_obs.*(1-gammai-sigmai) - VAi_obs)) >= 1e-9
    disp('check ind factors shares')
end

% condition (24) 20200206_Chapter3
DD = sum(mu.*(1-gammai)); %sum of Domar shares discounting intermediate shares
if abs(DD - 1) >= 1e-15
    disp('check parameters')
end



%logs: log(X) returns the natural logarithm ln(x) of each element in array X
lbeta = log(beta);
lmu = log(mu);
lGamma = log(Gamma);
lgammai = log(gammai);
lSigma = log(Sigma);
lsigmai = log(sigmai);

lalpha = log(alpha);
lalphai = log(alphai);

lki_obs = log(k_obs);
lei_obs = log(e_obs);

lK_obs = log(K_obs);
lE_obs = log(E_obs);



% PRICES

% domestic prices
p_data = zeros (N,M);
for i=1:length(M_countries) 
        filename = sprintf('WIOD_%s_2014_48.xls',M_countries(i));
p_data(:,i) = readmatrix(filename,'Sheet','SEA','Range','K2:K49');     
end

p_data(isnan(p_data))=1;


% real values and agg prices that would be generated with data prices
y_data = py_obs./p_data;
Y_data_agg = prod(y_data.^beta);
P_data = PY_obs./Y_data_agg;

%normalising prices to P=1
p_obs = p_data./P_data;
y_obs = py_obs./p_obs;
Y_obs = prod(y_obs.^beta); % REAL OUTPUT: Cobb-Douglas aggregator final good

P_obs = PY_obs./Y_obs;
if abs(P_obs - 1) >= 1e-14          
    disp('P not equal to 1')
end


% test real = nominal
if abs(Y_obs - PY_obs) >= 1e-07     
    disp('Cobb-Douglas dnt hold')
end


% test FOC final good
if max(abs(py_obs - beta.*PY_obs)) >= 1e-9
    disp('FOC final good dnt hold')
end


% logs
lp_obs = log(p_obs);


% test condition used to solve model
if abs(sum(beta.*lbeta)-sum(beta.*lp_obs)) >= 1e-14 % error = 3.9968e-15
    disp('Beta-Price condition dnt hold')
end


% other REAL VARIABLES
q_obs = pq_obs./p_obs;

d_obs = zeros (N,N,M);
for t = 1:M
d_obs(:,:,t)=pd_obs(:,:,t)./p_obs(:,t);
end


% world prices
wp_obs = readmatrix('WIOD_WP_IRL_KOR.xls','Sheet','WP','Range','D2:D49');

% other REAL VARIABLES

f_obs = zeros (N,N,M);
for t = 1:M
f_obs(:,:,t)=pf_obs(:,:,t)./wp_obs;
end


% logs
lq_obs = log(q_obs);
ld_obs = log(d_obs);
lf_obs = log(f_obs);
ly_obs = log(y_obs);

lwp_obs = log(wp_obs);






% OTHER TESTS %

alpha_formula = sum(alphai.*(1-gammai-sigmai).*mu)./(1-sum(sigmai.*mu));
if abs(alpha_formula-alpha) >= 1e-15
    disp('alpha not matching formula')
end



% backing out productivity
% sum(Gamma.*ld_obs,1): sum over rows, same column / (1X3)
lA = zeros (N,M);
for t = 1:M
lA(:,t) = lq_obs(:,t) - (1-gammai(:,t)-sigmai(:,t)).*(alphai(:,t).*lki_obs(:,t) + (1-alphai(:,t)).*lei_obs(:,t)) - sum(Gamma(:,:,t).*ld_obs(:,:,t),1)' - sum(Sigma(:,:,t).*lf_obs(:,:,t),1)';
end
A = exp(lA);



%% Analysis

industry = (1:N)';

% aggregate labour share
E_share = (1-alpha);

% industry labour shares
e_shares = (1-alphai);



%% Save data
% ALL VARIABLES ARE IN VECTORS (48 X 2)%

clear alpha_formula DD mu_test aux auxname inverse filename i m t

save('WIOD_IRL_KOR')



%% Export calculated data


% Gamma per country
for i=1:length(M_countries) 
        CC = table(Gamma(:,:,i));
        filepath = 'C:\Users\Veri\OneDrive - University of Essex\RESEARCH\PUBLISHING\Chapter3\Codes\Stata\Input'; 
        filename = sprintf('Gamma_%s_48.xls',M_countries(i));
        writetable(CC,fullfile(filepath,filename),'Sheet',1,'Range','A1','WriteVariableNames',true)
end
for i=1:length(M_countries) 
        CC = table(Gamma(:,:,i));
        filepath = 'C:\Users\Veri\OneDrive - University of Essex\RESEARCH\PUBLISHING\Chapter3\Codes\Gephi'; 
        filename = sprintf('Gamma_%s_48.csv',M_countries(i));
        writetable(CC,fullfile(filepath,filename),'WriteVariableNames',false)
end



% Sigma per country
for i=1:length(M_countries) 
        CC = table(Sigma(:,:,i));
        filepath = 'C:\Users\Veri\OneDrive - University of Essex\RESEARCH\PUBLISHING\Chapter3\Codes\Stata\Input'; 
        filename = sprintf('Sigma_%s_48.xls',M_countries(i));
        writetable(CC,fullfile(filepath,filename),'Sheet',1,'Range','A1','WriteVariableNames',true)
end
for i=1:length(M_countries) 
        CC = table(Sigma(:,:,i));
        filepath = 'C:\Users\Veri\OneDrive - University of Essex\RESEARCH\PUBLISHING\Chapter3\Codes\Gephi'; 
        filename = sprintf('Sigma_%s_48.csv',M_countries(i));
        writetable(CC,fullfile(filepath,filename),'WriteVariableNames',false)
end


% industrial data per country
for m=1:length(M_countries) 
        IND = table(industry,e_shares(:,m),mu(:,m),gammai(:,m),sigmai(:,m));
        filepath = 'C:\Users\Veri\OneDrive - University of Essex\RESEARCH\PUBLISHING\Chapter3\Codes\Stata\Input'; 
        filename = sprintf('Industries_%s_48.xls',M_countries(m));
        writetable(IND,fullfile(filepath,filename),'Sheet',1,'Range','A1','WriteVariableNames',true)
end


% agg data both countries
ALL = table(M_countries,E_share');
    filepath = 'C:\Users\Veri\OneDrive - University of Essex\RESEARCH\PUBLISHING\Chapter3\Codes\Stata\Input'; 
    filename = 'WIOD_2014_IRL_KOR.xls';
    writetable(ALL,fullfile(filepath,filename),'Sheet',1,'Range','A1','WriteVariableNames',true)
 
    



%% PLOTS main body
 


% Figure 3.5  

 % IO plot
 for i=1:length(M_countries)
imagesc(Gamma(:,:,i));
set(findall(gcf,'-property','FontSize'),'FontSize',12);
title(sprintf('Domestic IO Matrix - Share of input purchases in total cost (%s 2014)',M_countries(i)),'FontSize',14);
         x0=50, y0=50, width=900, height=600;
         set(gcf,'position',[x0,y0,width,height]);
xlabel('Buying industries'), ylabel('Selling industries');
colorbar;
colormap (flipud(pink));
axis xy;
xticks([1:48]), xtickangle(90);
yticks([1:48]);
filepath = 'C:\Users\Veri\OneDrive - University of Essex\RESEARCH\PUBLISHING\Chapter3\Codes\Exhibits';
filename = sprintf('%s_IO_48.png',M_countries(i));
saveas(gcf,fullfile(filepath,filename),'jpeg');
 end

 
 % II plot
 for i=1:length(M_countries)
imagesc(Sigma(:,:,i));
set(findall(gcf,'-property','FontSize'),'FontSize',12)
title(sprintf('Imported II Matrix - Share of imports purchases in total cost (%s 2014)',M_countries(i)),'FontSize',14);
         x0=50, y0=50, width=900, height=600;
         set(gcf,'position',[x0,y0,width,height])
xlabel('Buying industries'), ylabel('Selling industries');
colorbar;
colormap (flipud(pink));
axis xy;
xticks([1:48]), xtickangle(90)
yticks([1:48])
filepath = 'C:\Users\Veri\OneDrive - University of Essex\RESEARCH\PUBLISHING\Chapter3\Codes\Exhibits';
filename = sprintf('%s_M_48.png',M_countries(i));
saveas(gcf,fullfile(filepath,filename),'jpeg');
 end

 
%% Plots Appendix


% Figure 3.B.9 
y = [p_data';wp_obs'];
bar(industry,y','EdgeColor','none','BarWidth',1)
set(findall(gcf,'-property','FontSize'),'FontSize',12)
title('Domestic and International Prices','FontSize',14), xlabel('Industries'), ylabel('Index')
         x0=50, y0=50, width=900, height=600;
         set(gcf,'position',[x0,y0,width,height])
xticks([1:48]), xtickangle(90)
legend('IRL','KOR','World','Location','northwest')
%createfigure1(industry, y')
filepath = 'C:\Users\Veri\OneDrive - University of Essex\RESEARCH\PUBLISHING\Chapter3\Codes\Exhibits';
filename = sprintf('prices_data_IRL_KOR.png');
saveas(gcf,fullfile(filepath,filename));
