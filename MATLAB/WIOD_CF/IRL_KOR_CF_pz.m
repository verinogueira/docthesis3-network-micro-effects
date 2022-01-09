%% WIOD, 2014, 48 industries, 2 countries
% counterfactual: import price shock
% % Industry #1 A-B % % 
%%  
%clear everything
clear
close all
clc


%% LOAD DATA 

%clear everything
clear
close all
clc

load('WIOD_IRL_KOR')


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solving the counterfactual equilibrium:

% change all notation to _cf 
% % tests below are for the case of no shock % %

%List of things potentially changed in the counterfactual:
K_cf = K_obs;
E_cf = E_obs;
wp_cf = wp_obs;
wp_cf(1,:) = wp_cf(1,:)*1.01; % % Industry #1 A-B % % 

lK_cf = log(K_cf);
lE_cf = log(E_cf);
lwp_cf = log(wp_cf);


%1) solve for factor prices and Y
lwK_cf = zeros (1,M);
for t = 1:M
lwK_cf(:,t) = (alpha(:,t)-1).*(lK_cf(:,t)-lalpha(:,t)) + (1-alpha(:,t)).*(lE_cf(:,t)-log(1-alpha(:,t))) + ...
    sum(mu(:,t).*(1-gammai(:,t)-sigmai(:,t)).*alphai(:,t).*lalphai(:,t))./(sum(mu(:,t).*(1-gammai(:,t)-sigmai(:,t)))) + ...
    sum(mu(:,t).*(1-gammai(:,t)-sigmai(:,t)).*(1-alphai(:,t)).*log(1-alphai(:,t)))./(sum(mu(:,t).*(1-gammai(:,t)-sigmai(:,t)))) + ... 
    sum(mu(:,t).*(1-gammai(:,t)-sigmai(:,t)).*log(1-gammai(:,t)-sigmai(:,t)))./(sum(mu(:,t).*(1-gammai(:,t)-sigmai(:,t)))) + ...
    sum(mu(:,t).*lA(:,t))./(sum(mu(:,t).*(1-gammai(:,t)-sigmai(:,t)))) + ...
    sum(mu(:,t).*sum(Gamma(:,:,t).*log(Gamma(:,:,t)),1)')./(sum(mu(:,t).*(1-gammai(:,t)-sigmai(:,t)))) + ...
    sum(mu(:,t).*sum(Sigma(:,:,t).*log(Sigma(:,:,t)),1)')./(sum(mu(:,t).*(1-gammai(:,t)-sigmai(:,t)))) + ...
    (sum(beta(:,t).*lbeta(:,t)) - sum(mu(:,t).*sum(Sigma(:,:,t).*lwp_cf(:,1),1)'))./(sum(mu(:,t).*(1-gammai(:,t)-sigmai(:,t))));
end
wK_cf = exp(lwK_cf);

lwE_cf = zeros (1,M);
for t = 1:M
lwE_cf(:,t) = alpha(:,t).*(lK_cf(:,t)-lalpha(:,t)) + (-alpha(:,t)).*(lE_cf(:,t)-log(1-alpha(:,t))) + ...
    sum(mu(:,t).*(1-gammai(:,t)-sigmai(:,t)).*alphai(:,t).*lalphai(:,t))./(sum(mu(:,t).*(1-gammai(:,t)-sigmai(:,t)))) + ...
    sum(mu(:,t).*(1-gammai(:,t)-sigmai(:,t)).*(1-alphai(:,t)).*log(1-alphai(:,t)))./(sum(mu(:,t).*(1-gammai(:,t)-sigmai(:,t)))) + ... 
    sum(mu(:,t).*(1-gammai(:,t)-sigmai(:,t)).*log(1-gammai(:,t)-sigmai(:,t)))./(sum(mu(:,t).*(1-gammai(:,t)-sigmai(:,t)))) + ...
    sum(mu(:,t).*lA(:,t))./(sum(mu(:,t).*(1-gammai(:,t)-sigmai(:,t)))) + ...
    sum(mu(:,t).*sum(Gamma(:,:,t).*log(Gamma(:,:,t)),1)')./(sum(mu(:,t).*(1-gammai(:,t)-sigmai(:,t)))) + ...
    sum(mu(:,t).*sum(Sigma(:,:,t).*log(Sigma(:,:,t)),1)')./(sum(mu(:,t).*(1-gammai(:,t)-sigmai(:,t)))) + ...
    (sum(beta(:,t).*lbeta(:,t)) - sum(mu(:,t).*sum(Sigma(:,:,t).*lwp_cf(:,1),1)'))./(sum(mu(:,t).*(1-gammai(:,t)-sigmai(:,t))));
end
wE_cf = exp(lwE_cf);

lY_cf = zeros (1,M);
for t = 1:M
lY_cf(:,t) = alpha(:,t).*(lK_cf(:,t)-lalpha(:,t)) + (1-alpha(:,t)).*(lE_cf(:,t)-log(1-alpha(:,t))) - log(1-sum(sigmai(:,t).*mu(:,t)))+ ...
    sum(mu(:,t).*(1-gammai(:,t)-sigmai(:,t)).*alphai(:,t).*lalphai(:,t))./(sum(mu(:,t).*(1-gammai(:,t)-sigmai(:,t)))) + ...
    sum(mu(:,t).*(1-gammai(:,t)-sigmai(:,t)).*(1-alphai(:,t)).*log(1-alphai(:,t)))./(sum(mu(:,t).*(1-gammai(:,t)-sigmai(:,t)))) + ... 
    sum(mu(:,t).*(1-gammai(:,t)-sigmai(:,t)).*log(1-gammai(:,t)-sigmai(:,t)))./(sum(mu(:,t).*(1-gammai(:,t)-sigmai(:,t)))) + ...
    sum(mu(:,t).*lA(:,t))./(sum(mu(:,t).*(1-gammai(:,t)-sigmai(:,t)))) + ...
    sum(mu(:,t).*sum(Gamma(:,:,t).*log(Gamma(:,:,t)),1)')./(sum(mu(:,t).*(1-gammai(:,t)-sigmai(:,t)))) + ...
    sum(mu(:,t).*sum(Sigma(:,:,t).*log(Sigma(:,:,t)),1)')./(sum(mu(:,t).*(1-gammai(:,t)-sigmai(:,t)))) + ...
    (sum(beta(:,t).*lbeta(:,t)) - sum(mu(:,t).*sum(Sigma(:,:,t).*lwp_cf(:,1),1)'))./(sum(mu(:,t).*(1-gammai(:,t)-sigmai(:,t))));
end
Y_cf = exp(lY_cf);




%2) solve for quantities and factor usage of each industry

lki_cf = log(alphai) + log(1-gammai-sigmai) + log(mu) + lY_cf - lwK_cf;
lei_cf = log(1-alphai) + log(1-gammai-sigmai) + log(mu) + lY_cf - lwE_cf;

ki_cf = exp(lki_cf);
ei_cf = exp(lei_cf);



% calculate qi: 1) and 2)
% 1) define V vector 

V = zeros (N,M);
for t = 1:M
V(:,t) = lA(:,t) + lmu(:,t) - (1-gammai(:,t)).*log(1-sum(sigmai(:,t).*mu(:,t))) + ...
    lK_cf(:,t).*(sigmai(:,t).*alpha(:,t) + (1-gammai(:,t)-sigmai(:,t)).*alphai(:,t)) + ...
    lE_cf(:,t).*(sigmai(:,t).*(1-alpha(:,t)) + (1-gammai(:,t)-sigmai(:,t)).*(1-alphai(:,t))) + ...
    - sigmai(:,t).*(alpha(:,t).*lalpha(:,t) + (1-alpha(:,t)).*log(1-alpha(:,t))) + ...
    (1-gammai(:,t)-sigmai(:,t)).*alphai(:,t).*(lalphai(:,t)-lalpha(:,t)) + ...
    (1-gammai(:,t)-sigmai(:,t)).*(1-alphai(:,t)).*(log(1-alphai(:,t))-log(1-alpha(:,t))) + ...
    sigmai(:,t).*(sum(mu(:,t).*(1-gammai(:,t)-sigmai(:,t)).*alphai(:,t).*lalphai(:,t))./(sum(mu(:,t).*(1-gammai(:,t)-sigmai(:,t))))) + ...
    sigmai(:,t).*(sum(mu(:,t).*(1-gammai(:,t)-sigmai(:,t)).*(1-alphai(:,t)).*log(1-alphai(:,t)))./(sum(mu(:,t).*(1-gammai(:,t)-sigmai(:,t))))) + ... 
    sigmai(:,t).*(sum(mu(:,t).*(1-gammai(:,t)-sigmai(:,t)).*log(1-gammai(:,t)-sigmai(:,t)))./(sum(mu(:,t).*(1-gammai(:,t)-sigmai(:,t))))) + ...
    sigmai(:,t).*(sum(mu(:,t).*lA(:,t))./(sum(mu(:,t).*(1-gammai(:,t)-sigmai(:,t))))) + ...
    sigmai(:,t).*(sum(mu(:,t).*sum(Gamma(:,:,t).*log(Gamma(:,:,t)),1)')./(sum(mu(:,t).*(1-gammai(:,t)-sigmai(:,t))))) + ...
    sigmai(:,t).*(sum(mu(:,t).*sum(Sigma(:,:,t).*log(Sigma(:,:,t)),1)')./(sum(mu(:,t).*(1-gammai(:,t)-sigmai(:,t))))) + ...
    sigmai(:,t).*((sum(beta(:,t).*lbeta(:,t)) - sum(mu(:,t).*sum(Sigma(:,:,t).*lwp_cf(:,1),1)'))./(sum(mu(:,t).*(1-gammai(:,t)-sigmai(:,t))))) + ...
    (1-gammai(:,t)-sigmai(:,t)).*log(1-gammai(:,t)-sigmai(:,t)) + ...
    sum(Gamma(:,:,t).*log(Gamma(:,:,t)),1)' - sum(Gamma(:,:,t).*log(mu(:,t)),1)' + ...
    sum(Sigma(:,:,t).*log(Sigma(:,:,t)),1)' - sum(Sigma(:,:,t).*lwp_cf(:,1),1)';
end




% 2) find q = (1-Gamma')^(-1)*V can be written as
lq_cf = zeros (N,M);
for t = 1:M
lq_cf(:,t) = (eye(N)-Gamma(:,:,t)')\V(:,t);
end

q_cf = exp(lq_cf);


%calculate dji = (Gamma_ji * mu_i * q_j)/mu_j
d_cf = zeros (N,N,M);
for t = 1:M
d_cf(:,:,t) = Gamma(:,:,t).*mu(:,t)'./mu(:,t).*q_cf(:,t);
end

%calculate fji 
lf_cf = zeros (N,N,M);
for t = 1:M
lf_cf(:,:,t) = lSigma(:,:,t) + lmu(:,t)' + lY_cf(:,t) - lwp_obs(:,1);
end

f_cf = exp(lf_cf);



%calculate yi = (beta_i * q_i)/mu_i
y_cf = beta.*q_cf./mu;


% TEST market clearing
if q_cf < y_cf
    disp('q smaller than y')
end


% y_cf has to be positive
% sum(d_obs,2): total production II / sum over columns, same row 
y_cf_MC = zeros (N,M);
for t = 1:M
y_cf_MC(:,t) = q_cf(:,t) - sum(d_cf(:,:,t),2);
end

if min(y_cf_MC) <= 0
    disp('q smaller than d+f')
end

% TEST: difference close to zero
if max(abs(y_cf_MC - y_cf)) >= 1e-08 
    disp('check y')
end



%calculate pi: 1) and 2)
% 1) define W vector 

W = zeros (N,M);
for t = 1:M
W(:,t) = -lA(:,t) + ...
    lK_cf(:,t).*(1-gammai(:,t)-sigmai(:,t)).*(alpha(:,t) - alphai(:,t)) + ...
    lE_cf(:,t).*(1-gammai(:,t)-sigmai(:,t)).*((1-alpha(:,t)) - (1-alphai(:,t))) + ...
    -lalpha(:,t).*(1-gammai(:,t)-sigmai(:,t)).*(alpha(:,t) - alphai(:,t)) + ...
    -log(1-alpha(:,t)).*(1-gammai(:,t)-sigmai(:,t)).*((1-alpha(:,t)) - (1-alphai(:,t))) + ...
    -(1-gammai(:,t)-sigmai(:,t)).*alphai(:,t).*lalphai(:,t) + ...
    -(1-gammai(:,t)-sigmai(:,t)).*(1-alphai(:,t)).*log(1-alphai(:,t)) + ...
    (1-gammai(:,t)-sigmai(:,t)).*(sum(mu(:,t).*(1-gammai(:,t)-sigmai(:,t)).*alphai(:,t).*lalphai(:,t))./(sum(mu(:,t).*(1-gammai(:,t)-sigmai(:,t))))) + ...
    (1-gammai(:,t)-sigmai(:,t)).*(sum(mu(:,t).*(1-gammai(:,t)-sigmai(:,t)).*(1-alphai(:,t)).*log(1-alphai(:,t)))./(sum(mu(:,t).*(1-gammai(:,t)-sigmai(:,t))))) + ... 
    (1-gammai(:,t)-sigmai(:,t)).*(sum(mu(:,t).*(1-gammai(:,t)-sigmai(:,t)).*log(1-gammai(:,t)-sigmai(:,t)))./(sum(mu(:,t).*(1-gammai(:,t)-sigmai(:,t))))) + ...
    (1-gammai(:,t)-sigmai(:,t)).*(sum(mu(:,t).*lA(:,t))./(sum(mu(:,t).*(1-gammai(:,t)-sigmai(:,t))))) + ...
    (1-gammai(:,t)-sigmai(:,t)).*(sum(mu(:,t).*sum(Gamma(:,:,t).*log(Gamma(:,:,t)),1)')./(sum(mu(:,t).*(1-gammai(:,t)-sigmai(:,t))))) + ...
    (1-gammai(:,t)-sigmai(:,t)).*(sum(mu(:,t).*sum(Sigma(:,:,t).*log(Sigma(:,:,t)),1)')./(sum(mu(:,t).*(1-gammai(:,t)-sigmai(:,t))))) + ...
    (1-gammai(:,t)-sigmai(:,t)).*((sum(beta(:,t).*lbeta(:,t)) - sum(mu(:,t).*sum(Sigma(:,:,t).*lwp_cf(:,1),1)'))./(sum(mu(:,t).*(1-gammai(:,t)-sigmai(:,t))))) + ...
    -(1-gammai(:,t)-sigmai(:,t)).*log(1-gammai(:,t)-sigmai(:,t)) + ...
    -sum(Gamma(:,:,t).*log(Gamma(:,:,t)),1)' + ...
    -sum(Sigma(:,:,t).*log(Sigma(:,:,t)),1)' + sum(Sigma(:,:,t).*lwp_cf(:,1),1)';
end


% 2) find p = (1-Gamma')^(-1)*W can be written as
lp_cf = zeros (N,M);
for t = 1:M
lp_cf(:,t) = (eye(N)-Gamma(:,:,t)')\W(:,t);
end


p_cf = exp(lp_cf);

pq_cf = p_cf.*q_cf;


%% RESULTS

% calculate effects of 1% pz on qi

dqidpz_model = zeros (N,M);
for t = 1:M
dqidpz_model(:,t) = (eye(N)-Gamma(:,:,t)')\(-sigmai(:,t).*(sum(mu(:,t).*Sigma(1,:,t)')/sum(mu(:,t).*(1-gammai(:,t)-sigmai(:,t)))) - Sigma(1,:,t)');
end

% aggregate transmission (first term)
agg_trans = zeros (N,M);
for t = 1:M
agg_trans(:,t) = -sigmai(:,t).*(sum(mu(:,t).*Sigma(1,:,t)')/sum(mu(:,t).*(1-gammai(:,t)-sigmai(:,t))));
end
 
% pz industries shares - FIRST IMPORTED GOOD (A-B)
pz_shares = zeros (N,M);
for t = 1:M
pz_shares(:,t) = Sigma(1,:,t)';
end

% share Mz/Y - EACH IMPORTED GOOD
MzY_shares = zeros (N,M);
for t = 1:M
MzY_shares(:,t) = sum(mu(:,t).*Sigma(:,:,t)');
end
% share C/Y - AGGREGATE
CY_share = zeros (1,M);
for t = 1:M
CY_share(:,t) = sum(mu(:,t).*(1-gammai(:,t)-sigmai(:,t)));
end
% share Mz/C - EACH IMPORTED GOOD
MzC_shares = zeros (N,M);
for t = 1:M
MzC_shares(:,t) = (sum(mu(:,t).*Sigma(:,:,t)')/sum(mu(:,t).*(1-gammai(:,t)-sigmai(:,t))));
end


% counterfactual
dqidpz_cf = zeros (N,M);
for t = 1:M
dqidpz_cf(:,t) = ((q_cf(:,t) - q_obs(:,t))./q_obs(:,t))*100;
end

if max(abs(dqidpz_model - dqidpz_cf)) >= 1e-2
    disp('check dqidpz_model')
end


% weighted average purchases from A-B importing industries
AB_in_wi = zeros(N,N,M);
AB_in_weighted = zeros(N,M);

Sigma(isnan(Sigma))=0;
Gamma(isnan(Gamma))=0;

for m=1:M
        AB_in_wi(:,:,m)=(pz_shares(:,m)).*Gamma(:,:,m)./sum(Gamma(:,:,m));
        AB_in_weighted(:,m) = sum(AB_in_wi(:,:,m))';
end

AB_in_weighted(isnan(AB_in_weighted))=0;

industry = (1:N)';


%% Save data
% ALL VARIABLES ARE IN VECTORS (48 X 2)%

save('WIOD_IRL_KOR_pz')


%% Export calculated data


% industrial data per country
for m=1:length(M_countries) 
        IND = table(industry,dqidpz_model(:,m),dqidpz_cf(:,m),agg_trans(:,m),pz_shares(:,m),MzY_shares(:,m),MzC_shares(:,m),AB_in_weighted(:,m));
        filepath = 'C:\Users\Veri\OneDrive - University of Essex\RESEARCH\PUBLISHING\Chapter3\Codes\Stata\Input'; 
        filename = sprintf('Industries_%s_48_pz.xls',M_countries(m));
        writetable(IND,fullfile(filepath,filename),'Sheet',1,'Range','A1','WriteVariableNames',true)
end


% all countries agg data 
ALL = table(M_countries,alpha',Y_obs',Y_cf',wK_obs',wK_cf',wE_obs',wE_cf',CY_share');
    filepath = 'C:\Users\Veri\OneDrive - University of Essex\RESEARCH\PUBLISHING\Chapter3\Codes\Stata\Input'; 
    filename = 'WIOD_2014_IRL_KOR_pz.xls';
    writetable(ALL,fullfile(filepath,filename),'Sheet',1,'Range','A1','WriteVariableNames',true)
 

 
%% PLOTS EXTRA

% 
%  % Sigma rows
%  for m=1:M
%      for i=1:N
%         y = [Sigma(i,:,m)'];
%          bar(industry,y'), xticks([1:48]), xtickangle(90)
%          set(findall(gcf,'-property','FontSize'),'FontSize',12)
%          xlabel('Buying Industries'), ylabel('import shares')
%          title(sprintf('Import shares of good %d (rows of matrix Sigma %s)',industry(i),M_countries(m)),'FontSize',14)
%          x0=50, y0=50, width=900, height=600;
%          set(gcf,'position',[x0,y0,width,height])
%          fpath = 'C:\Users\V\Dropbox\veridiana\Matlab\WIOD_CF\Sigma';
%          filename = sprintf('%s_48_Sigma_%d',M_countries(m),industry(i));
%          saveas(gcf, fullfile(fpath, filename), 'jpeg');
%      end
%  end
%  
%  
