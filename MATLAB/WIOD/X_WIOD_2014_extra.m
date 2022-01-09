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

%% LOAD DATA 
%Load 2014 World Input Output Data (WIOD)
load WIOD.mat


%% Extra (not used in main body)


% weighted average purchases from capital intensive industries
H_in_wi = zeros(N,N,M);
H_in_weighted = zeros(N,M);

alphai(isnan(alphai))=0;
Gamma(isnan(Gamma))=0;

for m=1:M
        H_in_wi(:,:,m)=alphai(:,m).*Gamma(:,:,m)./sum(Gamma(:,:,m));
        H_in_weighted(:,m) = sum(H_in_wi(:,:,m))';
end

H_in_weighted(isnan(H_in_weighted))=0;


% weighted average imports from (domestically) high-skill intensive industries
H_in_M_wi = zeros(N,N,M);
H_in_M_weighted = zeros(N,M);



for m=1:M
    H_in_M_wi(:,:,m)=alphai(:,m).*Sigma(:,:,m)./sum(Sigma(:,:,m));
    H_in_M_weighted(:,m) = sum(H_in_M_wi(:,:,m))';
end

H_in_M_weighted(isnan(H_in_M_weighted))=0;


% Gamma per country
% for i=1:length(M_countries) 
%         CC = table(industry,Gamma(:,:,i));
%         filename = sprintf('Gamma_%s_2014.xls',M_countries(i));
%         writetable(CC,filename,'Sheet',1,'Range','B2')
% end
for i=1:length(M_countries) 
        CC = table(Gamma(:,:,i));
        filename = sprintf('Gamma_%s_2014.csv',M_countries(i));
        writetable(CC,filename,'WriteVariableNames',false)
end

for i=1:length(M_countries) 
        CC = table(Gamma(:,:,i));
        filename = sprintf('Gamma_%s_2014.xls',M_countries(i));
        writetable(CC,filename,'Sheet',1,'Range','A1','WriteVariableNames',true)
end


% Sigma per country
for i=1:length(M_countries) 
        CC = table(Sigma(:,:,i));
        filename = sprintf('Sigma_%s_2014.csv',M_countries(i));
        writetable(CC,filename,'WriteVariableNames',false)
end

for i=1:length(M_countries) 
        CC = table(Sigma(:,:,i));
        filename = sprintf('Sigma_%s_2014.xls',M_countries(i));
        writetable(CC,filename,'Sheet',1,'Range','A1','WriteVariableNames',true)
end



% industrial data per country
for m=1:length(M_countries) 
        IND = table(industry,alphai(:,m),beta(:,m),mu(:,m),gammai(:,m),sigmai(:,m),H_in_weighted(:,m),H_in_M_weighted(:,m),py_obs(:,m));
        filename = sprintf('Industries_%s_2014.xls',M_countries(m));
        writetable(IND,filename,'Sheet',1,'Range','A1','WriteVariableNames',true)
end

 
 % IO plot
 for i=1:length(M_countries)
imagesc(Gamma(:,:,i));
title(sprintf('Domestic IO Matrix - Share of input purchases in total cost (%s 2014)',M_countries(i)));
xlabel('Buying industries'), ylabel('Selling industries');
colorbar;
colormap (flipud(pink));
axis xy;
xticks([1:56])
yticks([1:56])
saveas(gcf,sprintf('%s_2014_IO.png',M_countries(i)));
end

 
 % II shares plot
 for i=1:length(M_countries)
imagesc(Sigma(:,:,i));
title(sprintf('Imported II Matrix - Share of imports purchases in total cost (%s 2014)',M_countries(i)));
xlabel('Buying industries'), ylabel('Selling industries');
colorbar;
colormap (flipud(pink));
axis xy;
xticks([1:56])
yticks([1:56])
saveas(gcf,sprintf('%s_2014_M.png',M_countries(i)));
end

 
 
