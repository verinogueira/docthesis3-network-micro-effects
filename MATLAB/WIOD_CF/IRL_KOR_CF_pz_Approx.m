%% WIOD, 2014, 48 industries, 2 countries
% approximations of the labour supply shock
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

load('WIOD_IRL_KOR_pz')

%% ANALYSIS

% AB import shares
pz_shares;

% shock vector
pz_shock = agg_trans-pz_shares;

% linkages effects
Gamma_effect = zeros(N,N,M); 
for m=1:M
    Gamma_effect(:,:,m) = inv((eye(N)-Gamma(:,:,m)')); 
end

full_pz_i = zeros(N,M); 
for m=1:M
    full_pz_i(:,m) = (eye(N)-Gamma(:,:,m)')\pz_shock(:,m); 
end



% approximation matrices
% rows have each order of approximation per industry
% 4th decimal convergence no later than 17th order
pz_effect = zeros(N,N,M);
for m=1:M
    for n = 1:N
        i=n-1;
        pz_effect(:,n,m)=AppEffect(i, N, Gamma(:,:,m), pz_shock(:,m));   
    end
end 
   

% test all 1s
pz_effect_test = zeros(N,M); 
for m=1:M
    pz_effect_test(:,m) = pz_effect(:,N,m)./full_pz_i(:,m);
end


% speed
% to capture for which industries convergence takes longer
pz_effect_speed = zeros(N,N,M);
for m=1:M
    for n = 1:N % each column is an order of approximation
        pz_effect_speed(:,n,m)=pz_effect(:,n,m)./full_pz_i(:,m);   
    end
end



% error per order
% to capture for which industries convergence takes longer
pz_effect_error = zeros(N,N,M);
for m=1:M
    for n = 1:N
        pz_effect_error(:,n,m)=(full_pz_i(:,m)-pz_effect(:,n,m))./full_pz_i(:,m);   
    end
end



industry = (1:N)';


 

%% PLOTS



newcolors = [.3451 .3451 .3451 % dark gray rgb(88, 88, 88)
             .8902 .4941 .0000 % gold rgb(227, 126, 0)
             .5647 .2078 .2314 % red rgb(144, 53, 59)
             .3333 .4588 .1843 % green rgb(85, 117, 47)
             .1020 .2784 .4353]; % navy rgb(26, 71, 111)

colororder(newcolors)


% plot approximation overlay colourful FINAL VERSION        
for m=1:M
     i=1
    for o = [47 3 2 1 0]
        x = (1:N)';
        c = o+1;
        y = [pz_effect(:,c,m)'];
         bar(x,y,'FaceColor','flat','EdgeColor','none')
            xticks([1:48]), xtickangle(90)
        %set(gca,'TickDir','out','XAxisLocation','top'); % The only other option is 'in'
            set(findall(gcf,'-property','FontSize'),'FontSize',12)
            xlabel('Industries'), ylabel('real output changes')
            x0=50, y0=50, width=900, height=600;
            set(gcf,'position',[x0,y0,width,height])
         hold on
         i=i+1;
    end 
    hold off
            %box off 
            %xlim([0 48]) 
            %ylim([0 1])
            %line([0 49],[-1 -1],'Color','k')
            %line([49 49],[0 -1],'Color','k')
 title(sprintf('Country: %s',M_countries(m)),'FontSize',14);
 legend({'Total Effect','Third Order','Second Order','First Order','Zeroth Order'},'Location','southeast','EdgeColor','none')
 filepath = 'C:\Users\Veri\OneDrive - University of Essex\RESEARCH\PUBLISHING\Chapter3\Codes\Exhibits'; 
 filename = sprintf('Approx_%s_48_overlay_pz2_colo_Alex',M_countries(m));
 saveas(gcf, fullfile(filepath, filename), 'epsc'); %eps format higher quality
 %saveas(gcf, fullfile(fpath, filename), 'jpeg'); 
end
