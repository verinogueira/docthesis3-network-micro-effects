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

load('WIOD_IRL_KOR')

%% ANALYSIS

% labour (employment) industries shares adjusted
e_shares_adj = (1-gammai).*(1-alphai);

% shock vector
e_shock = (1-gammai).*(1-alphai)+sigmai.*((1-alpha)-(1-alphai));

% linkages effects
Gamma_effect = zeros(N,N,M); 
for m=1:M
    Gamma_effect(:,:,m) = inv((eye(N)-Gamma(:,:,m)')); 
end

full_ES_i = zeros(N,M); 
for m=1:M
    full_ES_i(:,m) = (eye(N)-Gamma(:,:,m)')\e_shock(:,m); 
end


% weights
weights_matrix = zeros(N,N,M); 
for m=1:M
    weights_matrix(:,:,m) = Gamma_effect(:,:,m).*(1-gammai(:,m)'); 
end

weights_test = zeros(N,M); 
for m=1:M
    weights_test(:,m) = sum(weights_matrix(:,:,m),2); 
end


% approximation matrices
% rows have each order of approximation per industry
% 4th decimal convergence no later than 17th order
ES_effect = zeros(N,N,M);
for m=1:M
    for n = 1:N
        i=n-1;
        ES_effect(:,n,m)=AppEffect(i, N, Gamma(:,:,m), e_shock(:,m));   
    end
end 
    
% test all 1s
ES_effect_test = zeros(N,M); 
for m=1:M
    ES_effect_test(:,m) = ES_effect(:,N,m)./full_ES_i(:,m);
end


% speed
% to capture for which industries convergence takes longer
ES_effect_speed = zeros(N,N,M);
for m=1:M
    for n = 1:N % each column is an order of approximation
        ES_effect_speed(:,n,m)=ES_effect(:,n,m)./full_ES_i(:,m);   
    end
end



% error per order
% to capture for which industries convergence takes longer
ES_effect_error = zeros(N,N,M);
for m=1:M
    for n = 1:N
        ES_effect_error(:,n,m)=(full_ES_i(:,m)-ES_effect(:,n,m))./full_ES_i(:,m);   
    end
end



industry = (1:N)';

%% Export data



% speed = ratios % each column is an order of approximation
for m=1:M 
        SR = table(industry,ES_effect_speed(:,:,m));
        filepath = 'C:\Users\Veri\OneDrive - University of Essex\RESEARCH\PUBLISHING\Chapter3\Codes\Stata\Input'; 
        filename = sprintf('Speed_%s_48_E.xls',M_countries(m));
        writetable(SR,fullfile(filepath,filename),'Sheet',1,'Range','A1','WriteVariableNames',true)
end


 

%% PLOTS %%


% Figure 3.12 %
newcolors = [.3451 .3451 .3451 % dark gray rgb(88, 88, 88)
             .8902 .4941 .0000 % gold rgb(227, 126, 0)
             .5647 .2078 .2314 % red rgb(144, 53, 59)
             .3333 .4588 .1843 % green rgb(85, 117, 47)
             .1020 .2784 .4353]; % navy rgb(26, 71, 111)

colororder(newcolors)
 
 for m=1:M
     i=1
    for o = [47 3 2 1 0]
        x = (1:N)';
        c = o+1;
        y = [ES_effect(:,c,m)'];
         bar(x,y,'FaceColor','flat','EdgeColor','none')
            xticks([1:48]), xtickangle(90)
        set(gca,'TickDir','out'); % The only other option is 'in'
            set(findall(gcf,'-property','FontSize'),'FontSize',12)
            xlabel('Industries'), ylabel('real output changes')
            x0=50, y0=50, width=900, height=600;
            set(gcf,'position',[x0,y0,width,height])
         hold on
         i=i+1;
    end 
    hold off
            box off 
            %xlim([0 48]) 
            %ylim([0 1])
            line([0 49],[1 1],'Color','k')
            line([49 49],[0 1],'Color','k')
 title(sprintf('Country: %s',M_countries(m)),'FontSize',14);
 legend({'Total Effect','Third Order','Second Order','First Order','Zeroth Order'},'Location','southeast','EdgeColor','none')
 filepath = 'C:\Users\Veri\OneDrive - University of Essex\RESEARCH\PUBLISHING\Chapter3\Codes\Exhibits'; 
 filename = sprintf('Approx_%s_48_overlay_E_colo_Alex',M_countries(m));
 saveas(gcf, fullfile(filepath, filename), 'epsc'); %eps format higher quality
end




%% Not used 

% % Gamma_effect plot
% for m=1:M
% imagesc(Gamma_effect(:,:,m));
% xlabel('Columns','fontsize',18), ylabel('Rows','fontsize',18);
% colorbar;
% colormap (flipud(bone(64)));
% xticks([1:48]), xtickangle(90)
% yticks([1:48])
%          set(findall(gcf,'-property','FontSize'),'FontSize',12)
%          title(sprintf('IO effects: matrix [ I - Gamma^I]^{-1} (%s)',M_countries(m)),'FontSize',14)
%          x0=50, y0=50, width=900, height=600;
%          set(gcf,'position',[x0,y0,width,height])
%          fpath = 'C:\Users\V\Dropbox\veridiana\Matlab\WIOD_CF';
%          filename = sprintf('Gamma_effect_%s_48',M_countries(m));
%          saveas(gcf, fullfile(fpath, filename), 'jpeg');
% end
% 
% 
% % Weights matrix
% for m=1:M
% imagesc(weights_matrix(:,:,m));
% xlabel('Columns','fontsize',18), ylabel('Rows','fontsize',18);
% colorbar;
% colormap (flipud(bone(64)));
% xticks([1:48]), xtickangle(90)
% yticks([1:48])
%          set(findall(gcf,'-property','FontSize'),'FontSize',12)
%          title(sprintf('Weights matrix [ I - Gamma^I]^{-1}.(1-gamma) (%s)',M_countries(m)),'FontSize',14)
%          x0=50, y0=50, width=900, height=600;
%          set(gcf,'position',[x0,y0,width,height])
%          fpath = 'C:\Users\V\Dropbox\veridiana\Matlab\WIOD_CF';
%          filename = sprintf('weights_matrix_%s_48',M_countries(m));
%          saveas(gcf, fullfile(fpath, filename), 'jpeg');
% end
% 
%   
% 
% 
% % plot approximation
%  for m=1:M
% for n = 1:N
% x = (1:N)';
% y = [ES_effect(n,:,m)];
% bar(x,y)
%          xticks([1:48]), xtickangle(90)
%          set(findall(gcf,'-property','FontSize'),'FontSize',12)
%          xlabel('Order of approximation (-1)'), ylabel('Calculated effect of E CF')
%          title(sprintf('Approximation for Industry %d (%s)',n,M_countries(m)),'FontSize',14), 
%          x0=50, y0=50, width=900, height=600;
%          set(gcf,'position',[x0,y0,width,height])
%          fpath = 'C:\Users\V\Dropbox\veridiana\Matlab\WIOD_CF\Approx';
%          filename = sprintf('Approx_%s_48_%d',M_countries(m),n);
%          saveas(gcf, fullfile(fpath, filename), 'jpeg');
% end
%  end
% 
%  
%  
%  
%  
% % plot approximation overlay B&W
%  for m=1:M
%      i=1
% for o = [47 3 2 1 0]
% x = (1:N)';
% c = o+1;
% l = [1 1 1]
% y = [ES_effect(:,c,m)'];
%          bar(x,y,1.15-0.2*i,'FaceColor',l/c,'EdgeColor','none')
%          xticks([1:48]), xtickangle(90)
%          set(findall(gcf,'-property','FontSize'),'FontSize',12)
%             xlabel('Industries'), ylabel('real output changes')
%          x0=50, y0=50, width=900, height=600;
%          set(gcf,'position',[x0,y0,width,height])
%          hold on
%          i=i+1;
% end 
% hold off
%  title(sprintf('Country: %s',M_countries(m)),'FontSize',14);
%  legend({'Total Effect','Third Order','Second Order','First Order','Zeroth Order'},'Location','southeast')
%  fpath = 'C:\Users\Veri\Dropbox\veridiana\Matlab\WIOD_CF\Approx';
%  filename = sprintf('Approx_%s_48_overlay_E',M_countries(m));
%  saveas(gcf, fullfile(fpath, filename), 'epsc'); %eps format higher quality
%  end
%  
%  
%  
%  
% % plot approximation overlay colourful         
% for m=1:M
%      i=1
%     for o = [47 3 2 1 0]
%         x = (1:N)';
%         c = o+1;
%         y = [ES_effect(:,c,m)'];
%          bar(x,y,1.15-0.2*i,'FaceColor','flat','EdgeColor','none')
%             xticks([1:48]), xtickangle(90)
%             set(findall(gcf,'-property','FontSize'),'FontSize',12)
%             xlabel('Industries'), ylabel('real output changes')
%             x0=50, y0=50, width=900, height=600;
%             set(gcf,'position',[x0,y0,width,height])
%          hold on
%          i=i+1;
%     end 
%     hold off
%  title(sprintf('Country: %s',M_countries(m)),'FontSize',14);
%  legend({'Total Effect','Third Order','Second Order','First Order','Zeroth Order'},'Location','southeast')
%  fpath = 'C:\Users\Veri\Dropbox\veridiana\Matlab\WIOD_CF\Approx';
%  filename = sprintf('Approx_%s_48_overlay_E_colo',M_countries(m));
%  saveas(gcf, fullfile(fpath, filename), 'epsc'); %eps format higher quality
% end
% 
% 
%  
% % plot approximation overlay B&W
%  for m=1:M
%      i=1
% for o = [47 3 2 1 0]
% x = (1:N)';
% c = o+1;
% l = [0.9 0.9 0.9]
% y = [ES_effect(:,c,m)'];
%          bar(x,y,'FaceColor',l/c,'EdgeColor','none')
%          xticks([1:48]), xtickangle(90)
%          set(findall(gcf,'-property','FontSize'),'FontSize',12)
%             xlabel('Industries'), ylabel('real output changes')
%          x0=50, y0=50, width=900, height=600;
%          set(gcf,'position',[x0,y0,width,height])
%          hold on
%          i=i+1;
% end 
% hold off
%  title(sprintf('Country: %s',M_countries(m)),'FontSize',14);
%  legend({'Total Effect','Third Order','Second Order','First Order','Zeroth Order'},'Location','southeast')
%  fpath = 'C:\Users\Veri\Dropbox\veridiana\Matlab\WIOD_CF\Approx';
%  filename = sprintf('Approx_%s_48_overlay_E_Alex',M_countries(m));
%  saveas(gcf, fullfile(fpath, filename), 'epsc'); %eps format higher quality
%  end
%   
%  
% % plot speed convergence
%  for m=1:M
% for n = 1:12
%     i=n-1;
% x = (1:N)';
% y = [ES_effect_speed(:,n,m)];
% bar(x,y,'FaceColor','#ffda33','EdgeColor','#ffda33')
% xticks([1:48]), xtickangle(90)
% %xticklabels({'A-B','C10-C12','C13-C15','C16-C18','C19','C20-C23','C24-C25','C26-C27','C28','C29-C30','C31-C33','D-E','F','G45','G46','G47','H49-H53','I','J58-J60','J61','J62-J63','K','L',	'M-N','O','P','Q','R-S','T'})
% ylim([0 1])
% grid on
%          set(findall(gcf,'-property','FontSize'),'FontSize',12)
% xlabel('Industries'), ylabel('Approximated value over total effect')
%         title(sprintf('Speed of convergence E CF - (%d)th order (%s)',i,M_countries(m)),'FontSize',14), 
%          x0=50, y0=50, width=900, height=600;
%          set(gcf,'position',[x0,y0,width,height])
%          set(gcf,'color','w');
%          fpath = 'C:\Users\V\Dropbox\veridiana\Matlab\WIOD_CF\Speed';
%          filename = sprintf('Speed_%s_48_%d',M_countries(m),n);
%          saveas(gcf, fullfile(fpath, filename), 'jpeg');
%    end
%  end
% 
% 
% 
% % plot approximation error
%  for m=1:M
% for n = 1:N
% x = (1:N)';
% y = [ES_effect_error(n,:,m)];
% bar(x,y)
% xlim([0 17])
% ylim([0 1])
% xlabel('Order of approximation (-1)'), ylabel('Error relative to total effect of ES CF')
%          set(findall(gcf,'-property','FontSize'),'FontSize',12)
%          title(sprintf('Approximation error for Industry %d (%s)',n,M_countries(m)),'FontSize',14), 
%          x0=50, y0=50, width=900, height=600;
%          set(gcf,'position',[x0,y0,width,height])
%          set(gcf,'color','w');
%          fpath = 'C:\Users\V\Dropbox\veridiana\Matlab\WIOD_CF\Error';
%          filename = sprintf('Error_%s_48_%d',M_countries(m),n);
%          saveas(gcf, fullfile(fpath, filename), 'jpeg');
% end
%  end
% 
%  
%  