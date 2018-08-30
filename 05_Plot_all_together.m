%% Initialization and data loading

%% Calculation of decussation index for COMMIT output
cd('/home/auguser2016/Projects/Leicester_vs_LiFE/COMMIT')
mkdir('Plots')

% Researched groups and participants
achiasma={['hw91'];['ps94']};
albinism={['ce04'];['fe21'];['ib57'];['kw99'];['nh50'];['rx88'];['sj22'];['tq63'];['uh47']};
controls={['la21'];['lw37'];['nb30'];['ow93'];['uf97'];['xn78'];['xs62'];['ta14']};

% Analyzed models
all_models={'LiFE_Model', 'StickZeppelinBall_Model'} 

% Type of input tracks
%all_tracks={'CSD_seed';'CSD_selected';'DT_seed';'DT_selected'}
all_tracks={'CSD_seed'; 'DT_seed'}

%variant='CSD_seed'
%filtering='noLiFE'
all_combinations=cell(1,1)
for r=1:length(all_models)
    for s=1:length(all_tracks)
    all_combinations{1,(r-1)*length(all_tracks)+s} =strcat(all_models{r},'_',all_tracks{s})
    end
end


% create table with decussation indices for albinism and controls for every
% approach
con_dec_ind = array2table(nan(length(all_combinations), length(controls)));
con_dec_ind.Properties.VariableNames = controls;
con_dec_ind.Properties.RowNames = all_combinations';

alb_dec_ind = array2table(nan(length(all_combinations), length(albinism)));
alb_dec_ind.Properties.VariableNames = albinism;
alb_dec_ind.Properties.RowNames = all_combinations';

ach_dec_ind = array2table(nan(length(all_combinations), length(achiasma)));
ach_dec_ind.Properties.VariableNames = achiasma;
ach_dec_ind.Properties.RowNames = all_combinations';

for xx=1:length(all_models)
    for yy=1:length(all_tracks)
        
        variant=all_models{xx};
        filtering=all_tracks{yy};
        
        load(strcat('COMMIT_',variant,'_',filtering,'.mat'))

%% Plot all groups using gscatter

names = con_ips.Properties.VariableNames;

contra = table2array(con_ips({'Left to Right'},:)) + table2array(con_ips({'Right to Left'},:))
ipsi = table2array(con_ips({'Left to Left'},:)) + table2array(con_ips({'Right to Right'},:))

colors=zeros(size(con_ips,2),3);
group=zeros(1,size(con_ips,2))

for i=1:size(con_ips,2)
    if sum(strcmp(achiasma,names(i)))
        colors(i,:)=[0 0 1];
        group(1,i)=1
    end
    if sum(strcmp(albinism,names(i)))
        colors(i,:)=[1 0 0];
        group(1,i)=2
    end
    if sum(strcmp(controls,names(i)))
        colors(i,:)=[0 1 0];
        group(1,i)=3
    end
end

nonan_ipsi=~isnan(ipsi)
nonan_contra=~isnan(contra)

figure('units','normalized','outerposition',[0 0 1 1])

fig = gscatter(ipsi(nonan_ipsi),contra(nonan_contra), group(nonan_ipsi), 'kkk', 'ooo', 15, 'filled') %54,colors(nonan_contra,:),'filled')
set(fig(1), 'MarkerFaceColor', 'b')
set(fig(2), 'MarkerFaceColor', 'r')
set(fig(3), 'MarkerFaceColor', 'g')

fig=gcf
hold on

scale=1.1*max(max(contra),max(ipsi))

plot(linspace(0,scale), linspace(0,scale),'k-.')

axis([0 scale 0 scale])

title(strcat('Crossing for ', 32, variant,'_',filtering), 'Fontsize', 20, 'Interpreter', 'none', 'Fontsize', 40)
xlabel('Number of non-crossing streamlines', 'Fontsize', 30) 
ylabel('Number of crossing streamlines', 'Fontsize', 30)
legend({'Achiasma (n=2)', 'Albinism   (n=9)','Controls   (n=8)','50:50 ratio'},'Location','northwest', 'Fontsize', 30)
%legend( fig([1 4 ], {'OVGU Controls (n=8)', 'HCP Controls (n=7)'}, 'Location', 'northwest', 'Fontsize', 30))
pbaspect([1 1 1])

text(ipsi + scale/100, contra + scale/100, names, 'Fontsize', 12)

%saveas(fig, strcat('Plots/COMMIT_',variant,'_', filtering,'.jpg'))
close

%    end
%end

%% Plot deccusation indices in form of bar plots

% calculate decussation index
decussation_index=contra./(contra+ipsi)*100

% separate data into groups by creating true-false vector for each group
% and extraction of corresponding values and names
control_vec=[group==3]
control_name=names(control_vec)

achiasma_vec=[group==1]
achiasma_name=names(achiasma_vec)

albinism_vec=[group==2]
albinism_name=names(albinism_vec)

% create dynamic name to mark row to be filled
writeto=strcat(variant,'_',filtering)

% write values to corresponding tables
con_dec_ind({writeto},:)=array2table(decussation_index(control_vec))
alb_dec_ind({writeto},:)=array2table(decussation_index(albinism_vec))
ach_dec_ind({writeto},:)=array2table(decussation_index(achiasma_vec))

% plot groups
%bar_plot=[decussation_index(achiasma_vec) 0 decussation_index(control_vec) 0 decussation_index(albinism_vec) ]
%bar_names={names{achiasma_vec} '' names{control_vec} '' names{albinism_vec}}
    
%bar(bar_plot)

% add names
%set(gca, 'Xtick', 1:1:length(bar_names))
%set(gca, 'xticklabel', bar_names)

% color bars

% significance test (p-value)
%[h,p]=ttest2(decussation_index(control_vec), decussation_index(albinism_vec))
%[h,p]=ttest2(decussation_index(albinism_vec), decussation_index(achiasma_vec))
%[h,p]=ttest2(decussation_index(achiasma_vec), decussation_index(control_vec))

%saveas(gcf, strcat('barplot_',variant,'_', filtering_name,'.jpg'))
%close
%clf

    end
end

old_combinations=all_combinations
con_dec_ind_COMMIT = con_dec_ind
alb_dec_ind_COMMIT = alb_dec_ind
ach_dec_ind_COMMIT = ach_dec_ind


%% Calculation of decussation index for non-COMMIT output

cd('/home/auguser2016/Projects/Leicester_vs_LiFE')


% 3 different approaches
%all_variants={'CSD_seed';'CSD_selected';'DT_seed';'DT_selected'}%;'DT2_seed';'DT2_selected'}
%all_variants={'CSD_seed';'CSD_selected';'DT_seed';'DT_selected'}
all_variants={'CSD_seed', 'DT_seed'}
% Leicester or LiFE
all_filterings={'noLiFE', 'postLiFE'}

%variant='CSD_seed'
%filtering='noLiFE'
all_combinations=cell(1,1)
for r=1:length(all_variants)
    for s=1:length(all_filterings)
    %all_combinations{1,2*(r-1)+1} =strcat(all_variants{r},'_noLiFE')
    %all_combinations{1,2*r} =strcat(all_variants{r},'_postLiFE')
    
    all_combinations{1,(s-1)*length(all_variants)+r} =strcat(all_variants{r},'_',all_filterings{s})
    
    end
end


% create table with decussation indices for albinism and controls for every
% approach
con_dec_ind = array2table(nan(length(all_combinations), length(controls)));
con_dec_ind.Properties.VariableNames = controls;
con_dec_ind.Properties.RowNames = all_combinations';

alb_dec_ind = array2table(nan(length(all_combinations), length(albinism)));
alb_dec_ind.Properties.VariableNames = albinism;
alb_dec_ind.Properties.RowNames = all_combinations';

ach_dec_ind = array2table(nan(length(all_combinations), length(achiasma)));
ach_dec_ind.Properties.VariableNames = achiasma;
ach_dec_ind.Properties.RowNames = all_combinations';

for xx=1:length(all_variants)
    for yy=1:length(all_filterings)
        
        variant=all_variants{xx};
        filtering=all_filterings{yy};

        load(strcat(variant,'_',filtering,'.mat'))

%% Plot all groups using gscatter

names = con_ips.Properties.VariableNames;

contra = table2array(con_ips({'Left to Right'},:)) + table2array(con_ips({'Right to Left'},:))
ipsi = table2array(con_ips({'Left to Left'},:)) + table2array(con_ips({'Right to Right'},:))

colors=zeros(size(con_ips,2),3);
group=zeros(1,size(con_ips,2))

for i=1:size(con_ips,2)
    if sum(strcmp(achiasma,names(i)))
        colors(i,:)=[0 0 1];
        group(1,i)=1
    end
    if sum(strcmp(albinism,names(i)))
        colors(i,:)=[1 0 0];
        group(1,i)=2
    end
    if sum(strcmp(controls,names(i)))
        colors(i,:)=[0 1 0];
        group(1,i)=3
    end
end

nonan_ipsi=~isnan(ipsi)
nonan_contra=~isnan(contra)

figure('units','normalized','outerposition',[0 0 1 1])

fig = gscatter(ipsi(nonan_ipsi),contra(nonan_contra), group(nonan_ipsi), 'kkk', 'ooo', 15, 'filled') %54,colors(nonan_contra,:),'filled')
set(fig(1), 'MarkerFaceColor', 'b')
set(fig(2), 'MarkerFaceColor', 'r')
set(fig(3), 'MarkerFaceColor', 'g')

fig=gcf
hold on

scale=1.1*max(max(contra),max(ipsi))

plot(linspace(0,scale), linspace(0,scale),'k-.')

axis([0 scale 0 scale])

title(strcat('Crossing for ',32,variant,'_',filtering), 'Fontsize', 20, 'Interpreter', 'none', 'Fontsize', 40)
xlabel('Number of non-crossing streamlines', 'Fontsize', 30) 
ylabel('Number of crossing streamlines', 'Fontsize', 30)
legend({'Achiasma (n=2)', 'Albinism   (n=9)','Controls   (n=8)','50:50 ratio'},'Location','northwest', 'Fontsize', 30)
%legend( fig([1 4 ], {'OVGU Controls (n=8)', 'HCP Controls (n=7)'}, 'Location', 'northwest', 'Fontsize', 30))
pbaspect([1 1 1])

text(ipsi + scale/100, contra + scale/100, names, 'Fontsize', 12)

%saveas(fig, strcat(variant,'_', filtering_name,'.jpg'))
close
%clf

%    end
%end

%% Plot deccusation indices in form of bar plots

% calculate decussation index
decussation_index=contra./(contra+ipsi)*100

% separate data into groups by creating true-false vector for each group
% and extraction of corresponding values and names
control_vec=[group==3]
control_name=names(control_vec)

achiasma_vec=[group==1]
achiasma_name=names(achiasma_vec)

albinism_vec=[group==2]
albinism_name=names(albinism_vec)

% create dynamic name to mark row to be filled
writeto=strcat(variant,'_',filtering)

% write values to corresponding tables
con_dec_ind({writeto},:)=array2table(decussation_index(control_vec))
alb_dec_ind({writeto},:)=array2table(decussation_index(albinism_vec))
ach_dec_ind({writeto},:)=array2table(decussation_index(achiasma_vec))

% plot groups
%bar_plot=[decussation_index(achiasma_vec) 0 decussation_index(control_vec) 0 decussation_index(albinism_vec) ]
%bar_names={names{achiasma_vec} '' names{control_vec} '' names{albinism_vec}}
    
%bar(bar_plot)

% add names
%set(gca, 'Xtick', 1:1:length(bar_names))
%set(gca, 'xticklabel', bar_names)

% color bars

% significance test (p-value)
%[h,p]=ttest2(decussation_index(control_vec), decussation_index(albinism_vec))
%[h,p]=ttest2(decussation_index(albinism_vec), decussation_index(achiasma_vec))
%[h,p]=ttest2(decussation_index(achiasma_vec), decussation_index(control_vec))

%saveas(gcf, strcat('barplot_',variant,'_', filtering,'.jpg'))
%close
%clf

    end
end

% Now properly merge two tables

ACHIASMA= vertcat(ach_dec_ind, ach_dec_ind_COMMIT)
CONTROL= vertcat(con_dec_ind, con_dec_ind_COMMIT)
ALBINISM= vertcat(alb_dec_ind, alb_dec_ind_COMMIT)

COMBINATIONS=[all_combinations, old_combinations]


%% Statistics

addpath(pwd)
% Test normality
NORMALITY=array2table(zeros(length(COMBINATIONS),2))
NORMALITY.Properties.VariableNames={'Albinism','Control'}
NORMALITY.Properties.RowNames=COMBINATIONS
for xyz=1:length(COMBINATIONS)
    
    current = COMBINATIONS(xyz)
    NORMALITY{current,'Albinism'}= kstest(ALBINISM{current,:})
    NORMALITY{current,'Control'}= kstest(CONTROL{current,:})
end   
%T=NORMALITY
%uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,'RowName',T.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);

% Create table with statistical differences
STATISTICS=array2table(zeros(length(COMBINATIONS),2))
STATISTICS.Properties.VariableNames={'Control_vs_Albinism_Wilcoxon','Control_vs_Albinism_Student'}
STATISTICS.Properties.RowNames=COMBINATIONS
for xyz=1:length(COMBINATIONS)
    current = COMBINATIONS(xyz)
    STATISTICS{current,'Control_vs_Albinism_Wilcoxon'}= ranksum(CONTROL{current,:}, ALBINISM{current,:},'alpha',0.025)
    [~,STATISTICS{current,'Control_vs_Albinism_Student'}]= ttest2(CONTROL{current,:}, ALBINISM{current,:},'alpha',0.005, 'tail', 'both')
end
%T=STATISTICS
%uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,'RowName',T.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);

%[H,P,CI,STATS]=ttest2(CONTROL{'DT_seed_noLiFE',:}, ALBINISM{'DT_seed_noLiFE',:},'alpha',0.025)


%% Plotting all the data with the statistical information
for k=1:length(COMBINATIONS)
    
    subplot(4,length(COMBINATIONS)/4,k)
    hold on
    plot_variant=COMBINATIONS{k}
    
    % create array for x axis tick - their number must match number of
    % corresponding decussation_indices values
    xaxis_alb(1:size(ALBINISM({plot_variant},:),2))=1
    % scatter dots equally between x-0.5, x+0.5
    xaxis_alb=xaxis_alb + [-0.25:0.5/(length(xaxis_alb)-1):0.25]
    
    % calculate mean only including non NaN values
    temp_albdecind=table2array(ALBINISM({plot_variant},:))
    median_albdecind=median(temp_albdecind(~isnan(temp_albdecind)))
    
    % plot scatter plot
    h1=scatter(xaxis_alb, table2array(ALBINISM({plot_variant},:)),45, 'filled','r','MarkerEdgeColor',[0,0,0])
    % plot mean as a horizontal bar
    h2=plot([0, 4],[median_albdecind, median_albdecind],'--r', 'MarkerEdgeColor',[0,0,0], 'LineWidth', 1)
    
    %controls
    xaxis_con(1:size(CONTROL({plot_variant},:),2))=2
    xaxis_con=xaxis_con + [-0.25:0.5/(length(xaxis_con)-1):0.25]
    temp_condecind=table2array(CONTROL({plot_variant},:))
    median_condecind=median(temp_condecind(~isnan(temp_condecind)))
    
    h3=scatter(xaxis_con, table2array(CONTROL({plot_variant},:)),45, 'filled','g','MarkerEdgeColor',[0,0,0])
    h4=plot([0, 4],[median_condecind, median_condecind], '--g','LineWidth', 25)
    
    % achiasma
    xaxis_ach(1:size(ACHIASMA({plot_variant},:),2))=3
    xaxis_ach=xaxis_ach + [-0.25:0.5/(length(xaxis_ach)-1):0.25]
    temp_achdecind=table2array(ACHIASMA({plot_variant},:))
    median_achdecind=median(temp_achdecind(~isnan(temp_achdecind)))
    
    h5=scatter(xaxis_ach, table2array(ACHIASMA({plot_variant},:)),45, 'filled','b','MarkerEdgeColor',[0,0,0])
    h6=plot([0, 4],[median_achdecind, median_achdecind], '--b')
    
    axis([0 4 0 100])
    set(gca,'YTick',0:25:100)
    set(gca,'YTickLabel',str2mat('0','25','50','75','100'))
    
    %{
     if [k == 1]
          title('CSD_seed', 'Fontsize', 15, 'Interpreter', 'none')
     end   
     
     if [k == 2]
          title('CSD_selected', 'Fontsize', 15, 'Interpreter', 'none')
     end 
     
     if [k == 3]
          title('DT_seed', 'Fontsize', 15, 'Interpreter', 'none')
     end   
     
     if [k == 4]
          title('DT_selected', 'Fontsize', 15, 'Interpreter', 'none')
     end 
     
    if [k == 1]
          ylabel('No filtering', 'Fontsize', 15, 'Interpreter', 'none')
    end    
     
     if [k == 5]
          ylabel('LiFE by Franco Pestilli', 'Fontsize', 15, 'Interpreter', 'none')
    end    

     if [k == 9]
          ylabel('LiFE in COMMIT', 'Fontsize', 15, 'Interpreter', 'none')
     end 
    
      if [k == 13]
          ylabel('SZB in COMMIT', 'Fontsize', 15, 'Interpreter', 'none')
      end   
    %}

    if [k ==1 ]
      legend({'\color{red}Albinism','\color{red}Albinism-median','\color{green}Control','\color{green}Control-median','\color{blue}Achiasma','\color{blue}Achiasma-median'},'Location','northeast', 'Fontsize', 10)
      legend([h1 h3 h5],{'\color{red}Albinism (--- Median)','\color{green}Control (--- Median)','\color{blue}Achiasma (--- Median)'},'Location','northeast', 'Fontsize', 10)
      
      %ylabel('Decussation index', 'Fontsize', 15, 'Interpreter', 'none')
      ylabel('No filtering', 'Fontsize', 15, 'Interpreter', 'none')
    end
    
  
    if [k == 3]
          ylabel('LiFE by Franco Pestilli', 'Fontsize', 15, 'Interpreter', 'none')
    end    

     if [k == 5]
          ylabel('LiFE in COMMIT', 'Fontsize', 15, 'Interpreter', 'none')
     end 
    
      if [k == 7]
          ylabel('SZB in COMMIT', 'Fontsize', 15, 'Interpreter', 'none')
      end   
    
     if [k == 1]
          title('CSD_seed', 'Fontsize', 15, 'Interpreter', 'none')
     end   
     
     if [k == 2]
          title('DT_seed', 'Fontsize', 15, 'Interpreter', 'none')
     end     
    
    
    if [NORMALITY{plot_variant,'Albinism'}==0 || NORMALITY{plot_variant,'Control'}==0]
        if [STATISTICS{plot_variant,'Control_vs_Albinism_Wilcoxon'}<0.025 ]
            xlabel(strcat('Wilcoxon p= ',num2str(STATISTICS{plot_variant,'Control_vs_Albinism_Wilcoxon'}),32,'(Student p= ', num2str(STATISTICS{plot_variant, 'Control_vs_Albinism_Student'}),')'),'fontsize',10, 'FontWeight','bold')
        else
             xlabel(strcat('Wilcoxon p= ',num2str(STATISTICS{plot_variant,'Control_vs_Albinism_Wilcoxon'}),32,'(Student p= ', num2str(STATISTICS{plot_variant, 'Control_vs_Albinism_Student'}),')'),'fontsize',10)
        end
    else
        if [STATISTICS{plot_variant,'Control_vs_Albinism_Student'}<0.025 || STATISTICS{plot_variant,'Control_vs_Albinism_Wilcoxon'}<0.025 ]
            xlabel(strcat('Wilcoxon p= ',num2str(STATISTICS{plot_variant,'Control_vs_Albinism_Wilcoxon'}),32,'and Student p=',num2str(STATISTICS{plot_variant,'Control_vs_Albinism_Student'})),'fontsize',10, 'FontWeight','bold')
        else
             xlabel(strcat('Wilcoxon p= ',num2str(STATISTICS{plot_variant,'Control_vs_Albinism_Wilcoxon'}),32,'and Student p=',num2str(STATISTICS{plot_variant,'Control_vs_Albinism_Student'})),'fontsize',10)
        end
    end
    
    %xlabel(plot_variant, 'Fontsize', 15, 'Interpreter', 'none')
    set(gca,'xtick',[])
    set(findall(gca, 'Type', 'Line'),'LineWidth',3);
    

    %saveas(gcf,strcat('/media/auguser2016/Tiny/Plots_for_Lubeck/Best/',plot_variant,'.jpg'))
    
    %close
    
end
pause
saveas(gcf,strcat('/home/auguser2016/Projects/Leicester_vs_LiFE/OC_crossing_500_primary.jpg'))%',plot_variant,'.jpg'))
%print(gcf,'standard.png','-dpng','-r300');        
%saveas(gcf,'/media/auguser2016/Tiny/Plots_for_Lubeck/Best/dots.jpg')


%ACHIASMA{1,:}
%ACHIASMA.('hw91')
%ACHIASMA{'CSD_seed_noLiFE',:}

%{

%% Separate data into group - achiasma, albinism and controls

% Set labels

% Group data
no_er_controls=array2table(zeros(size(no_er,1),size(controls,1)))
con_ips_controls=array2table(zeros(size(con_ips,1),size(controls,1)))

no_er_controls.Properties.RowNames=no_er.Properties.RowNames
con_ips_controls.Properties.RowNames=con_ips.Properties.RowNames

for i=1:size(controls,1)
    name=controls(i)
    namestr=name{1}

    %no_er_controls.Properties.VariableNames(i)=name
    %no_er_controls.(namestr)=no_er.(namestr)
    
    con_ips_controls.Properties.VariableNames(i)=name
    con_ips_controls.(namestr)=con_ips.(namestr)  
end

no_er_achiasma=array2table(zeros(size(no_er,1),size(achiasma,1)))
con_ips_achiasma=array2table(zeros(size(con_ips,1),size(achiasma,1)))

no_er_achiasma.Properties.RowNames=no_er.Properties.RowNames
con_ips_achiasma.Properties.RowNames=con_ips.Properties.RowNames

for i=1:size(achiasma,1)
    name=achiasma(i)
    namestr=name{1}

    %no_er_achiasma.Properties.VariableNames(i)=name
    %no_er_achiasma.(namestr)=no_er.(namestr)
    
    con_ips_achiasma.Properties.VariableNames(i)=name
    con_ips_achiasma.(namestr)=con_ips.(namestr)  
end


no_er_albinism=array2table(zeros(size(no_er,1),size(albinism,1)))
con_ips_albinism=array2table(zeros(size(con_ips,1),size(albinism,1)))

no_er_albinism.Properties.RowNames=no_er.Properties.RowNames
con_ips_albinism.Properties.RowNames=con_ips.Properties.RowNames

for i=1:size(albinism,1)
    name=albinism(i)
    namestr=name{1}

    %no_er_albinism.Properties.VariableNames(i)=name
    %no_er_albinism.(namestr)=no_er.(namestr)
    
    con_ips_albinism.Properties.VariableNames(i)=name
    con_ips_albinism.(namestr)=con_ips.(namestr)  
end

% Plot each group separately

contra_control
ipsi_control

contra_achiasma
ipsi_achiasma

contra_albinism
ipsi_albinism

%


contra = table2array(con_ips_controls({'Left to Right'},:)) + table2array(con_ips_controls({'Right to Left'},:));
ipsi = table2array(con_ips_controls({'Left to Left'},:)) + table2array(con_ips_controls({'Right to Right'},:));
names = con_ips_controls.Properties.VariableNames;

fig=scatter(ipsi, contra, 120, 'ko')
set(fig(1),'MarkerFaceColor','r')
hold on

title('Resulting streamlines count for LiFE filtering in optic chiasm')
xlabel('Number of Ipsilateral Streamlines', 'Fontsize', 30)
ylabel('Number of Contralateral Streamlines', 'Fontsize', 30) 

axis([0 100 0 100])
plot(linspace(0,100), linspace(0,100),'k-.')
legend('Controls', '50:50 ratio','Location','northwest')
text(ipsi +2, contra +1, names, 'Fontsize', 8)



% Plot remaining streamlines vs rmse


gscatter(mean_error, number_remaining, group,'brg','.xo')
%scatter(mean_error, number_remaining, 54, colors, 'filled')

title('Correlation between number of fibers remaining after filtering and root-mean-square error for optic chiasm data from 3 subjects')
xlabel('Root-mean-square Error')
ylabel('Number of fibres remaining after filtering')
text(mean_error+0.25, number_remaining, names)
axis([0 25 0 420])

legend('achiasma', 'albinism','controls')

% Plot contralateral streamlines vs ipsilateral with fitted line from literature

% Define seed features of plot
text(ipsi + 2, contra + 2, names)
axis([0 240 0 240])

% Define parameters of our control curve

cross= 555850
cross_er = 31400
uncross= 504250
uncross_er = 16600

exp_param = cross/uncross
exp_param_er = exp_param * sqrt((cross_er/cross)^2 + (uncross_er/uncross)^2)

exp_param_plus=exp_param + exp_param_er
exp_param_minus = exp_param - exp_param_er

% Plot published values that we expect

plot(linspace(1,240),exp_param*linspace(1,240),'k.')
title('Streamlines in optic chiasm remaining after LiFE filtering')
xlabel('Ipsilateral Streamlines')
ylabel('Contralateral Streamlines')
legend('Expected ratio','Location','northwest')

hold on
plot(linspace(1,240),exp_param_plus*linspace(1,240),'k-.')
plot(linspace(1,240),exp_param_minus*linspace(1,240),'k-.')
legend('Expected ratio', 'Expected ratio uncertainty','Location','northwest')

% Plot data for controls
ipsi_control=zeros(1,size(controls,1))
contra_control

for i=1:size(con_ips,2)
    if sum(strcmp(achiasma,names(i)))
        colors(i,:)=[0 0 1];
        group(1,i)=1
    end
    if sum(strcmp(albinism,names(i)))
        colors(i,:)=[1 0 0];
        group(1,i)=2
    end
    if sum(strcmp(controls,names(i)))
        colors(i,:)=[0 1 0];
        group(1,i)=3
    end
end



% Plot data for weights
cd('/home/auguser2016/Projects/Results/Report_for_Pestilli')
load('weights.mat')

achiasma={['hw91'];['ps94']};
albinism={['fe21'];['ib57'];['kw99'];['rx88'];['sj22'];['tq63'];['uh47']};
controls={['ce04'];['la21'];['lw37'];['nb30'];['ow93'];['uf97'];['xn78'];['xs62'];['ta14']};

names=weights.Properties.VariableNames

for i=1:size(weights,2)
    if sum(strcmp(achiasma,names(i)))
        colors(i,:)=[0 0 1];
        group(1,i)=1
    end
    if sum(strcmp(albinism,names(i)))
        colors(i,:)=[1 0 0];
        group(1,i)=2
    end
    if sum(strcmp(controls,names(i)))
        colors(i,:)=[0 1 0];
        group(1,i)=3
    end
end

weights_contra = table2array(weights({'Weights left to right'},:)) + table2array(weights({'Weights right to left'},:))
weights_ipsi = table2array(weights({'Weights left to left'},:)) + table2array(weights({'Weights right to right'},:))


fig = gscatter(weights_ipsi,weights_contra, group, 'brg', '.+x', 20, 'filled') %54,colors(nonan_contra,:),'filled')
fig=gcf
hold on
plot(linspace(0,40), linspace(0,40),'k-')

axis([0 40 0 40])
pbaspect([1 1 1])

title('Sum of normalized weights for crossing and noncrossing streamlines', 'Fontsize', 28)
xlabel('Sum of normalized weights of noncrossing streamlines', 'Fontsize', 22)
ylabel('Sum of normalized weights of crossing streamlines', 'Fontsize', 22)
legend({'Achiasma', 'Albinism','Controls','50:50 ratio'},'Location','northwest', 'Fontsize', 18)

text(weights_ipsi + 0.25, weights_contra + 0.25, names, 'Fontsize', 12)

% Plot data for SNR

cd('/home/auguser2016/Projects/Results/Report_for_Pestilli')
load('snr.mat')

achiasma={['hw91'];['ps94']};
albinism={['fe21'];['ib57'];['kw99'];['rx88'];['sj22'];['tq63'];['uh47']};
controls={['ce04'];['la21'];['lw37'];['nb30'];['ow93'];['uf97'];['xn78'];['xs62'];['ta14']};

names=snr.Properties.VariableNames

for i=1:size(snr,2)
    if sum(strcmp(achiasma,names(i)))
        colors(i,:)=[0 0 1];
        group(1,i)=1
    end
    if sum(strcmp(albinism,names(i)))
        colors(i,:)=[1 0 0];
        group(1,i)=2
    end
    if sum(strcmp(controls,names(i)))
        colors(i,:)=[0 1 0];
        group(1,i)=3
    end
end

snr_tab=table2array(snr)

fig = gscatter(linspace(1,size(snr_tab,2),18), snr_tab, group, 'brg', '.+x', 20, 'filled')

fig=gcf
hold on

title('Mean signal-to-noise ratio in optic chiasm', 'Fontsize', 28)
xlabel('Participant', 'Fontsize', 22)
ylabel('Signal-to-noise ratio', 'Fontsize', 22)
legend({'Achiasma', 'Albinism','Controls','50:50 ratio'},'Location','northwest', 'Fontsize', 18)

axis([0 19 8 12])

text(linspace(1,size(snr_tab,2),18) + 0.25, snr_tab + 0.25, names, 'Fontsize', 12)


% Plot FA

cd('/home/auguser2016/Projects/Results/Report_for_Pestilli')
load('fa.mat')

achiasma={['hw91'];['ps94']};
albinism={['fe21'];['ib57'];['kw99'];['rx88'];['sj22'];['tq63'];['uh47']};
controls={['ce04'];['la21'];['lw37'];['nb30'];['ow93'];['uf97'];['xn78'];['xs62'];['ta14']};

names=snr.Properties.VariableNames

for i=1:size(fa,2)
    if sum(strcmp(achiasma,names(i)))
        colors(i,:)=[0 0 1];
        group(1,i)=1
    end
    if sum(strcmp(albinism,names(i)))
        colors(i,:)=[1 0 0];
        group(1,i)=2
    end
    if sum(strcmp(controls,names(i)))
        colors(i,:)=[0 1 0];
        group(1,i)=3
    end
end

fa=table2array(fa)

fig = gscatter(linspace(1,size(fa,2),18), fa, group, 'brg', '.+x', 20, 'filled')

fig=gcf
hold on

title('Mean fractional anisotropy in optic chiasm', 'Fontsize', 28)
xlabel('Participant', 'Fontsize', 22)
ylabel('Fractional anisotropy', 'Fontsize', 22)
legend({'Achiasma', 'Albinism','Controls','50:50 ratio'},'Location','northwest', 'Fontsize', 18)

text(linspace(1,size(fa,2),18) + 0.25, fa + 0.005, names, 'Fontsize', 12)

% Plot data for SNR in midline point

cd('/home/auguser2016/Projects/Results/Report_for_Pestilli')
load('snr_mid.mat')

achiasma={['hw91'];['ps94']};
albinism={['fe21'];['ib57'];['kw99'];['rx88'];['sj22'];['tq63'];['uh47']};
controls={['ce04'];['la21'];['lw37'];['nb30'];['ow93'];['uf97'];['xn78'];['xs62'];['ta14']};

names=snr_mid.Properties.VariableNames

for i=1:size(snr_mid,2)
    if sum(strcmp(achiasma,names(i)))
        colors(i,:)=[0 0 1];
        group(1,i)=1
    end
    if sum(strcmp(albinism,names(i)))
        colors(i,:)=[1 0 0];
        group(1,i)=2
    end
    if sum(strcmp(controls,names(i)))
        colors(i,:)=[0 1 0];
        group(1,i)=3
    end
end

snr_mid_tab=table2array(snr_mid)

fig = gscatter(linspace(1,size(snr_mid_tab,2),18), snr_mid_tab, group, 'brg', 'xxx', 20, 'filled')

fig=gcf
hold on

title('SNR in middle of the optic chiasm', 'Fontsize', 28)
xlabel('Ordering number of participant', 'Fontsize', 22)
ylabel('Approximate SNR in the middle optic chiasm', 'Fontsize', 22)
legend({'Achiasma', 'Albinism','Controls','50:50 ratio'},'Location','northwest', 'Fontsize', 18)

text(linspace(1,size(snr_mid_tab,2),18) + 0.25, snr_mid_tab + 0.25, names, 'Fontsize', 12)


% Plot FA in the middle of chiasm

cd('/home/auguser2016/Projects/Results/Report_for_Pestilli')
load('fa_mid.mat')

achiasma={['hw91'];['ps94']};
albinism={['fe21'];['ib57'];['kw99'];['rx88'];['sj22'];['tq63'];['uh47']};
controls={['ce04'];['la21'];['lw37'];['nb30'];['ow93'];['uf97'];['xn78'];['xs62'];['ta14']};

names=fa_mid.Properties.VariableNames

for i=1:size(fa_mid,2)
    if sum(strcmp(achiasma,names(i)))
        colors(i,:)=[0 0 1];
        group(1,i)=1
    end
    if sum(strcmp(albinism,names(i)))
        colors(i,:)=[1 0 0];
        group(1,i)=2
    end
    if sum(strcmp(controls,names(i)))
        colors(i,:)=[0 1 0];
        group(1,i)=3
    end
end

fa_mid=table2array(fa_mid)

fig = gscatter(linspace(1,size(fa_mid,2),18), fa_mid, group, 'brg', 'xxx', 14, 'filled')

fig=gcf
hold on

title('FA in middle of optic chiasm', 'Fontsize', 20)
xlabel('Ordering number of participant', 'Fontsize', 16)
ylabel('Approximate FA in the middle in optic chiasm', 'Fontsize', 16)
legend({'Achiasma', 'Albinism','Controls','50:50 ratio'},'Location','northwest', 'Fontsize', 14)

text(linspace(1,size(fa_mid,2),18) + 0.25, fa_mid + 0.005, names, 'Fontsize', 12)


% Plot eye vs eye

cd('/home/auguser2016/Projects/Results/Report_for_Pestilli')
load('conips2.mat')

achiasma={['hw91'];['ps94']};
albinism={['fe21'];['ib57'];['kw99'];['rx88'];['sj22'];['tq63'];['uh47']};
controls={['ce04'];['la21'];['lw37'];['nb30'];['ow93'];['uf97'];['xn78'];['xs62'];['ta14']};

names=con_ips.Properties.VariableNames

for i=1:size(con_ips,2)
    if sum(strcmp(achiasma,names(i)))
        colors(i,:)=[0 0 1];
        group(1,i)=1
    end
    if sum(strcmp(albinism,names(i)))
        colors(i,:)=[1 0 0];
        group(1,i)=2
    end
    if sum(strcmp(controls,names(i)))
        colors(i,:)=[0 1 0];
        group(1,i)=3
    end
end

%right_eye=[con_ips{1,:}; con_ips{2,:}]
%right_eye=array2table(right_eye)
%right_eye.Properties.VariableNames = {subjects(4:end).name};
%right_eye.Properties.RowNames = {'Non-crossing nerves from right eye','Crossing nerves from right eye'};
%right_eye

%left_eye=[con_ips{4,:}; con_ips{3,:}]
%left_eye=array2table(left_eye)
%left_eye.Properties.VariableNames = {subjects(4:end).name};
%left_eye.Properties.RowNames = {'Non-crossing nerves from right eye','Crossing nerves from right eye'};
%left_eye

fig = gscatter(con_ips{{'Right to Right'},:}, con_ips{{'Right to Left'},:}, group, 'brg', 'xxx', 14, 'filled')
fig=gcf
hold on

plot(linspace(1,150), linspace(1,150),'k-.')

title('Streamlines count in eye to eye comparison', 'Fontsize', 20)
xlabel('Non-crossing streamlines from eye', 'Fontsize', 16)
ylabel('Crossing streamlines from eye', 'Fontsize', 16)

text(con_ips{{'Right to Right'},:} + 0.5, con_ips{{'Right to Left'},:} + 0.5, names, 'Fontsize', 12)

gscatter(con_ips{{'Left to Left'},:}, con_ips{{'Left to Right'},:}, group, 'brg', 'ooo', 14, 'filled')
text(con_ips{{'Left to Left'},:} + 0.5, con_ips{{'Left to Right'},:} + 0.5, names, 'Fontsize', 12)

legend({'Achiasma, right eye', 'Albinism, right eye','Controls, right eye','50:50 ratio', 'Achiasma, left eye', 'Albinism, left eye','Controls, left eye'},'Location','northwest', 'Fontsize', 14)

% Plot ratio vs SNR
load('snr.mat')
load('conips1.mat')

fig = gscatter((con_ips{{'Left to Right'},:}+con_ips{{'Right to Left'},:})./(con_ips{{'Right to Right'},:}+con_ips{{'Left to Left'},:} +(con_ips{{'Left to Right'},:}+con_ips{{'Right to Left'},:}) ),snr{'SNR',:},  group, 'brg', '.+x', 20, 'filled')
fig=gcf
hold on

title('Signal-to-noise ratio vs ratio of crossing streamlines', 'Fontsize', 28)
ylabel('Signal-to-noise ratio', 'Fontsize', 22)
xlabel('Ratio of crossing streamlines', 'Fontsize', 22)

%text(snr{'SNR',:} + 0.05, (con_ips{{'Left to Right'},:}+con_ips{{'Right to Left'},:})./(con_ips{{'Right to Right'},:}+con_ips{{'Left to Left'},:} +(con_ips{{'Left to Right'},:}+con_ips{{'Right to Left'},:}) ) + 0.01, names, 'Fontsize', 18)
pbaspect([1 1 1])
axis([0 0.7 8 12])
legend({'Achiasma', 'Albinism','Controls'},'Location','northwest', 'Fontsize', 18)

%% Plot ratio vs FA

load('fa.mat')
load('conips1.mat')

fig = gscatter((con_ips{{'Left to Right'},:}+con_ips{{'Right to Left'},:})./(con_ips{{'Right to Right'},:}+con_ips{{'Left to Left'},:} +(con_ips{{'Left to Right'},:}+con_ips{{'Right to Left'},:}) ),fa{'FA',:}, group, 'brg', '.+x', 20, 'filled')
fig=gcf
hold on

title('Fractional anisotropy vs ratio of crossing streamlines', 'Fontsize', 28)
ylabel('Fractional anisotropy', 'Fontsize', 22)
xlabel('Ratio of crossing streamlines', 'Fontsize', 22)

text(fa{'FA',:} + 0.00125, (con_ips{{'Left to Right'},:}+con_ips{{'Right to Left'},:})./(con_ips{{'Right to Right'},:}+con_ips{{'Left to Left'},:} +(con_ips{{'Left to Right'},:}+con_ips{{'Right to Left'},:}) ) + 0.01, names, 'Fontsize', 12)
pbaspect([1 1 1])
legend({'Achiasma', 'Albinism','Controls'},'Location','northwest', 'Fontsize', 18)

%% Plot ratio vs SNR for mid part of chiasm
load('snr_mid.mat')
load('conips1.mat')

fig = gscatter(snr_mid{'SNR',:},(con_ips{{'Left to Right'},:}+con_ips{{'Right to Left'},:})./(con_ips{{'Right to Right'},:}+con_ips{{'Left to Left'},:} +(con_ips{{'Left to Right'},:}+con_ips{{'Right to Left'},:}) ), group, 'brg', 'xxx', 14, 'filled')
fig=gcf
hold on

title('Crossing/All streamlines ratio vs SNR for middle part of chiasm', 'Fontsize', 20)
xlabel('SNR', 'Fontsize', 16)
ylabel('Crossing/All streamlines ratio', 'Fontsize', 16)

text(snr_mid{'SNR',:} + 0.05, (con_ips{{'Left to Right'},:}+con_ips{{'Right to Left'},:})./(con_ips{{'Right to Right'},:}+con_ips{{'Left to Left'},:} +(con_ips{{'Left to Right'},:}+con_ips{{'Right to Left'},:}) ) + 0.01, names, 'Fontsize', 12)

legend({'Achiasma', 'Albinism','Controls'},'Location','northwest', 'Fontsize', 14)

% Plot ratio vs FA for middle part of chiasm

load('fa_mid.mat')
load('conips1.mat')

fig = gscatter(fa_mid{'FA',:},(con_ips{{'Left to Right'},:}+con_ips{{'Right to Left'},:})./(con_ips{{'Right to Right'},:}+con_ips{{'Left to Left'},:} +(con_ips{{'Left to Right'},:}+con_ips{{'Right to Left'},:}) ), group, 'brg', 'xxx', 14, 'filled')
fig=gcf
hold on

title('Crossing/All streamlines ratio vs FA for middle part of chiasm', 'Fontsize', 20)
xlabel('FA', 'Fontsize', 16)
ylabel('Crossing/All streamlines ratio', 'Fontsize', 16)

text(fa_mid{'FA',:} + 0.00125, (con_ips{{'Left to Right'},:}+con_ips{{'Right to Left'},:})./(con_ips{{'Right to Right'},:}+con_ips{{'Left to Left'},:} +(con_ips{{'Left to Right'},:}+con_ips{{'Right to Left'},:}) ) + 0.01, names, 'Fontsize', 12)

legend({'Achiasma', 'Albinism','Controls'},'Location','northwest', 'Fontsize', 14)



% Add additionally standard deviation obtained from XXX's paper

% Fit linear function, obtained goodness of fit and estimated crossing

% fit linear function to each group separately

%}