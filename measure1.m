%Measure 1
%Statistical testing of intra- and inter-individual correlation in eFC and
%pFC. Hypothesis is that intra-indivdal correlation is greater than
%inter-individual correlation. This suggests that pFC captures individuals
%effects. 
clear all
close all

addpath .\Violinplot-Matlab-master
addpath .\cbrewer
ct=cbrewer('qual', 'Set1', 8);

 
%Standardization of SC, pFC and eFC across subjects?
Standardize=0; 
%0: no standardization
%1: demean only
%2: zscore

%Cross-validation type 
CV=1;
%0: standardize independent of folds
%1: standardize while respecting folds

%Load SC, pFC and eFC matrices from Sarwar et al
%All stored as vectorized upper triangular elements 
load full_data_with_pFCs.mat
efc=FC_emp;    %empirical FC
sc=SC;         %structural connectome
nnfc=FC_NN;    %predicted FC
bmfc=FC_bm;    %biophysical model FC (not used here)
N=size(efc,1); %number of subjects
J=size(efc,2); %number of connections
sc_old=sc; 

%Cross-validation partition for z-scoring
%Replicate 10-fold structure from Sarwar et al 
sz=floor(N/10);
res=N-10*sz;
for i=1:10
    cv{i}=(i-1)*sz+1:i*sz;
end
cv{10}=[cv{10},cv{10}(end)+1:cv{10}(end)+res];

if CV==0 %no cross-validation
    if Standardize==2
        nnfc=zscore(nnfc); 
        bmfc=zscore(bmfc); 
        mean_efc=zeros(N,J);  %benchmark prediction
    elseif Standardize==1
        nnfc=nnfc-repmat(mean(nnfc),N,1);
        bmfc=bmfc-repmat(mean(bmfc),N,1);
        mean_efc=zeros(N,J); 
    elseif Standardize==0
        mean_efc=repmat(mean(efc),N,1); 
    end
else   %cross-validation based standardization
    mean_efc=zeros(N,J); 
    for i=1:length(cv)
        ind_train=setdiff(1:N,cv{i}); %swapped ind_test and ind_train (bug)
        ind_test=cv{i};
        mu_efc=mean(efc(ind_train,:)); std_efc=std(efc(ind_train,:));
        mu_nnfc=mean(nnfc(ind_train,:)); std_nnfc=std(nnfc(ind_train,:));
        mu_bmfc=mean(bmfc(ind_train,:)); std_bmfc=std(bmfc(ind_train,:));
        if Standardize==1 || Standardize==2
            %remove mean
            nnfc(ind_test,:)=nnfc(ind_test,:)-repmat(mu_nnfc,length(ind_test),1);
            bmfc(ind_test,:)=bmfc(ind_test,:)-repmat(mu_bmfc,length(ind_test),1);
            mean_efc(ind_test,:)=repmat(mean(efc(ind_train,:))-mu_efc,length(ind_test),1);
        end
        if Standardize==2
            nnfc(ind_test,:)=nnfc(ind_test,:)./repmat(std_nnfc,length(ind_test),1);
            nnfc(ind_test,std_nnfc==0);
            bmfc(ind_test,:)=bmfc(ind_test,:)./repmat(std_bmfc,length(ind_test),1);
            bmfc(ind_test,std_bmfc==0);
        end
        if Standardize==0
            mean_efc(ind_test,:)=repmat(mean(efc(ind_train,:)),length(ind_test),1);
        end
    end
end

%No cross-validation structure for eFC and SC
if Standardize==2
    sc=zscore(sc); efc=zscore(efc);
elseif Standardize==1
    sc=sc-repmat(mean(sc),N,1); efc=efc-repmat(mean(efc),N,1);
end

%Compute r_inter and r_intra, N x N matrices
%eFC with pFC 
r_pFC=corr(efc',nnfc');

%eFC with biophysical model (not used here)
r_bm=corr(efc',bmfc');

%eFC with SC 
r_SC=corr(efc',sc');

%eFC with mean efc  (benchmark 1)
r_mean=corr(efc',mean_efc'); 

%eFC with null (benchmark 2)
r_null=corr(efc',(repmat(std(efc),N,1).*randn(N,J)+mean_efc)'); 

%Extract intra-individual correlation
r_intra_pFC=diag(r_pFC); 
r_intra_SC=diag(r_SC); 
r_intra_bm=diag(r_bm); 
r_intra_null=diag(r_null); 
r_intra_mean=diag(r_mean); 

%Extract inter-individual correlation for each individual
r_inter_pFC=zeros(N,1);
r_inter_SC=zeros(N,1); 
r_inter_bm=zeros(N,1); 
r_inter_null=zeros(N,1);
r_inter_mean=zeros(N,1); 
for i=1:N
    r_inter_pFC(i)=mean(r_pFC(i,[1:i-1,i+1:N]));
    r_inter_SC(i)=mean(r_SC(i,[1:i-1,i+1:N]));
    r_inter_bm(i)=mean(r_bm(i,[1:i-1,i+1:N]));
    r_inter_null(i)=mean(r_null(i,[1:i-1,i+1:N]));
    r_inter_mean(i)=mean(r_mean(i,[1:i-1,i+1:N]));
end

%%
%TEST 1
%Combine predictions across all folds and conduct t-test
%This ignores dependence introduced through cross-validation 
%%
%Perform one-tailed, one-sample t-test
[~,p_pFC,~,stats_pFC]=ttest(r_intra_pFC-r_inter_pFC,0,'tail','right'); 
[~,p_SC,~,stats_SC]=ttest(r_intra_SC-r_inter_SC,0,'tail','right'); 
[~,p_bm,~,stats_bm]=ttest(r_intra_bm-r_inter_bm,0,'tail','right'); 
[~,p_null,~,stats_null]=ttest(r_intra_null-r_inter_null,0,'tail','right');
[~,p_mean,~,stats_mean]=ttest(r_intra_mean-r_inter_mean,0,'tail','right'); 
%%

%%
%TEST 2
%Prediction from each fold evaluated independently 
%Comparison of pFC prediction to benchmark prediction 
%Adjusted t-statistic as described in Section 3.2 of Bouckaert & Frank
%This accounts for the dependence structure between folds

%Select benchmark type
Benchmark=1;
%1=mean efc 
%2=mean efc + noise 

%Difference in performance between pFC prediction and benchmark prediction
%(r_intra-r_inter)_{pFC} - (r_intra-r_inter)_benchmark
d=zeros(length(cv),1); 
for n=1:length(cv) %loop across cv folds
    Ncv=length(cv{n}); %number of samples in test fold
    tmp=r_pFC(cv{n},cv{n}); %subsample
    if Benchmark==1
        tmp_null=r_mean(cv{n},cv{n});
    elseif Benchmark==2
        tmp_null=r_null(cv{n},cv{n});
    end
    r_intra_pFC_adj=diag(tmp);
    r_intra_null_adj=diag(tmp_null);
    r_inter_pFC_adj=zeros(Ncv,1);
    r_inter_null_adj=zeros(Ncv,1);
    for i=1:Ncv
       r_inter_pFC_adj(i)=mean(tmp(i,[1:i-1,i+1:Ncv]));
       r_inter_null_adj(i)=mean(tmp_null(i,[1:i-1,i+1:Ncv]));
    end
    %Store difference in performance between  model for each fold
    dk{n}=(r_intra_pFC_adj-r_inter_pFC_adj) - (r_intra_null_adj-r_inter_null_adj); 
    %Compute central tendency for t-statistc
    d(n)=median((r_intra_pFC_adj-r_inter_pFC_adj)  -  (r_intra_null_adj-r_inter_null_adj));
    %Distribution is skewed, and thus median used as tendency measure 
end
mu=mean(d);
v=var(d,1);
n1=N-sz; %number of training samples
n2=sz;   %number of testing samples
t=mu/sqrt((1/length(cv)+(n2/n1))*v); %See Bouckaert & Frank
p_adj = tcdf(t,10*length(cv)-1,'upper'); 
fprintf('Adjusted (r_intra-r_inter)_pFC - (r_intra-r_inter)_bench: p=%0.4f\n',p_adj); 

%Compute bootrapped intervals for each fold
Boots=5000; 
crit_val=zeros(1,length(cv)); 
for n=1:length(cv)
    str=zeros(Boots,1); 
    for i=1:Boots
        ind_rand=randi([1,length(dk{n})],length(dk{n}),1); 
        str(i)=median(dk{n}(ind_rand)); 
    end
    str=sort(str);
    crit=floor(0.05*Boots);  
    crit_val(n)=str(crit); 
end
hf=figure; hf.Color='w'; hf.Position=[100,100,1000,300]; 
for i=1:length(cv)
    subplot(2,5,i); [f,xi] = ksdensity(dk{i},linspace(-0.1,0.1,100)); plot(xi,f,'LineWidth',1.5); hold on; 
    ha=gca; ha.YLim=[0,4.4]; ha.FontSize=12; ha.XGrid='on'; 
    plot([median(dk{i}),median(dk{i})],[0,4],'k-'); 
    plot([crit_val(i),crit_val(i)],[0,4],'r--'); 
    if crit_val(i)>0
        title(sprintf('Fold %d, p<0.05',i)); 
    else
        title(sprintf('Fold %d, n.s.',i)); 
    end
end

%%

%%  
%TEST 3
%Permutations constrained to CV folds
%This ensures any dependence introduced through CV is preserved in the null
%distribution 
Perms=5000; t_null=zeros(Perms,1); 
for i=1:Perms
    r_pFCx=zeros(N,N); 
    for j=1:length(cv)
        tmp=r_pFC(:,cv{j}); 
        tmp=tmp(:,randperm(length(cv{j}))); %permute columns within fold
        r_pFCx(:,cv{j})=tmp; %insert into full similarity matrix
    end
    r_intra_pFCx=diag(r_pFCx); 
    r_inter_pFCx=zeros(N,1); 
    for j=1:N
        r_inter_pFCx(j)=mean(r_pFCx(j,[1:j-1,j+1:N]));
    end
    [~,~,~,stats_pFCx]=ttest(r_intra_pFCx-r_inter_pFCx,0,'tail','right');
    t_null(i)=stats_pFCx.tstat; 
end
p_pFC_perm=sum(stats_pFC.tstat<t_null)/Perms; 
hf=figure; hf.Color='w'; hf.Position=[100,100,400,300]; 
hist(t_null,100); h = findobj(gca,'Type','patch'); h.FaceColor = [1 0 0]; hold on; 
plot([stats_pFC.tstat,stats_pFC.tstat],[0,100],'k-','LineWidth',1.5);
ha=gca; ha.FontSize=16;
%%


rval=corr(r_intra_pFC,r_inter_pFC);
fprintf('pFC-eFC: Intra: r=%0.4f, Inter: r=%0.4f, p=%0.4f, perm=%0.4f, Cohen d=%0.2f\n',...
       mean(r_intra_pFC),mean(r_inter_pFC),p_pFC,p_pFC_perm,stats_pFC.tstat/(sqrt(stats_pFC.df) *sqrt(1-rval^2)) );

fprintf('Percent improvement=%0.2f%%\n', mean((r_intra_pFC-r_intra_null)./r_intra_pFC)*100 ); 
%Percent improvement measure is not meaningful when standardization used

rval=corr(r_intra_SC,r_inter_SC); 
fprintf('SC-eFC: Intra: r=%0.4f, Inter: r=%0.4f, p=%0.4f, Cohen d=%0.2f\n',...
        mean(r_intra_SC),mean(r_inter_SC),p_SC,stats_SC.tstat/(sqrt(stats_SC.df) *sqrt(1-rval^2))  ) ; 

rval=corr(r_intra_bm,r_inter_bm); 
fprintf('bm-eFC: Intra: r=%0.4f, Inter: r=%0.4f, p=%0.4f, Cohen d=%0.2f\n',...
        mean(r_intra_bm),mean(r_inter_bm),p_bm, stats_bm.tstat/(sqrt(stats_bm.df) *sqrt(1-rval^2))  ) ; 

rval=corr(r_intra_null,r_inter_null); 
fprintf('null-eFC: Intra: r=%0.4f, Inter: r=%0.4f, p=%0.4f, Cohen d=%0.2f\n',...
        mean(r_intra_null),mean(r_inter_null),p_null,stats_null.tstat/(sqrt(stats_null.df) *sqrt(1-rval^2))  ) ; 

rval=corr(r_intra_mean,r_inter_mean); 
fprintf('mean_eFC-eFC: Intra: r=%0.4f, Inter: r=%0.4f, p=%0.4f, Cohen d=%0.2f\n',...
        mean(r_intra_mean),mean(r_inter_mean),p_mean,stats_mean.tstat/(sqrt(stats_mean.df) *sqrt(1-rval^2))  ) ; 

%Generate figure
hf=figure; hf.Color='w'; hf.Position=[100,100,600,500]; 
vs{1} = Violin({r_intra_pFC},0.95,...
    'QuartileStyle','shadow',... % boxplot, none
    'DataStyle', 'none',... % scatter, histogram
    "HalfViolin","left","ViolinColor",{ct(2,:)},...
    'ShowMean',true);
vs{2} = Violin({r_inter_pFC},1.05,...
    'QuartileStyle','shadow',... % boxplot, none
    'DataStyle', 'none',... % scatter, histogram
    "HalfViolin","right","ViolinColor",{ct(1,:)},...
    'ShowMean',true);
vs{3} = Violin({r_intra_SC},1.95,...
    'QuartileStyle','shadow',... % boxplot, none
    'DataStyle', 'none',... % scatter, histogram
    'ShowMean',true,...
    "HalfViolin","left","ViolinColor",{ct(2,:)});
vs{4} = Violin({r_inter_SC},2.05,...
    'QuartileStyle','shadow',... % boxplot, none
    'DataStyle', 'none',... % scatter, histogram
    'ShowMean',true,...
    "HalfViolin","right","ViolinColor",{ct(1,:)});
vs{5} = Violin({r_intra_null},2.95,...
    'QuartileStyle','shadow',... % boxplot, none
    'DataStyle', 'none',... % scatter, histogram
    'ShowMean',true,...
    "HalfViolin","left","ViolinColor",{ct(2,:)});
vs{6} = Violin({r_inter_null},3.05,...
    'QuartileStyle','shadow',... % boxplot, none
    'DataStyle', 'none',... % scatter, histogram
    'ShowMean',true,...
    "HalfViolin","right","ViolinColor",{ct(1,:)});

ha=gca;
%ha.YLim=[-0.05,0.05];  
ha.XTick=[1,2,3]; 
ha.XTickLabel={'pFC','SC','Null'};
ha.YLabel.String='Correlation';  
ha.FontSize=16;
