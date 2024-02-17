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
        ind_test=setdiff(1:N,cv{i});
        ind_train=cv{i};
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

%Perform one-tailed, one-sample t-test
[~,p_pFC,~,stats_pFC]=ttest(r_intra_pFC-r_inter_pFC,0,'tail','right'); 
[~,p_SC,~,stats_SC]=ttest(r_intra_SC-r_inter_SC,0,'tail','right'); 
[~,p_bm,~,stats_bm]=ttest(r_intra_bm-r_inter_bm,0,'tail','right'); 
[~,p_null,~,stats_null]=ttest(r_intra_null-r_inter_null,0,'tail','right');
[~,p_mean,~,stats_mean]=ttest(r_intra_mean-r_inter_mean,0,'tail','right'); 
    
%Repeat this test using permutation
Perms=5000; null_pFC=zeros(Perms,1); null_SC=zeros(Perms,1); null_bm=zeros(Perms,1); 
null_null=zeros(Perms,1); null_mean=zeros(Perms,1); 
for i=1:Perms
    [~,~,~,tmp_pFC]=ttest((r_intra_pFC-r_inter_pFC).*sign(rand(N,1)-0.5)); 
    [~,~,~,tmp_SC]=ttest((r_intra_SC-r_inter_SC).*sign(rand(N,1)-0.5));
    [~,~,~,tmp_bm]=ttest((r_intra_bm-r_inter_bm).*sign(rand(N,1)-0.5));
    [~,~,~,tmp_null]=ttest((r_intra_null-r_inter_null).*sign(rand(N,1)-0.5));
    [~,~,~,tmp_mean]=ttest((r_intra_mean-r_inter_mean).*sign(rand(N,1)-0.5));
    null_pFC(i)=tmp_pFC.tstat; 
    null_SC(i)=tmp_SC.tstat; 
    null_bm(i)=tmp_bm.tstat; 
    null_null(i)=tmp_null.tstat; 
    null_mean(i)=tmp_mean.tstat; 
end
p_pFC_perm=sum(null_pFC>stats_pFC.tstat)/Perms; 
p_SC_perm=sum(null_SC>stats_SC.tstat)/Perms; 
p_bm_perm=sum(null_bm>stats_bm.tstat)/Perms; 
p_null_perm=sum(null_null>stats_null.tstat)/Perms; 
p_mean_perm=sum(null_mean>stats_mean.tstat)/Perms; 

rval=corr(r_intra_pFC,r_inter_pFC);
fprintf('pFC-eFC: Intra: r=%0.4f, Inter: r=%0.4f, p=%0.4f, perm=%0.4f, Cohen d=%0.2f\n',...
       mean(r_intra_pFC),mean(r_inter_pFC),p_pFC,p_pFC_perm,stats_pFC.tstat/(sqrt(stats_pFC.df) *sqrt(1-rval^2)) );

fprintf('Percent improvement=%0.2f%%\n', (mean(r_intra_pFC)-mean(r_inter_pFC))/mean(r_intra_pFC)*100 ); 
%Percent improvement measure is not meaningful when standardization used

rval=corr(r_intra_SC,r_inter_SC); 
fprintf('SC-eFC: Intra: r=%0.4f, Inter: r=%0.4f, p=%0.4f, perm=%0.4f, Cohen d=%0.2f\n',...
        mean(r_intra_SC),mean(r_inter_SC),p_SC,p_SC_perm,stats_SC.tstat/(sqrt(stats_SC.df) *sqrt(1-rval^2))  ) ; 

rval=corr(r_intra_bm,r_inter_bm); 
fprintf('bm-eFC: Intra: r=%0.4f, Inter: r=%0.4f, p=%0.4f, perm=%0.4f, Cohen d=%0.2f\n',...
        mean(r_intra_bm),mean(r_inter_bm),p_bm,p_bm_perm,stats_bm.tstat/(sqrt(stats_bm.df) *sqrt(1-rval^2))  ) ; 

rval=corr(r_intra_null,r_inter_null); 
fprintf('null-eFC: Intra: r=%0.4f, Inter: r=%0.4f, p=%0.4f, perm=%0.4f, Cohen d=%0.2f\n',...
        mean(r_intra_null),mean(r_inter_null),p_null,p_null_perm,stats_null.tstat/(sqrt(stats_null.df) *sqrt(1-rval^2))  ) ; 

rval=corr(r_intra_mean,r_inter_mean); 
fprintf('mean_eFC-eFC: Intra: r=%0.4f, Inter: r=%0.4f, p=%0.4f, perm=%0.4f, Cohen d=%0.2f\n',...
        mean(r_intra_mean),mean(r_inter_mean),p_mean,p_mean_perm,stats_mean.tstat/(sqrt(stats_mean.df) *sqrt(1-rval^2))  ) ; 

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
ha.YLim=[-0.05,0.05];  
ha.XTick=[1,2,3]; 
ha.XTickLabel={'pFC','SC','Null'};
ha.YLabel.String='Correlation';  
ha.FontSize=16;