%Measure 1
%Matching individuals based on eFC and pFC. Involves solution of a binary
%assignment problem using the Hungarian algorithm. Statistical testing of
%whether the number of individuals correctly assigned in significantly
%greater than two benchmark/null conditions.
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

vals=2:20;   %number of indidivudal comrpising a subset
Trials=5000; %number of random trails
%Matrices storing counts of correct assigments
cnt_pFC=zeros(Trials,length(vals));
cnt_SC=zeros(Trials,length(vals));
cnt_bm=zeros(Trials,length(vals)); 
cnt_null=zeros(Trials,length(vals)); 

%Select null type
NullType=1; 
%1: similarity matric replaced with a random matrix
%2: rows (or columns) of r_pFC permuted

frst=0; 
for t=1:Trials %average results over multiple trials
    ind_rand=randperm(N); %samples a random ordering of individuals
    if NullType==1
        r_null=randn(N,N); %generate a random similarity matching
    elseif NullType==2
        r_null=r_pFC; 
        r_null=r_null(:,randperm(N)); %permute columns (or rows) 
    end
    for i=1:length(vals)
        ind=ind_rand(1:vals(i)); %Sample a subset of individuals
        %Trim similarity matrices to individuals comprising the subset
        r_pFC_new=r_pFC(ind,ind);  
        r_SC_new=r_SC(ind,ind); 
        r_bm_new=r_bm(ind,ind); 
        r_null_new=r_null(ind,ind); 
        %Hungarian algorithm
        %Note that -1 is needed because cost is minimized
        assign_pFC=munkres(-1*r_pFC_new); 
        assign_SC=munkres(-1*r_SC_new);
        assign_bm=munkres(-1*r_bm_new); 
        assign_null=munkres(-1*r_null_new); 
        %Compute percentage of individuals correctly matched
        cnt_pFC(t,i)=sum(diag(assign_pFC))/length(ind)*100; 
        cnt_SC(t,i)=sum(diag(assign_SC))/length(ind)*100; 
        cnt_bm(t,i)=sum(diag(assign_bm))/length(ind)*100; 
        cnt_null(t,i)=sum(diag(assign_null))/length(ind)*100; 
    end
    show_progress(t,Trials,frst); frst=1; 
end

%Standard errors
cnt_pFC_se=std(cnt_pFC)/sqrt(Trials); cnt_SC_se=std(cnt_SC)/sqrt(Trials); 
cnt_null_se=std(cnt_null)/sqrt(Trials); cnt_bm_se=std(cnt_bm)/sqrt(Trials);
%Means
cnt_pFC_mean=mean(cnt_pFC); cnt_SC_mean=mean(cnt_SC); cnt_null_mean=mean(cnt_null);
cnt_bm_mean=mean(cnt_bm); 

hf=figure; hf.Color='w'; hf.Position=[400,100,400,400]; 
mn=[mean(cnt_pFC_mean(1:length(vals))),mean(cnt_SC_mean(1:length(vals))),mean(cnt_null_mean(1:length(vals)))];
b=bar([(mn(1)-mn(3))/mn(1)*100,(mn(2)-mn(3))/mn(2)*100]); 
ha=gca; 
ha.XTickLabel={'pFC','SC'}; 
ha.FontSize=16; 


hf=figure; hf.Color='w'; hf.Position=[400,100,400,400]; 
plot(vals,cnt_pFC_mean,'r'); hold on; plot(vals,cnt_pFC_mean+1.96*cnt_pFC_se,'r--'); plot(vals,cnt_pFC_mean-1.96*cnt_pFC_se,'r--');
plot(vals,cnt_SC_mean,'b'); plot(vals,cnt_SC_mean+1.96*cnt_SC_se,'b--'); plot(vals,cnt_SC_mean-1.96*cnt_SC_se,'b--');
plot(vals,cnt_null_mean,'k'); plot(vals,cnt_null_mean+1.96*cnt_null_se,'k--'); plot(vals,cnt_null_mean-1.96*cnt_null_se,'k--');

[h,pp_pFC]=ttest2(cnt_null,cnt_pFC); 
for i=1:length(vals)
    if h(i)
        text(vals(i),2,'$$\ast$$','Interpreter','latex'); 
    end
end
