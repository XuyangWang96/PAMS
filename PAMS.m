%load signals retrieved fron 268 mask
discovery=load('signal\discovery.mat');
external1=load('signal\external1.mat');
external2=load('signal\external2.mat');
discovery_signal=discovery.Signal.data;
external1_signal=external1.Signal.data;
external2_signal=external2.Signal.data;
data_path='E:\Xinan_Data\validations.xlsx';

%internal validation on discovery cohort
age1=xlsread(data_path,'discovery','D:D');
scores1=xlsread(data_path,'discovery','H:H');
[~,fc_mat]=W_pearsonfc(discovery_signal);
% set a series of thresholds, the q range for binarization and thresh for
% feature selection
range=0.8:0.01:0.9;
thresh = 0.001:0.001:0.01;
all_posr=zeros(length(range),length(thresh));
all_posp=zeros(length(range),length(thresh));
all_negr=zeros(length(range),length(thresh));
all_negp=zeros(length(range),length(thresh));
i=0;

for q=0.8:0.01:0.9
    
    i=i+1;
    accordance=zeros(size(discovery_signal,2),size(discovery_signal,2),size(discovery_signal,3));
    discordance=zeros(size(discovery_signal,2),size(discovery_signal,2),size(discovery_signal,3));
  
parfor sub=1:size(discovery_signal,3)
    sub
   [accordance(:,:,sub),discordance(:,:,sub)]=W_activeFC(discovery_signal(:,:,sub),q);

end


for thresh_i=1:length(thresh)
    cov=age1;
    all_mats= discordance;
    all_behav =scores1;
    % ---------------------------------------

        disp(['threshold P=',num2str(thresh(thresh_i))]);
        no_sub = size(all_mats,3);
        no_node = size(all_mats,1);
        behav_pred_pos = zeros(no_sub,1);
        behav_pred_neg = zeros(no_sub,1);
         %k-fold
        kfolds=size(all_mats,3);
        ksample=floor(no_sub/kfolds);
        randinds=randperm(no_sub);
        all_binary_posmask=zeros(no_node,no_node,kfolds);
        all_binary_negmask=zeros(no_node,no_node,kfolds);
        for leftout = 1:kfolds
        fprintf('\n Leaving out # %6.3f',leftout);
        if kfolds==no_sub%leave-one-out
            testinds=randinds(leftout);
            traininds=setdiff(randinds,testinds);
        else
            si=1+ksample*(leftout-1);
            fi=si+ksample-1;
            testinds=randinds(si:fi);
            traininds=setdiff(randinds,testinds);
        end

        % leave out subject from matrices and behavior
        train_mats = all_mats(:,:,traininds);
        train_vcts = reshape(train_mats,[],size(train_mats,3));
        train_behav = all_behav(traininds);
        train_cov=cov(traininds,:);
        [r_mat,p_mat]=partialcorr(train_vcts',train_behav,train_cov,'type','pearson');
     
        r_mat = reshape(r_mat,no_node,no_node);
        p_mat = reshape(p_mat,no_node,no_node);
        r_mat(isnan(r_mat))=0;
        p_mat(isnan(p_mat))=0;
       indPos = find(squareform(r_mat)>0 & squareform(p_mat) < thresh(thresh_i));
       indNeg = find(squareform(r_mat)<0 & squareform(p_mat) < thresh(thresh_i));

        binary_pos_mask = zeros(no_node,no_node);
        binary_neg_mask = zeros(no_node,no_node);

    
    binary_pos_edges = find(r_mat > 0 & p_mat < thresh(thresh_i));
    binary_neg_edges = find(r_mat < 0 & p_mat < thresh(thresh_i));
    
    binary_pos_mask(binary_pos_edges) = 1;
    binary_neg_mask(binary_neg_edges) = 1;
    
    all_binary_posmask(:,:,leftout)=binary_pos_mask;
    all_binary_negmask(:,:,leftout)=binary_neg_mask;
   
   
    train_sumpos = zeros(no_sub-ksample,1);
    train_sumneg = zeros(no_sub-ksample,1);
    train_posedges=zeros(no_sub-ksample,length(indPos));
    train_negedges=zeros(no_sub-ksample,length(indNeg));
    
    
    for ss = 1:size(train_posedges,1)
        submat=squareform(train_mats(:,:,ss));
        train_posedges(ss,:) = submat(indPos);
        train_negedges(ss,:) = submat(indNeg);
    end
    
    % build model on training set
    
   try
        [~,~,~,~,betaPLS_pos] = plsregress(train_posedges,train_behav,1);
    catch
        betaPLS_pos=zeros(size(train_posedges,2)+1,1)*nan;
    end
    try
        [~,~,~,~,betaPLS_neg] = plsregress(train_negedges,train_behav,1);
    catch
        betaPLS_neg=zeros(size(train_negedges,2)+1,1)*nan;
    end
   
    % run model on test set
    
    test_mat = squareform(all_mats(:,:,testinds));
    test_posedges=test_mat(indPos);
    test_negedges=test_mat(indNeg);
    
    behav_pred_pos(testinds) = [1 test_posedges]*betaPLS_pos;
    behav_pred_neg(testinds)=  [1 test_negedges]*betaPLS_neg;
   
    end

[R_pos, P_pos] = corr(behav_pred_pos,all_behav,'type','pearson');
[R_neg, P_neg] = corr(behav_pred_neg,all_behav,'type','pearson');
all_posr(i,thresh_i)=R_pos;
all_posp(i,thresh_i)=P_pos;
all_negr(i,thresh_i)=R_neg;
all_negp(i,thresh_i)=P_neg;
    end
end
%permutation test
no_iteration=10000;
discorPerm_r=zeros(no_iteration+1,1);
accorPerm_r=zeros(no_iteration+1,1);
pearPerm_r=zeros(no_iteration+1,1);
accorPerm_r(1)=R_accordance;% the true value between predicted scores and observed scores
pearPerm_r(1)=R_pearson;
discorPerm_r(1)=R_discordance;

for perm=2:no_iteration+1
    fprintf('\n Performing iteratioin %d ',perm-1);
    rng shuffle
    new_behav=all_behav(randperm(no_sub));
    [discorPerm_r(perm,1),~]=predict_behavior(discordance,new_behav,0.001,age1);
    [~,accorPerm_r(perm,1)]=predict_behavior(accordance,new_behav,0.001,age1);
    [~,pearPerm_r(perm,1)]=predict_behavior(fc_mat,new_behav,0.001,age1);
end

p_perm=(length(find(discorPerm_r>=discorPerm_r(1)))-1)/no_iteration
p_perm=(length(find(accorPerm_r>=accorPerm_r(1)))-1)/no_iteration
p_perm=(length(find(pearPerm_r>=pearPerm_r(1)))-1)/no_iteration

% make the unified predictive mask
u_negmask=zeros(length(all_binary_negmask),length(all_binary_negmask));
for node=1:size(all_binary_negmask,3)
    u_negmask=all_binary_negmask(:,:,node)+ u_negmask;
end
u_negmask(find(u_negmask<(size(all_binary_negmask,3)*0.95)))=0;
u_negmask(find(u_negmask>=(size(all_binary_negmask,3)*0.95)))=1;
u_posmask=zeros(length(all_binary_posmask),length(all_binary_posmask));
for node=1:size(all_binary_posmask,3)
    u_posmask=all_binary_posmask(:,:,node)+ u_posmask;
end
u_posmask(find(u_posmask<size(all_binary_posmask,3)*0.95))=0;
u_posmask(find(u_posmask>=size(all_binary_posmask,3)*0.95))=1;

%training the candidate predictive models on the whole discovery cohort

for sub=1:size(all_mats,3)
    all_edges(sub,:)=squareform(discordance(:,:,sub));
end
mask_disc=find(squareform(discorPos));
[~,~,~,~,beta_disc]=plsregress(all_edges(:,mask_disc),all_behav,1);

for sub=1:size(all_mats,3)
    all_edges(sub,:)=squareform(accordance(:,:,sub));
end
mask_acor=find(squareform(accorNeg));
[~,~,~,~,beta_acor]=plsregress(all_edges(:,mask_acor),all_behav,1);

for sub=1:size(all_mats,3)
    all_edges(sub,:)=squareform(fc_mat(:,:,sub));
end
mask_fc=find(squareform(pearsonNeg));
[~,~,~,~,beta_fc]=plsregress(all_edges(:,mask_fc),all_behav,1);


