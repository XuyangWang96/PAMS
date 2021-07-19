function [r_pos,r_neg]=predict_behavior(all_mats,all_behav,thresh,cov)
% this function is used in permutation test, repeat the same procedure in
% PAMS.m
disp(['threshold P=',num2str(thresh)]);
no_sub = size(all_mats,3);
no_node = size(all_mats,1);

behav_pred_pos = zeros(no_sub,1);
behav_pred_neg = zeros(no_sub,1);


parfor leftout = 1:no_sub
    fprintf('\n Leaving out subj # %6.3f',leftout);
    
    
    train_mats = all_mats;
    train_mats(:,:,leftout) = [];
    train_vcts = reshape(train_mats,[],size(train_mats,3));
    train_behav = all_behav;
    train_behav(leftout) = [];
    train_cov=cov;
    train_cov(leftout,:) = [];
    
    
    
    [r_mat,p_mat]=partialcorr(train_vcts',train_behav,train_cov,'type','pearson');
    
    r_mat = reshape(r_mat,no_node,no_node);
    p_mat = reshape(p_mat,no_node,no_node);
    r_mat(isnan(r_mat))=0;
    p_mat(isnan(p_mat))=0;
    indPos = find(squareform(r_mat)>0 & squareform(p_mat) < thresh);
    indNeg = find(squareform(r_mat)<0 & squareform(p_mat) < thresh);
    
    
    train_posedges=zeros(no_sub-1,length(indPos));
    train_negedges=zeros(no_sub-1,length(indNeg));
    
    for ss = 1:size(train_posedges,1)
        submat=squareform(train_mats(:,:,ss));
        train_posedges(ss,:) = submat(indPos);
        train_negedges(ss,:) = submat(indNeg);
    end
    
    % build model on TRAIN subs
   
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
    % run model on TEST sub
    test_mat = squareform(all_mats(:,:,leftout));
    test_posedges=test_mat(indPos);
    test_negedges=test_mat(indNeg);
    behav_pred_pos(leftout) = [1 test_posedges]*betaPLS_pos;
    behav_pred_neg(leftout)=  [1 test_negedges]*betaPLS_neg;
    
end
   
    r_pos=corr(all_behav,behav_pred_pos);
    r_neg=corr(all_behav,behav_pred_neg);
    
end