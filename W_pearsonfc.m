function [VecFC,FC,nums]=W_pearsonfc(data)
%Input: data: (timepoints,rois,subjects)

%Output: FC: (rois,rois,subjects)
%        nums: the number of connectivity betweeen each pair of nodes
%        VecFC:vectorize the upper triangle

FC_data=zeros(size(data,2),size(data,2),size(data,3));
FC=zeros(size(data,2),size(data,2),size(data,3));
for i=1:size(data,3)
    FC_data(:,:,i)=corr(data(:,:,i));
    FC(:,:,i)=0.5*log((1+FC_data(:,:,i))./(1-FC_data(:,:,i)));
end

%%矩阵的上三角拉直

nums=nchoosek(size(FC,2),2);
VecFC=zeros(nums,size(FC,3));
for v=1:size(data,3)
VecFC(:,v)=W_vec(FC(:,:,v))';
end
end