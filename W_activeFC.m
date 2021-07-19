function [accordance,discordance]=W_activeFC(data,q)
% data:times*rois
%q:quantile

rois=size(data,2);

dataNorm=(data-mean(data))./std(data);
accordance=zeros(rois,rois);
discordance=zeros(rois,rois);

for row=1:rois
    for col=1:rois
        [codeUp1,codeLow1]=W_code(dataNorm(:,row),q);
        [codeUp2,codeLow2]=W_code(dataNorm(:,col),q);
        sigma1=sqrt(dot(codeUp1,codeUp1)+dot(codeLow1,codeLow1));
        sigma2=sqrt(dot(codeUp2,codeUp2)+dot(codeLow2,codeLow2));
       
        if row < col
            accordance(row,col)=(dot(codeUp1,codeUp2)+dot(codeLow1,codeLow2))/(sigma1*sigma2);
            discordance(row,col)=(dot(codeUp1,codeLow2)+dot(codeLow1,codeUp2))/(sigma1*sigma2);
        end
       
    end
end
accordance=accordance+accordance';
discordance=discordance+discordance';
end



