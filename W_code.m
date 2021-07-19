function [codeUp,codeLow]=W_code(s,q)
%s:the normalized time series
UP=icdf('Normal',q,0,1);
Down=icdf('Normal',1-q,0,1);
codeUp=zeros(size(s));
codeLow=zeros(size(s));
    for i=1:length(s)
        if s(i)>=UP
            codeUp(i)=1;
        elseif s(i)<=Down
            codeLow(i)=1;
        end
    end
end