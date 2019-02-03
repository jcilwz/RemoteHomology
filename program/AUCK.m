function auc = AUCK(label,preval,N,direction)
%preval ÉıĞòÅÅÁĞ
    if nargin < 4
        N = 0;
        direction='ascend';
    end
    [~,I] = sort(preval,direction);
    
    tp = 0;
    fp = 0;
    auc = 0;
    S = label(I);
    len = length(S);
    
    if N == 0
        k = 1;
        while(k<=len)
            if S(k) == 1
                tp = tp + 1;
            else
                fp = fp + 1;
                auc = auc + tp;
            end
            k=k+1;
        end
    else
        k=1;
        while( k <= len)
            if S(k) == 1
                tp = tp + 1;
            else
                fp = fp + 1;
                auc = auc + tp;
            end
            k=k+1;
            if fp == N
                break;
            end
        end
    end
    
    %%
    if tp == 0
        auc = 0;
    else
        if fp == 0
            auc = 1;
        else
            auc = auc/(tp*fp);
        end
    end
end