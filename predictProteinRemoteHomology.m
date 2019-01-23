
function [family, prob] = predictProteinRemoteHomology(queryProteinFile)
    load 3187sequences_Grey21PSSM.mat %psepssm
    load 3187pfams.mat
    load 3187sequencesFamily.mat
    load pfamA-Clans.mat
    
    [queryheads,~]= fastaread(queryProteinFile);
    [trainheads,~] = fastaread('3187seqs.fasta');
        
    tlen = length(trainheads);
    
    if iscell(queryheads)
        qlen = length(queryheads);
        for i = 1 : qlen
            temp = split(queryheads{i});
            queryheads{i} = temp(1);
        end
    else
        qlen = 1;
        temp = split(queryheads);
        queryheads = temp(1);
    end
    
    for i = 1 : tlen
        temp = split(trainheads{i});
        trainheads{i} = temp(1);
    end
    
    family = cell(1,qlen);
    prob = cell(1,qlen);
    
    % 1. grey-PSSM
    % 1.1 get pssm of query proteins
    tp=blastpssm(queryProteinFile,'swissprot');
    if iscell(tp)
        p = tp;
    else
        p{1} = tp;
    end

    % 1.2 grey(2,1) PSSM
    queryPsepssm=greyPsePssm(p,2);
    % 2. functional domain set index
    dist_DSI = zeros(qlen,tlen);
    queryFams = buildFunctionDomainVect(queryProteinFile);
    for i = 1 : qlen
        for j = 1 : tlen
            dist_DSI = pfamcmp(queryFams{i},pfams{j},clanDict);
        end
    end
    
    % 3. merge distance of grey incidence alansys and domain set index
    for i = 1 : qlen
        dist_GIA = GreyIncidenceDegree(queryPsepssm(i,:),psepssm);
        dist = dist_GIA + dist_DSI(i,:)
        [B,I] = sort(dist,'descend');
        family{i} = familyId(I(1));
        prob{i} = B(I(1));
    end
   
