
function [family, prob] = predictProteinRemoteHomology(queryProteinFile)
    load 3187sequences_Grey21PSSM.mat %psepssm
    load 3187FunctionDomainMatch_0.1.mat
    load 3187sequencesFamily.mat

    x=[0.3 0.02 0.28 0.4];
    
    [queryheads,~]= fastaread(queryProteinFile);
    [trainheads,~] = fastaread('7329seqs.fasta');
    
    
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
    
    % psi-blast bit score
    cmd = ['psiblast -db 7329seqsdb -query ' queryProteinFile ' -out blastpairfile -outfmt 6'];
    system(cmd, '-echo')
    [bitscore,~] = pairVals('blastpairfile', queryheads, trainheads);
    
    % hmmer score
    cmd = ['jackhmmer  --tblout hmmerscoreout.txt ' queryProteinFile ' 7329seqs.fasta'];
    system(cmd, '-echo');
    hmmscore = ones(qlen,tlen);
    fid = fopen('hmmerscoreout.txt','r');
    tline = fgetl(fid);
    while (ischar(tline) && ~isempty(tline))
        if( startsWith(tline, '#'))
            tline = fgetl(fid);
            continue;
        else
            A = split(tline);
            targetProtein = A(1);
            queryProtein = A(3);
            i = findStringInCell(queryheads, queryProtein);
            j = findStringInCell(trainheads, targetProtein);
            hmmscore(i,j) = str2double(A(6));
            tline = fgetl(fid);
        end
    end
    fclose(fid);
    
    % get pssm of query proteins
    tp=blastpssm(queryProteinFile,'swissprot');
    if iscell(tp)
        p = tp;
    else
        p{1} = tp;
    end
    % glcm
    queryGlcm = getPSSM_GLCM(p);
    % grey(2,1) PSSM
    queryPsepssm=greyPsePssm(p,2);
    
    for i = 1 : qlen
        r = GreyIncidenceDegree(queryPsepssm(i,:),psepssm);
        d = GreyIncidenceDegree(queryGlcm(i,:), X);
        b = mapminmax(bitscore(i,:),0,1);
        h = mapminmax(hmmscore(i,:),0,1);
        dis=x(1)*r + x(2)*b + x(3)*d + x(4)*h;
        [B,I] = sort(dis,'descend');
        family{i} = familyId(I(1));
        prob{i} = B(I(1));
    end
   
