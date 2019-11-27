load 7329sequencesGrey21PSSM.mat %psepssm
load 7329pfams.mat
load 7329sequencesFamily.mat
load pfamA-Clans.mat
queryProteinFile = 'H:\\independ_test.fa';
[queryheads,~]= fastaread(queryProteinFile);
%[trainheads,~] = fastaread('7329seqs.fasta');

%tlen = length(trainheads);
tlen = 7329;
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

%for i = 1 : tlen
%    temp = split(trainheads{i});
%    trainheads{i} = temp(1);
%end

%family = cell(1,qlen);
%prob = cell(1,qlen);

% 1. grey-PSSM
% 1.1 get pssm of query proteins
tp = genPSSMFromFile(queryProteinFile);
if iscell(tp)
    p = tp;
else
    p{1} = tp;
end

% 1.2 grey(2,1) PSSM
queryPsepssm=greyPsePssm(p,2);
% 2. functional domain set index
dist_DSI = zeros(qlen,tlen);
dist = zeros(qlen,tlen);
queryFams = buildFunctionDomainSet(queryProteinFile);
for i = 1 : qlen
    for j = 1 : tlen
        dist_DSI(i,j) = pfamcmp(queryFams{i},pfams{j},clanDict);
    end
end

% 3. merge distance of grey incidence alansys and domain set index
for i = 1 : qlen
    dist_GIA = GreyIncidenceDegree(queryPsepssm(i,:),psepssm);
    dist(i,:) = (dist_GIA + dist_DSI(i,:))/2;
    %[B,I] = sort(dist,'descend');
    %family{i} = familyId(I(1));
    %prob{i} = B(I(1));
end
save independ_test_dist.mat dist
   
