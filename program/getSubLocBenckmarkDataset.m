load pfamA-Clans.mat
proteinFile = 'SubLoc.fasta';
[heads, seqs] = fastaread(proteinFile);
nb_seqs = length(heads);
pssm = blastpssm(proteinFile, 'swissprot');
Psepssm = greyPsePssm(p,2);
save('Psepssm.mat',Psepssm)

dist_DSI = zeros(nb_seqs, nb_seqs);
protFams = buildFunctionDomainSet(proteinFile);
for i = 1 : nb_seqs
    dist_DSI(i,i) = 1;
    for j = i+1 : nb_seqs
        dist_DSI(i,j) = pfamcmp(protFams{i}, protFams{j}, clanDict);
    end
end
save('DSI_dist.mat',dist_DSI)