load 3187FunctionDomainMatch_0.1.mat
load 3187sequences_Grey21PSSM.mat %psepssm
load 3187sequencesFamily.mat
Auc50 = 0;
Auc0 = 0;
for i = 1 : 3187
    ind = true(3187,1);
    ind(i) = false;
    
    trainY = familyId(ind);
    label_Y = strcmp(trainY,familyId{i});
    
    trainX_psepssm = psepssm(ind,:);
    
    dist1 = v(i,ind);
    dist2 = GreyIncidenceDegree(psepssm(i,:),trainX_psepssm);
    
    dist = dist1 + dist2;
    Auc50 = Auc50 + AUCK(label_Y,dist,50,'descend');
    Auc0 = Auc0 + AUCK(label_Y,dist,0,'descend');
end
disp(Auc50/3187);
disp(Auc0/3187);