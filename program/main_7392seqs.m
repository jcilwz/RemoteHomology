%function [Auc50, Auc1]=main_7392seqs()

load 7329sequencesGrey21PSSM.mat
%load 3187sequences_GLCMFeatures.mat
%load 3187sequences_blastscore.mat 
load 7329sequencesFamily.mat
%load 3187sequences_fuzadaFeatures.mat
load 7329FunctionDomainMatch_0.2.mat

Row = 7329;
Auc50 = 0;
Auc1 = 0;
x=[1,1];
for i = 1 : Row
    ind = true(Row,1);
    ind(i) = false;
    trainX_psepssm = psepssm(ind,:);
    
    trainY = familyId(ind);
    
    testX_psepssm = psepssm(i,:);
    testY = familyId{i};
    
    label_Y = strcmp(trainY,testY);
%     r = GreyIncidenceDegree([testX_psepssm testX_glcm],[trainX_psepssm trainX_glcm]);
    r = GreyIncidenceDegree(testX_psepssm,trainX_psepssm);
%     d = GreyIncidenceDegree(testX_glcm, trainX_glcm);    
%     b = mapminmax(bitscore(i,:),0,1);
%      f = distance(testX_fzd,trainX_fzd);
%     f = GreyIncidenceDegree(testX_fzd,trainX_fzd);
%    f = mapminmax(f',0,1);
%     dis=x(1)*r + x(2)*b + x(3)*d + x(4)*(1-f)+v(ind,:);
    dis=x(1)*r + x(2)*v(i,ind);
    Auc1 = Auc1 + AUCK(label_Y,dis,1,'descend');
    Auc50 = Auc50 + AUCK(label_Y,dis,50,'descend');
end
Auc50 = Auc50/Row;
Auc1 = Auc1/Row;
disp(Auc50)
disp(Auc1)
