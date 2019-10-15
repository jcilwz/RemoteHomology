load independ_test.mat
testFamilyId = familyId;
load 7329sequencesFamily.mat
trainFamilyId = familyId;
load independ_test_dist.mat
Row = 1846;
Auc50 = 0;
Auc1 = 0;
x=[1,1];
for i = 1 : Row
    testY = familyId{i};
    dis = dist(i,:);
    label_Y = strcmp(trainFamilyId,testY);
    Auc1 = Auc1 + AUCK(label_Y,dis,1,'descend');
    Auc50 = Auc50 + AUCK(label_Y,dis,50,'descend');
end
Auc50 = Auc50/Row;
Auc1 = Auc1/Row;
disp(Auc50)
disp(Auc1)
