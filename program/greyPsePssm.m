%greyPsePssm(pssmCell,n)由pssm矩阵生成伪pssm矩阵
%pssmCell是1*L的cell型数据，存放L个蛋白质序列的PSSM矩阵
%pssmCell{i}表示第i个蛋白质的PSSM矩阵。
%返回一个L*M的矩阵（M=60或80)，矩阵的每一行向量表示一个蛋白质序列的伪pssm向量，
%该向量分两部分：前20列是PSSM矩阵的按列方向的均值,
%后40或60列由灰色模型产生。
%IF n=1,则由GM(1,1)模型产生的灰色系数作为伪氨基酸成分
%IF n=2,则由GM(2,1)模型产生的灰色系数作为伪氨基酸成分
%具体方法是以PSSM矩阵的每一列为序列构建灰色模型，由于PSSM中元素的值-7~7，
%所以建模前以函数1/(1+exp(-x))处理。
%所得的伪PSSM向量没有标准化，也没有对灰色系数考虑设置权值，这些有使用者在函数外出来
function psepssm=greyPsePssm(pssmCell,n)

row = length(pssmCell);

if n == 1 %GM(1,1) generate pseudo 
    psepssm = zeros(row,60);
    for i = 1 : row
        pssm = 1./(1+exp(-pssmCell{i}));
        psepssm(i,1:20) = mean(pssm(:,1:20));
        for j = 1 : 20
           p = GMParam(pssm(:,j));
           psepssm(i,20+2*j-1:20+2*j) = [abs(p(1)) abs(p(2))];
        end
%         psepssm(i,:) = psepssm(i,:)/sum(psepssm(i,:));
    end
elseif n == 2 %GM(2,1) generate pseudo
    psepssm = zeros(row,80);
    for i = 1 : row
        pssm = 1./(1+exp(-pssmCell{i}));
        psepssm(i,1:20) = mean(pssm(:,1:20));
        for j = 1 : 20
           p = GM21Param(pssm(:,j));
           psepssm(i,20+3*j-2:20+3*j) = [abs(p(1)) abs(p(2)) abs(p(3))];
        end
%         psepssm(i,:) = psepssm(i,:)/sum(psepssm(i,:));
    end
end