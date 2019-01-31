
%--X: X={x1,x2,...,xm}为灰关联因子集.类型Matrix
%--y: y={y(0),y(1),...,y(n)}参考序列
%--xi:xi={xi(0),xi(1),...,xi(n)}比较序列
function R =GreyIncidenceDegree(y,X)
    m = size(X,1);
    n = length(y);
    D = zeros(m,n);%距离矩阵
    z = zeros(m,n);%关联系数矩阵
    R = zeros(1,m);%关联度序列
    p=0.5;%分辨系数
    %求各序列的初值像
    %y = y/y(1);
    %for k = 1:m
    %    a = X(k,1);
    %	X(k,:) = X(k,:)/a;
    %end
    %计算关联度
    for k = 1:m
        D(k,:)=abs(y-X(k,:));
    end
    dmin = min(min(D));
    dmax = max(max(D));
    
    for i = 1:m
        for j = 1:n
            z(i,j) = (dmin + p*dmax)/(D(i,j)+p*dmax);
        end
    end
    
    for i = 1:m
          R(i) = sum(z(i,:))/n;
    end
