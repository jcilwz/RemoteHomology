load pfamA-Clans.mat
load 7329pfams.mat

N = 7329;
v = zeros(N,N);

for i = 1 : N
    i
    for j = i+1 : N 

        if j ~= i
            s = pfamcmp(pfams{i},pfams{j},clanDict);
            v(i,j) = s;
            v(j,i) = s;
        end
    end
end
save('7329FunctionDomainMatch_0.2.mat','v');
