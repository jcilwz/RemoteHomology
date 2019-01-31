% pf1: 功能域集合1
% pf2: 功能域集合2
% clanDict: 功能域宗族归属字典，cell类型，序号为功能域的序号。
%           每一项是一个字符串cell,内容依次为：clan_id,clan_name,clan_abb,
%           clan_annotation
function s = pfamcmp(pf1, pf2, clanDict)
s = 0;
if isempty(pf1) || isempty(pf2)
    %两个功能域中有一个为空集，返回0
    return
end
%比较两个功能域的类似程度
intset = intersect(pf1,pf2);
uniset = union(pf1,pf2);
if ~isempty(intset)
    %两个功能域有交集
    s = size(intset,1)/size(uniset,1);
    return;
else
    %两个功能域没有交集，判断是否有相同的宗族
    for p1 = pf1
        m = extractBetween(p1,3,7);
        m = str2num(char(m));
        clan1 = clanDict{m};
        for p2 = pf2
            m = extractBetween(p2,3,7);
            m = str2num(char(m));
            clan2 = clanDict{m};
            if strcmp(clan1{1},'NONE') || strcmp(clan2{1}, 'NONE')
                s = 0;
            else
                %两个功能域都有宗族注释
                if strcmp(clan1{1},clan2{1})
                    s = 0.2;
                end
            end
        end
    end
    return
end
            