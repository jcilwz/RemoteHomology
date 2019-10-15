foldname = './SCOP167-superfamily/';
filelist = dir([foldname   'pos-train*']);
N = length(filelist);

for i = 1 : N
    inputfile = [foldname, filelist(i).name];
    s1 = 'cd-hit -i ';
    s2 = ' -o nr80 -c 0.8 -n 5 -d 0 -M 16000 -T 16';
    s3 = ' nr80 -o nr60 -c 0.6 -n 4 -d 0 -M 16000 -T 16';
    s4 = [' nr60 -o ', 'nr40-', filelist(i).name, ' -c 0.4 -n 2 -d 0 - M 16000 -T 16'];
    cmd_cdhit_1 = [s1, inputfile, s2];
    cmd_cdhit_2 = [s1, s3];
    cmd_cdhit_3 = [s1, s4];
    [~,~]=system(cmd_cdhit_1);
    [~,~]=system(cmd_cdhit_2);
    [~,~]=system(cmd_cdhit_3);
end

filelist = dir([foldname   'nr40*']);
[pos_h,pos_s] = fastaread([foldname   filelist(1).name]);

for i = 2 : length(filelist)
    [ht,st] = fastaread([foldname   filelist(i).name]);
    pos_h = union(pos_h,ht);
    pos_s = union(pos_s,st);
end