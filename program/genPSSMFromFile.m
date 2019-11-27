function p=genPSSMFromFile(fastaFile)
% version 2 modified on 10/8/2015.
% Add the blosum62 substitution matrix for protein that cann't be
% alignmented in db.



blosumMatrix = blosum(62);% get blosum62 matrix
aaletters = 'ARNDCQEGHILKMFPSTWYVBZX*';

% cmd = ['psiblast -db ' db ' -query input.fasta -num_iterations 3 -threshold 0.01 -word_size 2 -evalue 10 -out_ascii_pssm pssmresult'];%pssmresult'];
[heads,seqs]=fastaread(fastaFile);

if ~ischar(heads)
    row = length(heads);
    p = cell(1,row);

    for k = 1 : row
        disp(k);
%         fid_in = fopen('input.fasta','w');
%         fprintf(fid_in,seqs{k});
%         fclose(fid_in);
%         
%         %if pssmresult is existed, delete it.
%         if exist('pssmresult','file') > 0
%             delete('pssmresult');
%         end
        
        %perform blast program
%         [~,~]=system(cmd);
        
        M = zeros(length(seqs{k}),40);
        acc = split(heads{k});
        pssmfile=strjoin({'H:\\Mywork\\pssm\\',acc{1},'pssm.txt'},'');
        fid_out = fopen(pssmfile,'r');
        if fid_out == -1
            sequence = upper(seqs{k});
            for i = 1 : length(sequence)
                t = strfind( aaletters, sequence(i));
                M(i,1:20) = blosumMatrix(t,1:20);
            end
            p{k} = M;
        else
            for i = 1 : 4
                tline=fgetl(fid_out);
            end
            i = 1;
            while ischar(tline) && ~isempty(tline)
                A = sscanf(tline,'%d %s %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d');
                M(i,:) = A(3:42)';
                i=i+1;
                tline=fgetl(fid_out);
            end
            p{k}=M;
            fclose(fid_out);
        end
    end
else %just only one sequence
    p = zeros(length(seqs),40);
    
    fid_in = fopen('input.fasta','w');
    fprintf(fid_in,seqs);
    fclose(fid_in);
    
    %if pssmresult is existed, delete it.
    if exist('pssmresult','file') > 0
        delete('pssmresult');
    end
    
    %perform blast program
    system(cmd);
      
    fid_out = fopen('pssmresult','r');
    
    if fid_out ~= -1 %if pssmresult is existed
        for i = 1 : 4
            tline = fgetl(fid_out);
        end
        i = 1;
        while ischar(tline) && ~isempty(tline)
            A = sscanf(tline,'%d %s %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d');
            p(i,:) = A(3:42)';
            i = i + 1;
            tline = fgetl(fid_out);
        end
        fclose(fid_out);
    else %if pssmresult isn't existed
        for i = 1 : length(seqs)
            t = strfind(aaletters, seqs(i));
            p(i,1:20) = blosumMatrix(t,1:20);
        end
    end
end