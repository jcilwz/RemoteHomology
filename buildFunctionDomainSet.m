function pfams = buildFunctionDomainVect(dataFile)
	% load mapObj
	load('uniprot_seqs_dict.mat')
	
    [dheads,dseqs] = fastaread(dataFile);    
    pfams = cell(1,length(dheads));
    
    hmmscanCMD = 'hmmscan -o out.txt --tblout fmout.tbl --acc --noali Pfam-A.hmm input.fasta';
    for i = 1 : length(dheads)
        i
        pfam = [];
        h = getHomoProteinsByHMM(dheads{i},dseqs{i});
        if isempty(h)
            fprintf('%s has not homology proteins\n', dheads{i});
            continue;
        end
        for j = 1 : length(h)
            pid = char(h(j));%homology protein
            pseq = mapObj(pid);
            if exist('out.txt','file') > 0
                delete('out.txt');
            end
            if exist('fmout.tbl','file') > 0
                delete('fmout.tbl');
            end
            if exist('input.fasta','file') > 0
                delete('input.fasta');
            end
            fastawrite('input.fasta',pid,pseq);
            [~,~] = system(hmmscanCMD);
            
            fid = fopen('fmout.tbl','r');
            flag = 1;
            while flag 
                tline = fgetl(fid);
                if contains(tline,'[ok]')%the end of file
                    break;
                end
                if startsWith(tline,'#')
                    if flag == 1
                        continue;
                    else
                        break
                    end
                else
                    flag = 0;
                    s = split(tline);
                    pf = extractBetween(s(2),1,7);
                    pfam = [pfam; pf];
                end
            end
            fclose(fid);
            pfams{i} = pfam;
        end
    end



            
            
            
        
