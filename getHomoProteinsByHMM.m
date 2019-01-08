function homoProteins=getHomoProteinsByHMM(head, seq)
    hmmbuildCMD = 'hmmbuild input.hmm input.fasta' ;
    hmmsearchCMD = 'hmmsearch input.hmm e:/uniprot_sprot.fasta > input.out';
 
    if exist('input.hmm','file') > 0
        delete('input.hmm');
    end
    if exist('input.fasta','file') > 0
        delete('input.fasta');
    end
    %if exist('input.out','file') > 0
    %    delete('input.out');
    %end
    
    fastawrite('input.fasta',head,seq);
    [~,~] = system(hmmbuildCMD);
    [~,~] = system(hmmsearchCMD);%output homology proteins
    
    homoProteins=[];
    fid = fopen('input.out','r');
    for i = 1 : 14
        fgetl(fid);
    end
    i = 1;
    
    while i <= 10
        tline = fgetl(fid);
        if isempty(tline)
            break;
        end
        r = strfind(tline,'----');%no homology proteins
        if ~isempty(r)
            break;
        end
        s = textscan(tline,'%f %f %f %f %f %f %f %d %s %s');
        homoProteins=[homoProteins;s{9}];
        i = i + 1;
    end
    fclose(fid);

        