function [outputdata] = noname_elo(data,tokill)
%FONCTION POUR NORMALIZE LES NOMS DES ANIMAUX DANS DIFFERENT DATA SET
% et ENLEVER TEST_HM

try
    % REMOVING TEST_HM
    allmk=unique(data.PlayingSubject);
    if sum(strcmp('Test_HM',allmk))
        % REMOVING TEST_HM
        mask=strcmp('Test_HM',data.PlayingSubject);
        data(mask,:)=[];
        allmk=unique(data.PlayingSubject);
    end
    
catch
end

try
    % REMOVING TEST_HM
    allmk=unique( [ data.Monkey data.ConflictingMonkey ] );
    if sum(strcmp('Test_HM',allmk))
        % REMOVING TEST_HM
        mask=strcmp('Test_HM',data.Monkey) ||  strcmp('Test_HM',data.ConflictingMonkey) ;
        data(mask,:)=[];
        allmk=unique( [ data.Monkey data.ConflictingMonkey ] );
    end
catch
end

%%%%%%%%%%%%%%%%%%%%%%
% % %      NAME BANK%%%%
Z=1;
realmk{Z}={'abr', 'Abricot', ' Abr','abricot','Abr'}; Z=Z+1;
realmk{Z}={'ala', 'Alaryc', 'Alr', ' Ala',' Alr','alaryc','Ala'}; Z=Z+1;
realmk{Z}={'alv', 'Alvin', ' Alv','alvin','Alv'}; Z=Z+1;
realmk{Z}={'anu', 'Anubis', ' Anu', 'anubis','Anu'}; Z=Z+1;
realmk{Z}={'bar', 'Baranabe', 'Barnabe', ' Bar', 'Barnabé','Barnab�','barnabe','Bar','BarnabÃ©'}; Z=Z+1;
realmk{Z}={'ber', 'Berenice', ' Ber', 'Bérénice','B�r�nice','berenice','Ber','BÃ©rÃ©nice'}; Z=Z+1;
realmk{Z}={'ces', 'Ceasar', ' Ces', 'César','C�sar','Cesar','Ces','CÃ©sar'}; Z=Z+1;
realmk{Z}={'dor', 'Dory', ' Dor','Dor', 'dory'}; Z=Z+1;
realmk{Z}={'eri', 'Eric', ' Eri','Eri'}; Z=Z+1;
realmk{Z}={'jea', 'Jeanne', ' Jea','jeanne','Jea'}; Z=Z+1;
realmk{Z}={'lad', 'Lady', ' Lad','lady','Lad'}; Z=Z+1;
realmk{Z}={'las', 'Lassa', ' Las','lassa','Las'}; Z=Z+1;
realmk{Z}={'nem', 'Nema', ' Nem', 'Néma', 'N�ma','nema','Nem','NÃ©ma'}; Z=Z+1;
realmk{Z}={'nen', 'Nenno', 'Nen'}; Z=Z+1;
realmk{Z}={'ner', 'Nereis', ' Ner', 'Néréis','N�r�is','nereis', 'Ner','NÃ©rÃ©is'}; Z=Z+1;
realmk{Z}={'ola', 'Olaf', ' Ola','olaf', 'Ola'}; Z=Z+1;
realmk{Z}={'olg', 'Olga', ' Olg','olga', 'Olg' }; Z=Z+1;
realmk{Z}={'oli', 'Olli', ' Oll', 'Oll','ollie', 'Oli','oll'}; Z=Z+1;
realmk{Z}={'pac', 'Patchouli', ' Pac','patchouli', 'Pac'}; Z=Z+1;
realmk{Z}={'pat', 'Patsy', ' Pat','ptsy','patsy', 'Pat'}; Z=Z+1;
realmk{Z}={'wal', 'Wallace', ' Wal', 'Wal'}; Z=Z+1;
realmk{Z}={'wat', 'Walt', ' Wat', 'Wat'}; Z=Z+1;
realmk{Z}={'wot', 'Wotan', ' Wot', 'Wot'}; Z=Z+1;
realmk{Z}={'yan', 'Yang', ' Yan', 'Yan'}; Z=Z+1;
realmk{Z}={'yak', 'Yannick', ' Yak', 'Yak'}; Z=Z+1;
realmk{Z}={'yin', ' Yin','yin','Yin', 'yi'}; Z=Z+1;
realmk{Z}={'yoh', ' Yoh','yoh','Yoh'}; Z=Z+1;
realmk{Z}={'ult', 'Ulysse', ' Uly','ulysse', 'Uly'}; Z=Z+1;
realmk{Z}={'fic', 'Ficelle', 'Fic'}; Z=Z+1;
realmk{Z}={'gan', 'Gandhi'}; Z=Z+1;
realmk{Z}={'hav', 'Havana', 'Havane','Havanna'}; Z=Z+1;
realmk{Z}={'con', 'Controle'}; Z=Z+1;
realmk{Z}={'gai', 'Gaia'}; Z=Z+1;
realmk{Z}={'her', 'Hercules'}; Z=Z+1;
realmk{Z}={'han', 'Hanouk'}; Z=Z+1;
realmk{Z}={'hor', 'Horus'}; Z=Z+1;
realmk{Z}={'hor', 'Horus'}; Z=Z+1;
realmk{Z}={'imo', 'Imoen'}; Z=Z+1;
realmk{Z}={'ind', 'Indigo'}; Z=Z+1;
realmk{Z}={'iro', 'Iron'}; Z=Z+1;
realmk{Z}={'jip', 'Jipsy'}; Z=Z+1;

%Rhesus
realmk{Z}={'qui', 'Quinoa'}; Z=Z+1;
realmk{Z}={'sam', 'Samael'}; Z=Z+1;
realmk{Z}={'spl', 'Spliff'}; Z=Z+1;
realmk{Z}={'the', 'Theoden'}; Z=Z+1;
realmk{Z}={'spl', 'Spliff'}; Z=Z+1;
realmk{Z}={'ulr', 'Ulysse_Rhesus'}; Z=Z+1;
realmk{Z}={'vol', 'Volga'}; Z=Z+1;
realmk{Z}={'vla', 'Vladimir'}; Z=Z+1;
realmk{Z}={'any', 'Anyanka'}; Z=Z+1;
realmk{Z}={'arw', 'Arwen'}; Z=Z+1;
realmk{Z}={'baa', 'Baal'}; Z=Z+1;
realmk{Z}={'bor', 'Boromir'}; Z=Z+1;
realmk{Z}={'djo', 'Djocko'}; Z=Z+1;
realmk{Z}={'eow', 'Eowyn'}; Z=Z+1;
realmk{Z}={'far', 'Faramir'}; Z=Z+1;
realmk{Z}={'kab', 'Kabuki'}; Z=Z+1;
realmk{Z}={'nat', 'Natasha'}; Z=Z+1;
realmk{Z}={'yel', 'Yelena'}; Z=Z+1;
realmk{Z}={'yva', 'Yvan'}; Z=Z+1;

realmk=realmk';

disp('...normalazing Monkeys names')

try
    % TRANSFORM DATA NAME MONKEY
    for mk = 1:size(realmk,1)
        for n = 2:numel(realmk{mk})
            mask=strcmp(realmk{mk}{n},data.PlayingSubject);
            data.PlayingSubject(mask)=realmk{mk}(1);
        end
    end
    
    %%%%% REMOVING MONKEYS FROM DATASET
    for mk=1:numel(tokill)
        disp(['... removing monkey: ' tokill{mk}])
        mask=strcmp(tokill{mk},data.PlayingSubject);
        data(mask,:)=[];
    end
    
    try
        % TRANSFORM DATA NAME MONKEY
        for mk = 1:size(realmk,1)
            for n = 2:numel(realmk{mk})
                mask=strcmp(realmk{mk}{n},data.ConflictingSubject);
                data.ConflictingSubject(mask)=realmk{mk}(1);
            end
        end
        
        %%%%% REMOVING MONKEYS FROM DATASET
        for mk=1:numel(tokill)
            disp(['... removing monkey: ' tokill{mk}])
            mask=strcmp(tokill{mk},data.ConflictingSubject);
            data(mask,:)=[];
        end
    catch
    end
    
catch
end


% try data=cell2mat(data); catch; end

try
    % TRANSFORM DATA NAME MONKEY
    for mk = 1:size(realmk,1)
        for n = 2:numel(realmk{mk})
            mask=strcmp(realmk{mk}{n},data);
            data(mask)=realmk{mk}(1);
        end
    end
    
    
    %%%%% REMOVING MONKEYS FROM DATASET
    for mk=1:numel(tokill)
        disp(['... removing monkey: ' tokill{mk}])
        mask=strcmp(tokill{mk},data) ;
        data(mask,:)=[];
    end
catch
    
end

try
% TRANSFORM DATA NAME MONKEY
for mk = 1:size(realmk,1)
    for n = 2:numel(realmk{mk})
        mask=strcmp(realmk{mk}{n},data.Monkey);
        data.Monkey(mask)=realmk{mk}(1);
        mask=strcmp(realmk{mk}{n},data.ConflictingMonkey);
        data.ConflictingMonkey(mask)=realmk{mk}(1);
    end
    
end

catch

end


outputdata=data;


end

