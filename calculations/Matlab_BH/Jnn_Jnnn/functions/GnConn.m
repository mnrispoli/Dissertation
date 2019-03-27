%function [gnConn, gnDisConn] = GnConn(varargin)
function [gnConn] = GnConn(varargin)

% find connected vs disconnected correlation functions 
% monte carlo data
dataRaw=varargin{1};
inds=varargin{2};
NInd=length(inds);

if NInd==0
    gnConn=1;
else

% the basis produced for unity filling is actually a convenient way to find
% all unique integer partitions. Will fail for very large system sizes due
% to tryting to make the enormous basis
intParts=flipud( ...
                unique( ...
                    sort(  BasisMake(NInd,NInd),2,'descend'  ) ...
                ,'rows') ...
                );

% initialize disconnected and connected correlators
gnDisConn=GnCorrDatFx(dataRaw,inds);
gnConn=gnDisConn;

if size(intParts,1)>1
% find number of unique integer partitions. we skip the first one below
% because it is already the disconnected part+connected part
intParts=intParts(2:end,:);
NIntPrt=size(intParts,1);

indsRef=1:NInd;

IndsSave=cell(1,NIntPrt);
KeySave=cell(1,NIntPrt);

for ni=1:NIntPrt
    
    KeyIn=[];
    IndsIn=[];
    intPrt=intParts(ni,:);
    inds=indsRef;
    
    for ki=1:NInd
        if intPrt(ki)>0
        if size(KeyIn,1)>0
            
            tmpKey=[];
            tmpInds=[];
            
            for jj=1:size(KeyIn,1)
                
                inds=setdiff(indsRef, ....
                            IndsIn(jj,1:sum(intPrt(1:ki-1))) ...
                            );
                        
                [KeyOut,IndsOut]=permInds(inds,intPrt(ki));
                prefixKey=ones(size(KeyOut,1),1)*[KeyIn(jj,1:ki-1)];
                tmpKey=[tmpKey; prefixKey KeyOut];
                
                
                prefixInds=ones(size(IndsOut,1),1)*[IndsIn(jj,1:sum(intPrt(1:ki-1)))];
                tmpInds=[tmpInds; prefixInds IndsOut];
            end
            
            KeyIn=tmpKey;
            IndsIn=tmpInds;
            
        else    
            [KeyOut,IndsOut]=permInds(inds,intPrt(ki));
            KeyIn=[KeyIn KeyOut];
            IndsIn=[IndsIn IndsOut];
        end
        end
    end
    
    [unqKey,unqInds]=unique(sort(KeyIn,2),'rows');
    KeySave{ni}=unqKey;
    IndsSave{ni}=IndsIn(unqInds,:);
    
    
    for kk=1:size(IndsSave{ni},1)
        
        tmpInds=IndsSave{ni};
        tmpProd=1;
        
        pinds=find(intParts(ni,:)>0);
        indsBnds=[0; tril(ones(NInd))*intParts(ni,:)'];
        
        for ppi=1:length(pinds)
            tmpProd=tmpProd.*GnConn( ...
                                    dataRaw, ...
                                    tmpInds(kk,indsBnds(ppi)+1:indsBnds(ppi+1)) ...
                                    );
        end
        
        gnConn=gnConn-tmpProd;
    end

end
end
end

end