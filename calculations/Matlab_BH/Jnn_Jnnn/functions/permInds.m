function [KeyOut,IndsOut] = permInds(inds,intN)

    %intParts is a vector of the integer partitions that will be used to
    %organize the list of indices "inds"
    
    %inds and intParts should both be vectors
     KeyOut=(10.^(intN-1:-1:0)*nchoosek(inds,intN)')';
     IndsOut=nchoosek(inds,intN);

end

%{
function [prodOut] = permIndsProd(dataIn,prodIn,knum,inds)
    
    if size(inds,2)~=0
        nkInds=nchoosek(inds,knum(1));
        
        for nk=1:size(nkInds,1)            
            indsOut=setdiff(inds,nkInds(nk,:));
            prodOut=GnCorrDatFx(dataIn,nkInds(nk,:))*permIndsProd(dataIn,prodIn,knum(2:end),indsOut);
        end
        
    else
        prodOut=prodIn;
    end

end
%}