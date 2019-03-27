function [CAB] = cabFx(varargin)

%take in all inputs and make definitions
NVal= varargin{1};
ND=NVal(1); NW=NVal(2); NB=NVal(3); NT=NVal(4);
NPNS= varargin{2};
NPart=NPNS(1); NSites=NPNS(2);
basis = varargin{3};
psiAllSave = varargin{4};

% subsystems only the half system subsystem
subss=[0];
ss=1;

% create bounds and basis
halfS=floor(size(basis,2)/2)+subss(ss);
fullS=floor(size(basis,2));
HilbD=size(basis,1);
expon=(fullS:-1:1)-1;


% loop for finding all correlators CAB

for ww=1:NW
for dd=1:ND
for bb=1:NB
    
%find rho for given parameters of W,D,B
temp_rho=abs(psiAllSave{ww,dd,bb})'.^2;

    for nn=1:NPart+1
        
        %declare bases for two halfs of system to find separability
        %distribution for different particle sectors "nn"
        bases=ones(size(basis,1),1)*(10.^(expon));
        baseKey=sum(basis.*bases,2);
        unqKey=unique(baseKey);
        halfS=floor(size(basis,2)/2)+subss(ss);

        %left half base
        baseHalfKeyL=sum(basis(:,1:halfS).*bases(:,1:halfS),2);
        baseHalfKeyR=sum(basis(:,halfS+1:end).*bases(:,halfS+1:end),2);
        halfNKey=sum(basis(:,1:halfS),2);

        %find indices of contributing states in particle number basis
        inds=find(halfNKey==(nn-1));
        binTemp=zeros(length(inds),3,NT);

        %find keys for looking at joint versus marginal distribution
        rhoAgB=zeros(length(unique(baseHalfKeyL(inds))),NT);
        unqAKey=unique(baseHalfKeyL(inds));
        unqBKey=unique(baseHalfKeyR(inds));
        
        %particle number normalization
        normP=sum(temp_rho(:,inds),2);

        for ii=1:length(inds)
            %joint distribution
            binTemp(ii,1,:)=reshape(temp_rho(:,inds(ii)),[1 NT])./normP';
            
            %marginal dist 1
            inds2=find(baseHalfKeyL(inds(ii))==baseHalfKeyL);
            binTemp(ii,2,:)=reshape(sum(temp_rho(:,inds2),2),[1 NT])./normP';
            
            %marginal dist 2
            inds3=find(baseHalfKeyR(inds(ii))==baseHalfKeyR);
            binTemp(ii,3,:)=reshape(sum(temp_rho(:,inds3),2),[1 NT])./normP';
        end
        
        % test separability of joint distribution
        cab_elem=reshape((binTemp(:,1,:)-binTemp(:,2,:).*binTemp(:,3,:)),[size(binTemp,1),NT]);
        
        % compute correlator for different number sectors
        for nt=1:NT
            testCorrSaveTH(nn,nt)=normP(nt)'.*sum(abs(cab_elem(:,nt)));
        end
        
    end

    corrThStruct{ww,dd,bb}=testCorrSaveTH;

end
end
end



% just reorganization for output of CAB
tempAllTh=zeros(NT,NPart+1,NW,NB,ND);
CAB=zeros(NT,NB,NW);


for ww=1:NW
for bb=1:NB
for dd=1:ND    
    tempAllTh(:,:,ww,bb,dd)=corrThStruct{ww,dd,bb}';  
end
end
end

%sum across different particle sectors and average over disorders
CAB=reshape(mean(sum(tempAllTh,2),5),[NT NW NB]);

end

