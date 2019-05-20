function [AvgTrapzOut, TrapzOut, TrapzErr, DWout] = TrapzIntErr(dat)

%find the trapezoid area of some points and compute their error
% data should be in column form of [X, Y, E]

%they should be sorted but let's go ahead and do that incase it is wrong

if size(dat,2)==3
    [ii,jj]=sort(dat(:,1));

    dat=[dat(jj,1), dat(jj,2), dat(jj,3)];

    TZdat = [dat(2:end,2)+dat(1:end-1,2)]./2;
    TZdatErr = sqrt([dat(2:end,3).^2+dat(1:end-1,3).^2]./4);
    DW = [dat(2:end,1)-dat(1:end-1,1)];

    TrapzOut = tril(ones(size(TZdat,1)))*(TZdat.*DW);
    TrapzErr = sqrt( ...
        tril(ones(size(TZdat,1)))* ...
        ((TZdatErr.*DW)).^2 ...
        );

    AvgTrapzOut = TrapzOut./(dat(end,1)-dat(1,1));

    %DWout = DW./2+dat(1:end-1,1);
    DWout = dat(2:end,1);

else
    
    [ii,jj]=sort(dat(:,1));

    dat=[dat(jj,1), dat(jj,2)];

    TZdat = [dat(2:end,2)+dat(1:end-1,2)]./2;
    %TZdatErr = sqrt([dat(2:end,3).^2+dat(1:end-1,3).^2]./4);
    DW = [dat(2:end,1)-dat(1:end-1,1)];

    TrapzOut = tril(ones(size(TZdat,1)))*(TZdat.*DW);
    TrapzErr = zeros(size(TrapzOut));
    
    AvgTrapzOut = TrapzOut./(dat(end,1)-dat(1,1));

    %DWout = DW./2+dat(1:end-1,1);
    DWout = dat(2:end,1);
    
end

end