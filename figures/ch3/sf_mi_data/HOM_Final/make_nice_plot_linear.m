function make_nice_plot(datax, datay, erry_l, erry_u, theoryx, theoryy, varargin)
% taking clue from here:
% http://blogs.mathworks.com/loren/2007/12/11/making-pretty-graphs/

if numel(varargin)>=3
    plottitle =     varargin{1};
    yaxislabel =    varargin{2};
    ymax =          varargin{3};
else
    plottitle =     'blank';
    yaxislabel =    'blank';
    ymax =          max(max(datay));
end

if numel(varargin)>3
    leg={varargin{4:end}};
end
% load nice_plot_test.mat

%--------------------------------------
% width = 9; % inch;
% height = width*0.7;
width = 9; % inch;
height = width*0.7;
clf
if ishandle(88)==0
    hFig = figure(88);
    set(hFig,'Units','inches','Resize','off','Position',[5,2,width,height])
end

% colors:

% compute upper and lower bounds on y from errors

% erry_lower = erry_l;
% erry_upper = erry_u;
erry_lower = datay - erry_l;
erry_upper = datay + erry_u;

% co = get(gca,'ColorOrder');
% set(gca, 'ColorOrder', co([2 3 6],:),'NextPlot','add');

%order, x, y, {lower x, upper x}, {lower y, upper y}
hdata = ploterr(datax, ...
    datay, ...
    [],...
    {erry_lower,erry_upper}, ...
    'o',...
    'abshhx', 0.03, 'abshhy', 0.035); %'hhxy',0.25);
%legend(leg,'Location','Northwest')
hold on
htheory = plot(theoryx,theoryy,'r'); % 'LineWidth',1

set(gcf,'color','w');
xlim([-0.01 0.5])
ylim([0 ymax])

% hTitle  = title (plottitle   );
hXLabel = xlabel('J t');
hYLabel = ylabel(yaxislabel  );

%hLegend = legend('site 1','site 2','full system')
%-------------------------------------
Ncurves = size(datax,2)
% datapoints
set(hdata(1:Ncurves)                            , ...
   'LineWidth'       , 2.5           , ...
  'Marker'          , 'o'         , ...
  'MarkerSize'      , 11          , ...
  'MarkerFaceColor' , [1 1 1] );

% y error & x error
set(hdata(Ncurves+1:2*Ncurves)                            , ...
    'LineWidth'       , 2.5            );
% set(hdata(2*Ncurves+1:3*Ncurves)                      , ...
%     'LineWidth'       , 1.5            );
%set(hdata(2*Ncurves+1:3*Ncurves)                            , ...
% theory curves
set(htheory                        , ...
  'LineWidth'           , 1.5    );

% reorder the layering of the curves
%uistack(findobj(gca, 'Color', 'r'), 'top')
uistack(findobj(gca, 'Color', 'r'), 'top')
uistack(findobj(gca, 'Color', 'b'), 'top')
% uistack(htheory, 'bottom')

%-------------------------------------
% fonts and text size
set( gca                       , ...
    'FontName'   , 'Helvetica' );
set([gca]             , ...
    'FontSize'   , 18           ); % axis numbers
set([hXLabel, hYLabel]  , ...
    'FontSize'   , 20          );
% set( hTitle                    , ...
%     'FontSize'   , 20          , ...
%     'FontWeight' , 'bold'      );

% axis, ticks
set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'      , ...
  'XTick'       , 0:0.125:0.5, ...
  'YTick'       , 0:0.2:ymax, ...
  'XTickLabels' , {'0','1/8','1/4','3/8','1/2'},...
  'LineWidth'   , 1.5         );


% save to eps
set(gcf, 'PaperPositionMode', 'auto');
print('-depsc2', ['nice_plots/' plottitle '.eps'])
% close;

end