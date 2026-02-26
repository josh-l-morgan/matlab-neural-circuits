function hh = myPie(p)

%% modified to plot onto particular position on chart
% p.vals = vals
% p.color = [1 0 1; 0 1 0; 0 0 0]/3;
% p.center = [meany meanx];
% p.radius = 300;
%   
%
%PIE    Pie chart.
%   PIE(X) draws a pie plot of the data in the vector X.  The values in X
%   are normalized via X/SUM(X) to determine the area of each slice of pie.
%   If SUM(X) <= 1.0, the values in X directly specify the area of the pie
%   slices.  Only a partial pie will be drawn if SUM(X) < 1.
%
%   PIE(X,EXPLODE) is used to specify slices that should be pulled out from
%   the pie.  The vector EXPLODE must be the same size as X. The slices
%   where EXPLODE is non-zero will be pulled out.
%
%   PIE(...,LABELS) is used to label each pie slice with cell array LABELS.
%   LABELS must be the same size as X and can only contain strings.
%
%   PIE(AX,...) plots into AX instead of GCA.
%
%   H = PIE(...) returns a vector containing patch and text handles.
%
%   Example
%      pie([2 4 3 5],{'North','South','East','West'})
%
%   See also PIE3.

%   Clay M. Thompson 3-3-94
%   Copyright 1984-2005 The MathWorks, In
%   $Revision: 1.16.4.11 $  $Date: 2011/07/25 03:49:31 $


%%
% Parse possible Axes input

cax = [];

if isfield(p,'vals')
    x = p.vals(:); % Make sure it is a vector
else 
    x = 1;
end

if isfield(p,'center')
    center = p.center;
else
    center = [0 0];
end

if isfield(p, 'color')
    colors = p.color;
else
    colors =jet(length(x));
end

if isfield(p, 'text')
    txtlabels = p.text;
    useLabels = 1;
else
    txtlabels = '';
    useLabels = 0;
end

if isfield(p,'explode')
    explode = p.explode;
else
    explode = x * 0;
end


if isfield(p,'radius')
    circleScale = p.radius;
else
    circleScale = 1;
end

% nonpositive = (x <= 0);
% if all(nonpositive)
%     error(message('MATLAB:pie:NoPositiveData'));
% end
% if any(nonpositive)
%   warning(message('MATLAB:pie:NonPositiveData'));
%   %x(nonpositive) = [];
% end
xsum = sum(x);
if xsum > 1+sqrt(eps), x = x/xsum; end

% % Look for labels
% if nargs>1 && iscell(args{end})
%   txtlabels = args{end};
%   if any(nonpositive)
%     txtlabels(nonpositive) = [];
%   end
%   args(end) = [];
% else
%   for i=1:length(x)
%     if x(i)<.01,
%       txtlabels{i} = '< 1%';
%     else
%       txtlabels{i} = sprintf('%d%%',round(x(i)*100));
%     end
%   end
% end

% Look for explode
   explode = zeros(size(x)); 
% 
% if isempty(args),
% else
%    explode = args{1};
%    if any(nonpositive)
%      explode(nonpositive) = [];
%    end
% end
% 
% explode = explode(:); % Make sure it is a vector
% 
% if length(txtlabels)~=0 && length(x)~=length(txtlabels),
%   error(message('MATLAB:pie:StringLengthMismatch'));
% end
% 
% if length(x) ~= length(explode),
%   error(message('MATLAB:pie:ExploreLengthMismatch'));
% end

% 
% cax = gcf;
% if isempty(cax)
% end

%%
%cax = newplot(gca);
cax = newplot(cax);
next = lower(get(cax,'NextPlot'));
hold_state = ishold(cax);


shiftTheta = x(1)*pi*2;
theta0 = pi/2-shiftTheta;
maxpts = 100;
inside = 0;
h = [];


%% outside
n = max(1,ceil(maxpts*1));
r = [ones(n+1,1)];
theta = pi/2 + [1*(0:n)'/n]*2*pi;
if inside,
    [xtext,ytext] = pol2cart(theta0 + 1*pi,.5);
else
    [xtext,ytext] = pol2cart(theta0 + 1*pi,1.2);
end
[xx,yy] = pol2cart(theta,r);
xx = xx * circleScale;
yy = yy * circleScale;


xx = xx + center(2);
yy = yy + center(1);

i = length(x);
outsideColor = (colors(1,:) * double(x(1)>0)) + (colors(2,:) * double(x(2)>0));
%outsideColor(outsideColor==0) = .5;
% if sum(outsideColor>0) < 3;
%     outsideColor = outsideColor + .5;
% end
outsideColor = outsideColor+.3;
outsideColor(outsideColor>1) = 1;

h = [h,patch('XData',xx,'YData',yy,'CData',i*ones(size(xx)), ...
    'FaceColor',outsideColor,'FaceAlpha',0,'EdgeColor',outsideColor,'LineWidth',1,'parent',cax) ];

%%




for i=1:length(x)
    if x(i)
  n = max(1,ceil(maxpts*x(i)));
  r = [0;ones(n+1,1);0];
  theta = theta0 + [0;x(i)*(0:n)'/n;0]*2*pi;
  if inside,
    [xtext,ytext] = pol2cart(theta0 + x(i)*pi,.5);
  else
    [xtext,ytext] = pol2cart(theta0 + x(i)*pi,1.2);
  end
  [xx,yy] = pol2cart(theta,r);
  xx = xx * circleScale;
  yy = yy * circleScale;
  if explode(i),
    [xexplode,yexplode] = pol2cart(theta0 + x(i)*pi,.1);
    xtext = xtext + xexplode;
    ytext = ytext + yexplode;
    xx = xx + xexplode;
    yy = yy + yexplode;
  end
  theta0 = max(theta);
  
  
  xx = xx + center(2);
  yy = yy + center(1);
  if useLabels
  h = [h,patch('XData',xx,'YData',yy,'CData',i*ones(size(xx)), ...
               'FaceColor',colors(i,:),'parent',cax), ...
         text(xtext,ytext,txtlabels{i},...
              'HorizontalAlignment','center','parent',cax)];
  else
        h = [h,patch('XData',xx,'YData',yy,'CData',i*ones(size(xx)), ...
               'FaceColor',colors(i,:),'EdgeColor',[1 1 1],'EdgeAlpha',0,'LineWidth',.001,'parent',cax) ];
      
  end
    end
end






%
% if ~hold_state,
%   view(cax,2); set(cax,'NextPlot',next);
%   axis(cax,'equal','off',[-1.2 1.2 -1.2 1.2])
% end

if nargout>0, hh = h; end
% 
% % Register handles with m-code generator
% if ~isempty(h)
%   mcoderegister('Handles',h,'Target',h(1),'Name','pie');
% end


%}

