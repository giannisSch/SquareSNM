%   Author : Ioannis Schoinochoritis
%   Email : g.sxoin@hotmail.com
%   Version : 1.0, November. 2018
%   In this code i try to calculate the side of largest square between 2 curves.
%   The side of largest square is static noise margin of a SRAM memory 

%load data from  csv files 
curve1=csvread('C:/Users/gsxoi/curve1.csv');
curve2=csvread('C:/Users/gsxoi/curve2.csv');


%if Write SNM change value for writeOperation = 1 ,else 0
writeOperation=0;

%if data of curve2 is similar with curve 1 uncomment
%following line to rotate the curve2  
%curve2=rotateAtXAxis(curve2);

%extract more Data 
% if you have enough data please comment following lines
curve1=extractNewData(curve1);
curve1=extractNewData(curve1);
curve1=extractNewData(curve1);
curve1=extractNewData(curve1);
curve2=extractNewData(curve2);
curve2=extractNewData(curve2);
curve2=extractNewData(curve2);
curve2=extractNewData(curve2);

%load data into vectors
curveX1= curve1(:, 1);
curveY1 = curve1(:, 2);
curveX2 = curve2(:, 1);
curveY2 = curve2(:, 2);

if writeOperation ==0 
    %Hold or Read SNMs
    P = InterX([curveX1';curveY1'],[curveX2';curveY2']); % find intersection points
    plot(curveX1,curveY1,curveX2,curveY2,P(1,:),P(2,:),'ro') %plot curves and points 
    hold on
    %-------------Up ----------------------
    %split graph in middle(P(2,2)) and find the largest square in Up graph (eye) 
    fprintf(" Up Square Sides \n");
    [newCurve,newCurveR]=splitCurves(curve1,curve2,P(2,2),'up');   
    curveX1= newCurve(:, 1);
    curveY1 = newCurve(:, 2);
    curveX2 = newCurveR(:, 1);
    curveY2 = newCurveR(:, 2);
    bestFitSquare(curveX1,curveY1,curveX2,curveY2);
    
    %----------Down------------
    %split graph in middle(P(2,2)) and find the largest square in Down graph (eye) 
    fprintf(" Down Square Sides \n");
    [newCurve,newCurveR]=splitCurves(curve1,curve2,P(2,2),'down');   
    n=length(newCurve);
    curveX1= newCurve(:, 1);
    curveY1 = newCurve(:, 2);
    curveX2 = newCurveR(:, 1);
    curveY2 = newCurveR(:, 2);
    bestFitSquare(curveX1,curveY1,curveX2,curveY2);
else 
    %we are now Write Static Noise Margins
    plot(curveX1,curveY1,curveX2,curveY2)
    hold on
    bestFitSquare(curveX1,curveY1,curveX2,curveY2);
end
hold off

function bestFitSquare(curveX1,curveY1,curveX2,curveY2)
    winnersArray=[];
    sizeOfFirstPlot=length(curveX1);
    sizeOfSecondPlot=length(curveX2);
     
    for i=1:sizeOfFirstPlot   
        x1=curveX1(i);
        y1=curveY1(i);
        for j=1:sizeOfSecondPlot 
            x2=curveX2(j);
            y2=curveY2(j);
            if checkAngle(x1,y1,x2,y2) %if angle between 2 points is 45 degrees
                [xSquare,ySquare]=createSquare(x1,y1,x2,y2); %create  square containing this 2 points
                interraction1 = InterX([curveX1';curveY1'],[xSquare;ySquare]); 
                interraction2 = InterX([curveX2';curveY2'],[xSquare;ySquare]) ;
                contacts1=length(interraction1);
                contacts2=length(interraction2);
                totalContacts=contacts1+contacts2;

                if (totalContacts>=2 && totalContacts<=4)    
                    %possible Square
                    distanceOfSide=abs(x1-x2);
                    winnersArray=[xSquare ySquare distanceOfSide ; winnersArray]; 
                end
            end
             
        end
    end
    [xWinnerSquare,yWinnerSquare]=findWinner(winnersArray);
    plot(xWinnerSquare,yWinnerSquare);
    
end

function [xWinnerSquare,yWinnerSquare]=findWinner(winnersArray)
    xWinnerSquare=[];
    yWinnerSquare=[];
   [~,columns] = size(winnersArray);
   while true
    disp('Please press something to Plot the Square.')
    disp('Press ESC to continue');

    w = waitforbuttonpress; 
    switch w 
        case 1 % (keyboard press) 
          key = get(gcf,'currentcharacter'); 
              switch key
                  case 27 % 27 is the escape key
                      disp('User pressed the escape key.')
                      break % break out of the while loop
                  otherwise 
                      [rows,columns] = size(winnersArray);
                      index=-1;
                      maxDistance=-2;
                        
                      for i=1:rows
                            if (winnersArray(i,11)>=maxDistance)
                                index=i;
                                maxDistance=winnersArray(i,11);
                                xWinnerSquare=[winnersArray(i,1) winnersArray(i,2) winnersArray(i,3) winnersArray(i,4) winnersArray(i,5)];
                                yWinnerSquare=[winnersArray(i,6) winnersArray(i,7) winnersArray(i,8) winnersArray(i,9) winnersArray(i,10)];
                            end;
                      end;
                      
                      winnersArray([index],:)=[];
                      plot(xWinnerSquare,yWinnerSquare);
                      
              end;
    end;
   end
    
   
    
    fprintf("Side x of Square is %f\n",abs(xWinnerSquare(1,1)-xWinnerSquare(1 ,3)));
    fprintf("Side y of Square is %f\n",abs(yWinnerSquare(1,1)-yWinnerSquare(1,3)));
    fprintf("Side of Square is %f\n",abs((xWinnerSquare(1,1)+yWinnerSquare(1,1))/2-(xWinnerSquare(1,3)+yWinnerSquare(1,3))/2));
        
end

function [moreDataPlot] = extractNewData(oldPlot)
    %extract new data from average value between 2 points
    moreDataPlot = [];
    size=length(oldPlot);
    for i=1:size-1
        x=(oldPlot(i,1)+oldPlot(i+1,1))/2;
        y=(oldPlot(i,2)+oldPlot(i+1,2))/2;
        existingData=[oldPlot(i,1) oldPlot(i,2)];
        newData=[x y];
        helpArray = [existingData ; newData];
        moreDataPlot = [moreDataPlot;helpArray];
    end
    moreDataPlot = [moreDataPlot;oldPlot(size,1) oldPlot(size,2)];
end        

function output = checkAngle(x1,y1,x2,y2)
    rad=atan2(y2-y1,x2-x1);
    degree = rad * (180/pi);
    accuracy=4.8;
    output= (degree>40+accuracy & degree<50-accuracy)|(degree<-40-accuracy & degree>-50+accuracy)|(degree>130+accuracy & degree<140-accuracy)|(degree<-130-accuracy & degree>-140+accuracy);
end

function [x,y] =createSquare(x1,y1,x2,y2)

    x = [x1, x1, x2, x2, x1];
    y = [y1, y2, y2, y1, y1];
    
end

function [newCurve,newCurveR] = splitCurves(oldPlot,oldPlotR,yLimit,select)
    %graphP,graphR prokeitai gia ta arxeia poy 8a kanoyn return 2 plots
    %panw h katw apo to limit.Ena limit mou arkei
    n=length(oldPlot);
    upCounter=1;
    upCounterR=1;
    downCounter=1;
    downCounterR=1;
    for i =1:n
        if (oldPlot(i,2)>=yLimit)   %belong to Up Plot
            newUpPlot(upCounter,1)=oldPlot(i,1);
            newUpPlot(upCounter,2)=oldPlot(i,2);
            upCounter=upCounter+1;
        else %belong to Down Plot
            newDownPlot(downCounter,1)=oldPlot(i,1);
            newDownPlot(downCounter,2)=oldPlot(i,2);
            downCounter=downCounter+1;
        end;
        if oldPlotR(i,2)>=yLimit   %belong to Up Plot
            newUpPlotR(upCounterR,1)=oldPlotR(i,1);
            newUpPlotR(upCounterR,2)=oldPlotR(i,2);
            upCounterR=upCounterR+1;
        else
            newDownPlotR(downCounterR,1)=oldPlotR(i,1);
            newDownPlotR(downCounterR,2)=oldPlotR(i,2);
            downCounterR=downCounterR+1;
        end
    end;
    %------Add one extra element for Up Plots so they will have meeting
    %point
    %for Normal plot we have to add Point At the End
    newUpPlot(upCounter,1)=oldPlot(upCounter,1);
    newUpPlot(upCounter,2)=oldPlot(upCounter,2);
    %for Rotate plot we have to add Point At the Begin
    indexToAdd=n+1-upCounterR;
    newUpPlotR=[oldPlotR(indexToAdd,1) oldPlotR(indexToAdd,2);newUpPlotR];
    %Add one extra element for Down Plots so they will have meeting point
    %for Normal plot we have to add Point At the Begin
    indexToAdd=n+1-downCounter;
    newDownPlot=[oldPlot(indexToAdd,1) oldPlot(indexToAdd,2);newDownPlot];
    %for Rotate plot  we have to add Point At the End
    newDownPlotR(downCounterR,1)=oldPlotR(downCounterR,1);
    newDownPlotR(downCounterR,2)=oldPlotR(downCounterR,2);
    
    if strcmp(select,'up')
         newCurve=newUpPlot;
         newCurveR=newUpPlotR;
    elseif strcmp(select,'down')
        newCurve=newDownPlot;
        newCurveR=newDownPlotR;
    end
    
    
end
    

function [plotR] = rotateAtXAxis(plotR)
    n=length(plotR);
    for i=1:n
        temp=plotR(:, 1);
        plotR(:, 1)=plotR(:, 2);
        plotR(:, 2)=temp;
    end
end

function P = InterX(L1,varargin)
%INTERX Intersection of curves
%   P = INTERX(L1,L2) returns the intersection points of two curves L1 
%   and L2. The curves L1,L2 can be either closed or open and are described
%   by two-row-matrices, where each row contains its x- and y- coordinates.
%   The intersection of groups of curves (e.g. contour lines, multiply 
%   connected regions etc) can also be computed by separating them with a
%   column of NaNs as for example
%
%         L  = [x11 x12 x13 ... NaN x21 x22 x23 ...;
%               y11 y12 y13 ... NaN y21 y22 y23 ...]
%
%   P has the same structure as L1 and L2, and its rows correspond to the
%   x- and y- coordinates of the intersection points of L1 and L2. If no
%   intersections are found, the returned P is empty.
%
%   P = INTERX(L1) returns the self-intersection points of L1. To keep
%   the code simple, the points at which the curve is tangent to itself are
%   not included. P = INTERX(L1,L1) returns all the points of the curve 
%   together with any self-intersection points.
%   
%   Example:
%       t = linspace(0,2*pi);
%       r1 = sin(4*t)+2;  x1 = r1.*cos(t); y1 = r1.*sin(t);
%       r2 = sin(8*t)+2;  x2 = r2.*cos(t); y2 = r2.*sin(t);
%       P = InterX([x1;y1],[x2;y2]);
%       plot(x1,y1,x2,y2,P(1,:),P(2,:),'ro')

%   Author : NS
%   Version: 3.0, 21 Sept. 2010

%   Two words about the algorithm: Most of the code is self-explanatory.
%   The only trick lies in the calculation of C1 and C2. To be brief, this
%   is essentially the two-dimensional analog of the condition that needs
%   to be satisfied by a function F(x) that has a zero in the interval
%   [a,b], namely
%           F(a)*F(b) <= 0
%   C1 and C2 exactly do this for each segment of curves 1 and 2
%   respectively. If this condition is satisfied simultaneously for two
%   segments then we know that they will cross at some point. 
%   Each factor of the 'C' arrays is essentially a matrix containing 
%   the numerators of the signed distances between points of one curve
%   and line segments of the other.

    %...Argument checks and assignment of L2
    narginchk(1,2);
    if nargin == 1
        L2 = L1;    hF = @lt;   %...Avoid the inclusion of common points
    else
        L2 = varargin{1}; hF = @le;
    end
       
    %...Preliminary stuff
    x1  = L1(1,:)';  x2 = L2(1,:);
    y1  = L1(2,:)';  y2 = L2(2,:);
    dx1 = diff(x1); dy1 = diff(y1);
    dx2 = diff(x2); dy2 = diff(y2);
    
    %...Determine 'signed distances'   
    S1 = dx1.*y1(1:end-1) - dy1.*x1(1:end-1);
    S2 = dx2.*y2(1:end-1) - dy2.*x2(1:end-1);
    
    C1 = feval(hF,D(bsxfun(@times,dx1,y2)-bsxfun(@times,dy1,x2),S1),0);
    C2 = feval(hF,D((bsxfun(@times,y1,dx2)-bsxfun(@times,x1,dy2))',S2'),0)';

    %...Obtain the segments where an intersection is expected
    [i,j] = find(C1 & C2); 
    if isempty(i),P = zeros(2,0);return; end
    
    %...Transpose and prepare for output
    i=i'; dx2=dx2'; dy2=dy2'; S2 = S2';
    L = dy2(j).*dx1(i) - dy1(i).*dx2(j);
    i = i(L~=0); j=j(L~=0); L=L(L~=0);  %...Avoid divisions by 0
    
    %...Solve system of eqs to get the common points
    P = unique([dx2(j).*S1(i) - dx1(i).*S2(j), ...
                dy2(j).*S1(i) - dy1(i).*S2(j)]./[L L],'rows')';
              
    function u = D(x,y)
        u = bsxfun(@minus,x(:,1:end-1),y).*bsxfun(@minus,x(:,2:end),y);
    end
end

