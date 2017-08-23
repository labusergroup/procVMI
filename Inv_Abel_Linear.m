% modified on 2015-05-19,   
% input central coord and image/text, output density distribution and
% integrated yield.
function [Ring, AngIntegrated] = Inv_Abel_Linear(Centre, varargin)
% Ring: array of densities
% AngIntegrated: angularly integrated yields

% ref. page 143 and 144 on Arthur's notebook I. z is the projection axis.
% In this code, x is laser polarizatoin/symmetry axis.Rings are labeled by {y,x} tuples.

CentralRow = Centre(1);
CentralColoum = Centre(2);
if ceil(CentralRow) - CentralRow >= 0.25 && CentralRow - floor(CentralRow) >= 0.25 && ceil(CentralColoum) - CentralColoum >= 0.25 && CentralColoum - floor(CentralColoum) >= 0.25
    Lattice_Centered = 1 %lattice-centered
    % lower-right pixel of true centre in case of lattice-centered
    yc = ceil(CentralRow);
    xc = ceil(CentralColoum);
else
    Lattice_Centered = 0; %pixel-centered
    yc = round(CentralRow);
    xc = round(CentralColoum);
end


% tic
% to invert image or coordinates depending the input
% try
%     isjpg = strcmp(imfinfo(varargin{1}).Format,'jpg');
% catch
%     isjpg = false;
% end
% 
isjpg = 1;
if isjpg % image input
    Image = varargin{1};
    ImageOrig = Image;
    DimOrig = size(Image);
    %     Dim = 360;   % size of AOI = 360 pixel
    RadiusMax = min(yc-1, DimOrig(1)-yc);   % max y < Dim/2
    XRange = min(xc-1,DimOrig(2)-xc);   % < max y, cut edge
    
    %   to average Image, fold to one quadrant/half. The inversion only
    %   deals with the forth quadrant
    FlagAvg = 1;
    if FlagAvg == 1
    if Lattice_Centered == 1
        if XRange == 0 % 1-D case
            for y = 0:RadiusMax
                Image(yc+y,1) = (Image(yc+y,1)+Image(yc-1-y,1))/2;
            end
        else
            for x=0:XRange %quadrant size used
                for y = 0:RadiusMax
                    Image(yc+y,xc+x) = (Image(yc+y,xc+x)+Image(yc-1-y,xc+x)+Image(yc+y,xc-1-x)+Image(yc-1-y,xc-1-x))/4;
                end
            end
        end
    else % pixel-centered
        if XRange == 0 % 1-D case
            for y = 0:RadiusMax
                Image(yc+y,1) = (Image(yc+y,1)+Image(yc-y,1))/2;
            end
        else
            for x=1:XRange %quadrant size used
                for y = 1:RadiusMax
                    % only lower half
%                     Image(yc+y,xc+x) = (Image(yc+y,xc+x)+Image(yc+y,xc-x))/2; %only lower half
%                     Image(yc+y,xc+x) = (Image(yc+y,xc+x)+Image(yc-y,xc+x))/2; %only right half
                    % all quadrants
                    Image(yc+y,xc+x) = (Image(yc+y,xc+x)+Image(yc-y,xc+x)+Image(yc+y,xc-x)+Image(yc-y,xc-x))/4; %all quadrants
                    % only left half
%                     Image(yc+y,xc+x) = (Image(yc+y,xc-x)+Image(yc-y,xc-x))/2; %only lower half
%                     Image(yc+y,xc+x) = (Image(yc+y,xc-x)+Image(yc-y,xc+x))/2; %only right half
                    
                    
                    Image(yc,xc+x) = (Image(yc,xc+x)+Image(yc,xc-x))/2;
                    Image(yc+y,xc) = (Image(yc+y,xc)+Image(yc-y,xc))/2;
                end
            end
        end
    end
    end
else
    
    %     coord input, not in use
    % fragment coordinate .txt file
    %     Frag = uigetfile('*.txt');
%     Coord = load(varargin{1});
end

flag_plot = 0; % set to 1 to make plots, otherwise no plots
 

%%
% set up a series of transformation matrices, containing weights distributions, which only depends on
% number of rings, 'NumRing', in a slice (x = constant).
% Assume a circular AOI (area of interest).
% Given x and AOI (which is defined by 'RadiusMax'), 'NumRing'
% is determined by NumRing =(RadiusMax^2-x^2)^0.5. And so is the linear maps
% between projected pixel values, 'Proj(y,x)', and ring densities, 'Ring(y,x)'.
% Proj(NumRing,y) = Area(NumRing,y,k)*Ring(NumRing,k), summed over k>=y.
% Ring(NumRing,y) = Weight(NumRing,y,k)*Proj(NumRing,k), summed over k<=y.
% 'Area' and 'Weight' are pre-define matrix-valued functions of NumRing st. Area*Weight = 1.

% Here to define 'Area' and its inverse 'Weight' in one structure array 'Distr':
SizeFactor = 1; %refinement for better resolution
Size = RadiusMax*SizeFactor; %size of 'Distr' 
Distr = struct('Area',{},'Weight',{});
    
if Lattice_Centered == 1    % lattice-centered
    for NumRing = 1:Size
        Distr(NumRing).Area = zeros(NumRing,NumRing);
        Distr(NumRing).Weight = zeros(NumRing,NumRing);
        for y = 1:NumRing
            Sum = 0;
            for k = y:NumRing
                Distr(NumRing).Area(y,k) = k^2*(acos((y-1)/k)-acos(y/k))+y*(k^2-y^2)^.5-(y-1)*(k^2-(y-1)^2)^.5-Sum;
                Sum = Sum + Distr(NumRing).Area(y,k);
            end
        end
        Distr(NumRing).Weight = inv(Distr(NumRing).Area);
    end
end
if Lattice_Centered == 0     % pixel-centered
    for NumRing = 1:Size
        Distr(NumRing).Area = zeros(NumRing,NumRing);
        Distr(NumRing).Weight = zeros(NumRing,NumRing);
        for y = 1:NumRing
            Sum = 0;
            for k = y:NumRing
                Distr(NumRing).Area(y,k) = (k-.5)^2*(acos((y-1.5)/(k-.5))-acos((y-.5)/(k-.5)))+(y-.5)*((k-.5)^2-(y-.5)^2)^.5-(y-1.5)*((k-.5)^2-(y-1.5)^2)^.5-Sum;
                Sum = Sum + Distr(NumRing).Area(y,k);
            end
        end
%         Distr(NumRing).Area(1,:) = Distr(NumRing).Area(1,:)/2; %centre ring
        Distr(NumRing).Weight = inv(Distr(NumRing).Area);
    end
end
%%
% inversion options:image input or coordinate input

% Image input
if isjpg == 1  
    %load Image, choose the appropriate distribution function and construct 3-d Image.
    Ring = zeros(DimOrig(1),DimOrig(2));
    AngIntegrated = zeros(DimOrig(1),DimOrig(2));
    
    %lattice-centered
    if Lattice_Centered == 1    
        %         Ring = zeros(yc+RadiusMax,xc+XRange);
        for x = 0:XRange %x and y are distance here
            NumRing = RadiusMax - 1;
            Ring(yc:(yc+NumRing-1),x+xc) = Distr(NumRing).Weight*Image(yc:(yc+NumRing-1),x+xc);
            %             Ring(yc:(yc+NumRing-1),x+xc) = Dist_Weight(NumRing).Weight*Image(yc:(yc+NumRing-1),x+xc);
            %pi-angular integral and cross-section
            for y = 0:NumRing
                AngIntegrated(y+yc,x+xc) = Ring(y+yc,x+xc)*pi/2*((y+1)^2-y^2);
                AngIntegrated(-y+yc-1,x+xc) = AngIntegrated(y+yc,x+xc);
                AngIntegrated(y+yc,-x+xc-1) = AngIntegrated(y+yc,x+xc);
                AngIntegrated(-y+yc-1,-x+xc-1) = AngIntegrated(y+yc,x+xc);
            end
        end
    end
    
    %pixel-centered
    if Lattice_Centered == 0         
        for x = 0:XRange %x and y are distance here
            NumRing = RadiusMax - 1;
            Ring(yc:(yc+NumRing-1),x+xc) = Distr(NumRing).Weight * Image(yc:(yc+NumRing-1),x+xc);
            %             Ring(yc:(yc+NumRing-1),x+xc) = Dist_Weight(NumRing).Weight*Image(yc:(yc+NumRing-1),x+xc);
%             Ring(yc,x+xc) = Ring(yc,x+xc)/2; %centre rEing
            % pi-angular integral and cross-section
            AngIntegrated(yc,x+xc) = Ring(yc,x+xc)*pi*.5^2;
            AngIntegrated(yc,-x+xc) = AngIntegrated(yc,x+xc);   
            Ring(yc,-x+xc) = Ring(yc,x+xc);
            for y = 1:NumRing
                AngIntegrated(y+yc,x+xc) = Ring(y+yc,x+xc)*pi/2*((y+.5)^2-(y-.5)^2);
                AngIntegrated(-y+yc,x+xc) = AngIntegrated(y+yc,x+xc);
                AngIntegrated(y+yc,-x+xc) = AngIntegrated(y+yc,x+xc);
                AngIntegrated(-y+yc,-x+xc) = AngIntegrated(y+yc,x+xc);

                Ring(-y+yc,x+xc) = Ring(y+yc,x+xc);
                Ring(y+yc,-x+xc) = Ring(y+yc,x+xc);
                Ring(-y+yc,-x+xc) = Ring(y+yc,x+xc);
            end
        end
    end
end


% coord input
if isjpg == 0
    Dim = DimOrig*SizeFactor;
    XRange = XRange*SizeFactor; %x max
    RadiusMax = RadiusMax*SizeFactor; %y max
    Ring = zeros(Dim(1),Dim(2));
    AngIntegrated = zeros(Dim(1),Dim(1));
    xc = xc*SizeFactor;
    yc = yc*SizeFactor;
    Coord = Coord*SizeFactor;
    
    %lattice centered
    if Lattice_Centered == 1 
        Ring = zeros(Dim(1),Dim(2));
        Image = zeros(Dim(1),Dim(2));
        for j = 1:size(Coord,1)
            if Coord(j,1) < yc-0.5 
                Coord(j,1) = yc + (yc - Coord(j,1)) - 1;
            end
            if Coord(j,2) < xc-0.5 
                Coord(j,2) = xc + (xc - Coord(j,2)) - 1;
            end
            if Coord(j,1) < Dim(1) && Coord(j,2) < Dim(2)
                Image(floor(Coord(j,1)),floor(Coord(j,2))) = Image(floor(Coord(j,1)),floor(Coord(j,2))) + 1;
            end
        end
        %         Ring = zeros(yc+RadiusMax,xc+XRange);
        for x = 0:XRange %x and y are distance here
            NumRing = RadiusMax - 1;
            Ring(yc:(yc+NumRing-1),x+xc) = Distr(NumRing).Weight*Image(yc:(yc+NumRing-1),x+xc);
            %             Ring(yc:(yc+NumRing-1),x+xc) = Dist_Weight(NumRing).Weight*Image(yc:(yc+NumRing-1),x+xc);
            %pi-angular integral and cross-section
            for y = 0:NumRing-1
                AngIntegrated(y+yc,x+xc) = Ring(y+yc,x+xc)*pi/2*((y+1)^2-y^2);
                AngIntegrated(-y+yc-1,x+xc) = AngIntegrated(y+yc,x+xc);
                AngIntegrated(y+yc,-x+xc-1) = AngIntegrated(y+yc,x+xc);
                AngIntegrated(-y+yc-1,-x+xc-1) = AngIntegrated(y+yc,x+xc);
            end
        end
    end
    %pixel-centered
    if Lattice_Centered == 0   
                Ring = zeros(Dim(1),Dim(2));
        Image = zeros(Dim(1),Dim(2));
        for j = 1:size(Coord,1)
            if Coord(j,1) < yc 
                Coord(j,1) = yc + (yc - Coord(j,1));
            end
            if Coord(j,2) < xc
                Coord(j,2) = xc + (xc - Coord(j,2));
            end
            if Coord(j,1) < Dim(1) && Coord(j,2) < Dim(2)
                Image(floor(Coord(j,1)),floor(Coord(j,2))) = Image(floor(Coord(j,1)),floor(Coord(j,2))) + 1;
            end
        end
        for x = 0:XRange %x and y are distance here
            NumRing = RadiusMax - 1;
            Ring(yc:(yc+NumRing-1),x+xc) = Distr(NumRing).Weight*Image(yc:(yc+NumRing-1),x+xc);
            %             Ring(yc:(yc+NumRing-1),x+xc) = Dist_Weight(NumRing).Weight*Image(yc:(yc+NumRing-1),x+xc);
%             Ring(yc,x+xc) = Ring(yc,x+xc)/2; %centre ring
            % pi-angular integral and cross-section
            AngIntegrated(yc,x+xc) = Ring(yc,x+xc)*pi*.5^2;
            AngIntegrated(yc,-x+xc) = AngIntegrated(yc,x+xc);
            
            
            
            for y = 1:NumRing
                AngIntegrated(y+yc,x+xc) = Ring(y+yc,x+xc)*pi/2*((y+.5)^2-(y-.5)^2);
                AngIntegrated(-y+yc,x+xc) = AngIntegrated(y+yc,x+xc);
                AngIntegrated(y+yc,-x+xc) = AngIntegrated(y+yc,x+xc);
                AngIntegrated(-y+yc,-x+xc) = AngIntegrated(y+yc,x+xc);
            end
        end
    end
end

%%
% plots

% PtoE = 0.0887;
% PtoE = 0.0074/600*800; % closer focus, since 2014-09-19, 0.0074 for -600/-420
if flag_plot == 1
    figure,imagesc(AngIntegrated); title('cross section');
    figure,imagesc(Ring); title('Ring');
    
    %%
    % lineout
    lnt(:,:) = radial_lineoutVer3P(AngIntegrated,[yc xc],10);
    figure
    % plot(PtoE*(1:x).^2, squeeze(lnt(3,:)+lnt(8,:)));
    x=size(lnt(3,:),2);
    plot(PtoE*((1:x)-0.5).^2, squeeze(lnt(3,1:x)+lnt(8,1:x))/(2*(PtoE))./((1:x)-0.5)/max((squeeze(lnt(3,50:x)+lnt(8,50:x)))/(2*(PtoE))./((50:x)-0.5)));
%     
%     
%     
%     % project and compare difference
%     Proj = zeros(DimOrig(1),DimOrig(2));
%     if Lattice_Centered == 1
%         %lattice centered
%         for x = 0:XRange %x and y are distance here
%             NumRing = round((RadiusMax^2-x^2)^.5);
%             Proj(yc:(yc+NumRing-1),x+xc) = Distr(NumRing).Area*Ring(yc:(yc+NumRing-1),x+xc);
%             Proj(yc-1:-1:(yc-NumRing),x+xc) = Proj(yc:(yc+NumRing-1),x+xc);
%             Proj((yc-NumRing):(yc+NumRing-1),-x+xc-1) = Proj((yc-NumRing):(yc+NumRing-1),x+xc);
%         end
%     end
%     if Lattice_Centered == 0
%         %pixel-centered
%         for x = 0:XRange %x and y are distance here
%             NumRing = round((RadiusMax^2-x^2)^.5);
%             Proj(yc:(yc+NumRing-1),x+xc) = Distr(NumRing).Area*Ring(yc:(yc+NumRing-1),x+xc);
%             Proj(yc:-1:(yc-NumRing+1),x+xc) = Proj(yc:(yc+NumRing-1),x+xc);
%             Proj((yc-NumRing+1):(yc+NumRing-1),-x+xc) = Proj((yc-NumRing+1):(yc+NumRing-1),x+xc);
%         end
%     end
%     if isjpg == 1
%         Diff = ImageOrig - Proj;
%         figure
%         imagesc(Diff);
%         title('difference');
%     end
%     figure
%     imagesc(Proj); title('projection')
% end
% 
% % if flag_plot == 1
% %     lnt(:,:) = radial_lineoutVer3P(AngIntegrated,[yc xc],10);
% %     plot(PtoE*(1:x).^2, squeeze(lnt(3,:)+lnt(8,:)));
% %     x=size(lnt(3,:),2);
% %     plot(PtoE*((1:x)-0.5).^2, squeeze(lnt(3,1:x)+lnt(8,1:x))/(2*(PtoE))./((1:x)-0.5)/max((squeeze(lnt(3,50:x)+lnt(8,50:x)))/(2*(PtoE))./((50:x)-0.5)));
% % Dominik's
% %     lnt(:,:) = radial_lineoutVer3P(AngIntegrated,[yc xc],1);
% % end
% % toc
end