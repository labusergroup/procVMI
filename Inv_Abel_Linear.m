% This function inverts the raw VMI (velocity map imaging) images.
% Input: center of the image, raw images/coordinate file. Image MUST have the x-axis as the symmetry axis.
% Output: 3D density distribution and integrated yield.

% ref: The Journal of Chemical Physics 147, 013922 (2017); doi: http://dx.doi.org/10.1063/1.4981917

% NB: Only process image input atm.

%% Variables and definitions:

% x is the laser polarizatoin/symmetry axis. z is the projection axis.
% Ring: an array of 3D distribution densities (panel a in Fig.4 in the ref), labeled by {y,x} tuples.
% AngIntegrated: angularly integrated yields (panel c).

function [Ring, AngIntegrated] = Inv_Abel_Linear(Centre, varargin)

%% (1) Checke the center (pixel- or lattice site-centered). Check if the input is a coordinate list (EXACTLY 2 coloums) or an image (otherwise).

% Determeine whether pixel- or lattice-centered
CentralRow = Centre(1);
CentralCol = Centre(2);
if ceil(CentralRow) - CentralRow >= 0.25 && CentralRow - floor(CentralRow) >= 0.25 && ceil(CentralCol) - CentralCol >= 0.25 && CentralCol - floor(CentralCol) >= 0.25
    Lattice_Centered = 1; %lattice-centered
    % lower-right pixel of true centre in case of lattice-centered
    yc = ceil(CentralRow);
    xc = ceil(CentralCol);
else
    Lattice_Centered = 0; %pixel-centered
    yc = round(CentralRow);
    xc = round(CentralCol);
end

% % Determine whether image or coordinates, and then fold the image
% if size(varargin(1),2) == 2
%     isImg = false;
% else
%     isImg = true;
% end

isImg = 1;
if isImg  % image input
    Image = varargin{1};
    ImageOrig = Image;
    DimOrig = size(Image);
    %     Dim = 360;   % size of AOI = 360 pixel
    RadiusMax = min(yc-1, DimOrig(1)-yc);   % max y < Dim/2
    XRange = min(xc-1,DimOrig(2)-xc);   % < max y, cut edge
    
    % To average Image, fold to one quadrant/half. The inversion only operates on the forth quadrant
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
%                         Image(yc+y,xc+x) = (Image(yc+y,xc+x)+Image(yc+y,xc-x))/2; %only lower half
%                         Image(yc+y,xc+x) = (Image(yc+y,xc+x)+Image(yc-y,xc+x))/2; %only right half

                        % all quadrants
                        Image(yc+y,xc+x) = (Image(yc+y,xc+x)+Image(yc-y,xc+x)+Image(yc+y,xc-x)+Image(yc-y,xc-x))/4; %all quadrants
                        
                        % only left half
%                         Image(yc+y,xc+x) = (Image(yc+y,xc-x)+Image(yc-y,xc-x))/2; %only lower half
%                         Image(yc+y,xc+x) = (Image(yc+y,xc-x)+Image(yc-y,xc+x))/2; %only right half

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
 
%% (2) Construct Abel inversion matrix, cf Eq. (2.3) in Arthur's thesis.
% Set up a series of transformation matrices, containing weights distributions.
% Each matrix depends on its size or number of rings, 'NumRing', in a slice (x = constant).
% Assume a circular AOI (area of interest).
% Given x and AOI (which is defined by 'RadiusMax'), 'NumRing'
% is determined by NumRing =(RadiusMax^2-x^2)^0.5. And so is the linear maps
% between projected pixel values, 'Proj(y,x)', and ring densities, 'Ring(y,x)'.
% Proj(NumRing,y) = Area(NumRing,y,k)*Ring(NumRing,k), summed over k>=y.
% Ring(NumRing,y) = Weight(NumRing,y,k)*Proj(NumRing,k), summed over k<=y.
% 'Area' and 'Weight' are pre-defined matrix-valued functions of NumRing st. Area*Weight = 1.

% Define 'Area' and its inverse 'Weight' in one structure array 'Distr':
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
%% (3) Invert the image

% Image input
if isImg == 1  
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


% coord input, NOT IN USE
if isImg == 0
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

end