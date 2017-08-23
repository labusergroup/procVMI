% update 2015-5-19: only handle integral Centre
% update 2014-3-19: boundary pixel is dealt, both radially and angularly. But one still have to assume pixel-center --> next time, perhaps.

% This function takes an image and return its radial lineout.
% Input: monochromatic 2D image, its center, number of sectors starting from positive y-axis, counterclockwise.
% Output: a lineout variable histo = [num of sectors, num of pixel].

function histo = Lineout_Radial(pic,center,Nsectors)

%format for center: [y x] = [row column]
[dimy,dimx] = size(pic);
tmpx = ((1 : dimx) - center(2)).^2; %distance squared
tmpy = ((1 : dimy) - center(1)).^2;
% Rmax = ceil(max([sqrt(tmpx(1)+tmpy(1)), sqrt(tmpx(1)+tmpy(dimy)), sqrt(tmpx(dimx)+tmpy(1)), sqrt(tmpx(dimx)+tmpy(dimy))]));
Rmax = ceil(max([tmpx(1)^.5, tmpx(end)^.5,tmpy(1)^.5,tmpy(end)^.5])); %lineout length
histo = zeros(Nsectors, Rmax+1);
A = histo;  %number of hits
%sector vector
dS = 2*pi/Nsectors;
for indx = 1 : dimx
    for indy = 1 : dimy
        %distance
        %         R =  floor(sqrt(tmpx(indx) + tmpy(indy)))+1;
        R = sqrt(tmpx(indx) + tmpy(indy));
        if R <= Rmax
            R1 = floor(R);
            R2 = ceil(R);
            alpha = angle((center(1)-indy)+1i*(center(2)-indx));
            alpha = mod(alpha,2*pi);
            indS = max(1,ceil(alpha/dS)); %sector
            indS1 = indS - 1;
            indS2 = indS + 1;
            if indS == 1
                indS1 = Nsectors;
            end
            if indS == Nsectors
                indS2 = 1;
            end
            alphaRem = rem(alpha,dS); %angle to sector boundary
%             alphaD = min(alphaRem, dS-alphaRem); 
% assign values, fractional if necessary
            if R == 0 %center pixel
                histo(:,1) = pic(indy,indx)/Nsectors;
            elseif R*alphaRem <= 0.5    %if close to the previous sector
                histo(indS1,R1) = histo(indS1,R1) + (0.5-R*alphaRem)*round(R2-R)*(R2-R-0.5)*pic(indy,indx); %previous pixel in the lineout
                A(indS1,R1) = A(indS1,R1) + (0.5-R*alphaRem)*round(R2-R)*(R2-R-0.5);
                histo(indS1,R2) = histo(indS1,R2) + (0.5-R*alphaRem)*(1-abs(R-R1-0.5))*pic(indy,indx);  % current pixel
                A(indS1,R2) = A(indS1,R2) + (0.5-R*alphaRem)*(abs(R-R1-0.5)+0.5);
                histo(indS1,R2+1) = histo(indS1,R2+1) + (0.5-R*alphaRem)*round(R-R1)*(R-R1-0.5)*pic(indy,indx); %next pixel
                A(indS1,R2+1) = A(indS1,R2+1) + (0.5-R*alphaRem)*round(R-R1)*(R-R1-0.5);
                
                histo(indS,R1) = histo(indS,R1) + (0.5+R*alphaRem)*round(R2-R)*(R2-R-0.5)*pic(indy,indx);
                A(indS,R1) = A(indS,R1) + (0.5+R*alphaRem)*round(R2-R)*(R2-R-0.5);
                histo(indS,R2) = histo(indS,R2) + (0.5+R*alphaRem)*(1-abs(R-R1-0.5))*pic(indy,indx);
                A(indS,R2) = A(indS,R2) + (0.5+R*alphaRem)*(abs(R-R1-0.5)+0.5);
                histo(indS,R2+1) = histo(indS,R2+1) + (0.5+R*alphaRem)*round(R-R1)*(R-R1-0.5)*pic(indy,indx);
                A(indS,R2+1) = A(indS,R2+1) + (0.5+R*alphaRem)*round(R-R1)*(R-R1-0.5);
            elseif R*(dS-alphaRem) <= 0.5   %if close to the next sector
                histo(indS2,R1) = histo(indS2,R1) + (0.5-R*(dS-alphaRem))*round(R2-R)*(R2-R-0.5)*pic(indy,indx);
                A(indS2,R1) = A(indS2,R1) + (0.5-R*(dS-alphaRem))*round(R2-R)*(R2-R-0.5);
                histo(indS2,R2) = histo(indS2,R2) + (0.5-R*(dS-alphaRem))*(1-abs(R-R1-0.5))*pic(indy,indx);
                A(indS2,R2) = A(indS2,R2) + (0.5-R*(dS-alphaRem))*(abs(R-R1-0.5)+0.5);
                histo(indS2,R2+1) = histo(indS2,R2+1) + (0.5-R*(dS-alphaRem))*round(R-R1)*(R-R1-0.5)*pic(indy,indx);
                A(indS2,R2+1) = A(indS2,R2+1) + (0.5-R*(dS-alphaRem))*round(R-R1)*(R-R1-0.5);
                
                histo(indS,R1) = histo(indS,R1) + (0.5+R*(dS-alphaRem))*round(R2-R)*(R2-R-0.5)*pic(indy,indx);
                A(indS,R1) = A(indS,R1) + (0.5+R*(dS-alphaRem))*round(R2-R)*(R2-R-0.5);
                histo(indS,R2) = histo(indS,R2) + (0.5+R*(dS-alphaRem))*(1-abs(R-R1-0.5))*pic(indy,indx);
                A(indS,R2) = A(indS,R2) + (0.5+R*(dS-alphaRem))*(abs(R-R1-0.5)+0.5);
                histo(indS,R2+1) = histo(indS,R2+1) + (0.5+R*(dS-alphaRem))*round(R-R1)*(R-R1-0.5)*pic(indy,indx);
                A(indS,R2+1) = A(indS,R2+1) + (0.5+R*(dS-alphaRem))*round(R-R1)*(R-R1-0.5);
            else % well within a sector
                histo(indS,R1) = histo(indS,R1) + round(R2-R)*(R2-R-0.5)*pic(indy,indx);
                A(indS,R1) = A(indS,R1) + round(R2-R)*(R2-R-0.5);
                histo(indS,R2) = histo(indS,R2) + (1-abs(R-R1-0.5))*pic(indy,indx);
                A(indS,R2) = A(indS,R2) + (abs(R-R1-0.5)+0.5);
                histo(indS,R2+1) = histo(indS,R2+1) + round(R-R1)*(R-R1-0.5)*pic(indy,indx);
                A(indS,R2+1) = A(indS,R2+1) + round(R-R1)*(R-R1-0.5);
            end
        end;            
    end;
end;

for indS = 1 : Nsectors
    A(indS,:) = (A(indS,:)==0) + A(indS,:);
%     histo(indS,:) = histo(indS,:) ./A(indS,:) .*(1:Rmax) ;
end

end