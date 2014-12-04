function [DISPARITY_X, DISPARITY_Y, NPeak] = regionDisparity(REGION1, REGION2, CONFIDENCE_INTERVAL)
%THIS FUNCTION CALCULATES THE DISPARITY BETWEEN TWO CLOSELY MATCHING
%IMAGES(LIKE AFTER DWO OR DEFORM IN THE SECOND PASS WHEN PARTICLES IN IMAGE
%PAIR IS CLOSE TO EACH OTHER). 

%INPUT:
%REGION1 = IMAGE1, WITHOUT ZERO MEAN OPERATION IN PRANA AND PREFERRABLY SPATIAL FILTER
%APPLIED.
%REGION2 = IMAGE2, WITHOUT ZERO MEAN OPERATION IN PRANA AND PREFERRABLY SPATIAL FILTER
%APPLIED.
% OUTPUT:-
% DELTAX=UNCERTAINTY IN X USING DISPARITY DISTRIBUTION DISPX
% DELTAY=UNCERTAINTY IN Y USING DISPARITY DISTRIBUTION DISPY

%THIS IS BASED ON THE PRINCIPLE DISCUSSED IN  "PIV uncertainty
%quantification by image matching", MST. 2013.BY Andrea Sciacchitano1, BernhardWieneke2 and Fulvio Scarano1

%TWO MAIN THINGS TO VARY ARE "sr" AND "thresh" DEPENDING ON IMAGES

%SEARCH RADIUS IN PIXELS WITHIN WHICH INDIVIDUAL IMAGE PARTICLE PEAK IS EXPECTED TO LIE FROM THE CENTER OF THE PRODUCT OF THE INTENSITIES OF IMAGE PAIRS
sr = 1;

improduct=((REGION1).*(REGION2)); %PRODUCT OF INTENSITIES
 
rp3 = improduct(improduct~=min(improduct(:)));%GETTING RID OF MINIMUM
 
thresh= 1 * rms(rp3(:));%SELECTING THRESHOLD
improduct(improduct<thresh)=0;%THRESHOLDING IMAGE
peakmat=imregionalmax(improduct,8);
  
cnt=peakmat(peakmat==1);
NPeak=length(cnt);%NUMBER OF PEAKS DETECTED

 k=1;

while NPeak<=1 && k<=8
    improduct=((REGION1).*(REGION2));
    rp3=improduct(improduct~=min(improduct(:)));
    
    thresh=(1/2^k)*rms(rp3(:));%IF NO PEAK FOUND HALF THE THRESHOLD
    
    improduct(improduct<thresh)=0;
    peakmat=imregionalmax(improduct,8);
    cnt=peakmat(peakmat==1);
    NPeak=length(cnt);
    
    k=k+1;

end

%This offset is added such that gaussian fit does not fail for neighbouring
%0 or negative intensity (as mean subtraction is there);
REGION1 = 0 + REGION1-min(min(REGION1));
REGION2 = 0 + REGION2-min(min(REGION2));

%peakmat(peakmat(:)~=0)=1;
[Indylist,Indxlist] =find(peakmat(:,:)==1);
improduct2=diag(improduct(Indylist,Indxlist));
cw=(improduct2).^(0.5);%WEIGHTING FUNCTION
dispx=zeros(NPeak,1);
dispy=zeros(NPeak,1);

[Ny,Nx]=size(REGION1);
lflag='false';

%FOR EACH IMAGE PRODUCT INTENSITY PEAK
for ji = 1:NPeak
    %6590
    Indy=Indylist(ji);
    Indx=Indxlist(ji);
    
    %IF CORNER POINTS SELECTED THEY ARE REMOVED TO AVOID TROUBLE IN
    %GAUSSIAN FITTING
    if Indy<=3 || Indy>=(Ny-3) || Indx<=3 || Indx>=(Nx-3)
        dispx(ji)=-100;
        dispy(ji)=-100;
        lflag='true';
        continue;
    end
    %FORM THE WINDOW ON WHICH TO SEARCH    
    sradiusy=Indy-sr:1:Indy+sr;
    sradiusx=Indx-sr:1:Indx+sr;
    %SELECT THE REGION FROM INPUT IMAGES WHICH CORRESPONDS TO THE WINDOW
    rA=REGION1(sradiusy,sradiusx);
    rB=REGION2(sradiusy,sradiusx);
    
%PEAKS IN INDIVIDUAL IMAGE CORRESPONDING TO SEARCH WINDOW
    [im1y,im1x] =find(rA==max(rA(:)));
    [im2y,im2x] =find(rB==max(rB(:)));
 
    %iF MULTIPLE PEAKS FOUND SELECT THE PEAK CLOSEST TO WINDOW CENTER
    if length(im1y)>1 && length(im1x)>1
        %keyboard;
        d1=((im1y-(sr+1)).^2+(im1x-(sr+1)).^2).^0.5;
        in1=find(d1==min(d1));
        in1=in1(1);
        im1y=im1y(in1);
        im1x=im1x(in1);

    end
    if length(im2y)>1 && length(im2x)>1
        
        d2=((im2y-(sr+1)).^2+(im2x-(sr+1)).^2).^0.5;
        in2=find(d2==min(d2));
        in2=in2(1);
        im2y=im2y(in2);
        im2x=im2x(in2);
    end    
        
    %GETTING PEAK COORDINATE W.R.T. IMAGE COORDINATE SYSTEM
    im1y=Indy+im1y-(sr+1);
    im1x=Indx+im1x-(sr+1);
    im2y=Indy+im2y-(sr+1);
    im2x=Indx+im2x-(sr+1);
    
    %DOING SUBPIXEL FIT FOR EACH PEAK EACH DIRECTION
    lCm11 = log(REGION1(im1y-1,im1x));
    lC001 = log(REGION1(im1y,im1x));
    lCp11 = log(REGION1(im1y+1,im1x));
    if (2*(lCm11+lCp11-2*lC001)) == 0 % if three point Gaussian fit fail, set the particle diameter to a deafault guess value
        shift_erry1 = 0;
    else
        shift_erry1 = (lCm11-lCp11)/(2*(lCm11+lCp11-2*lC001));
    end
    clear lCm11 lC001 lCp11;
    
    lCm11 = log(REGION1(im1y,im1x-1));
    lC001 = log(REGION1(im1y,im1x));
    lCp11 = log(REGION1(im1y,im1x+1));
    if (2*(lCm11+lCp11-2*lC001)) == 0 % if three point Gaussian fit fail, set the particle diameter to a deafault guess value
        shift_errx1 = 0;
    else
        shift_errx1 = (lCm11-lCp11)/(2*(lCm11+lCp11-2*lC001));
    end
    clear lCm11 lC001 lCp11;
    
    lCm11 = log(REGION2(im2y-1,im2x));
    lC001 = log(REGION2(im2y,im2x));
    lCp11 = log(REGION2(im2y+1,im2x));
    if (2*(lCm11+lCp11-2*lC001)) == 0 % if three point Gaussian fit fail, set the particle diameter to a deafault guess value
        shift_erry2 = 0;
    else
        shift_erry2 = (lCm11-lCp11)/(2*(lCm11+lCp11-2*lC001));
    end
    clear lCm11 lC001 lCp11;
    
    lCm11 = log(REGION2(im2y,im2x-1));
    lC001 = log(REGION2(im2y,im2x));
    lCp11 = log(REGION2(im2y,im2x+1));
    if (2*(lCm11+lCp11-2*lC001)) == 0 % if three point Gaussian fit fail, set the particle diameter to a deafault guess value
        shift_errx2 = 0;
    else
        shift_errx2 = (lCm11-lCp11)/(2*(lCm11+lCp11-2*lC001));
    end
    clear lCm11 lC001 lCp11;
    
    if abs(shift_errx1)>1 || isnan(shift_errx1)
        shift_errx1=0;
    end
    if abs(shift_errx2)>1 || isnan(shift_errx2)
        shift_errx2=0;
    end
    if abs(shift_erry1)>1 || isnan(shift_erry1)
        shift_erry1=0;
    end
    if abs(shift_erry2)>1 || isnan(shift_erry2)
        shift_erry2=0;
    end
    if isnan(shift_errx1)
        shift_errx1=0;
    end
    if isnan(shift_errx2)
        shift_errx2=0;
    end
    if isnan(shift_erry1)
        shift_erry1=0;
    end
    if isnan(shift_erry2)
        shift_erry2=0;
    end
    
    %PARTICLE CENTER LOCATIONS UPTO SUBPIXEL APPROXIMATION
    Indy1=im1y+shift_erry1;
    Indx1=im1x+shift_errx1;
    Indy2=im2y+shift_erry2;
    Indx2=im2x+shift_errx2;
    
    %DISPARITY BETWEEN PARTICLE PAIR FOR EACH PAIR OF MATCHED PARTICLE
    dispx(ji)=Indx2-Indx1;
    dispy(ji)=Indy2-Indy1;
    
end

%DISCARDING POINTS WHICH WERE CLOSE TO THE BOUNDARY
if strcmp(lflag,'true')
    %keyboard;
    ind=find(dispx(:)==-100);
    dispx(ind)=[];
    dispy(ind)=[];
    cw(ind)=[];
    NPeak=NPeak-length(ind);

end

%CALCULATING WEIGHTED MEAN AND VARIANCE OF DISPARITY DISTRIBUTION FOR THIS
%WINDOW PAIRS AND SAVING THE UNCERTAINTY AS OUTPUT.
mewx=(1/sum(cw))*sum(cw.*dispx);
sigx=(1/sum(cw))*sum(cw.*((dispx-mewx).^2));

t_inv = tinv(0.5 + CONFIDENCE_INTERVAL/2, NPeak-1);

DISPARITY_X = t_inv * sqrt(mewx^2 + (sigx/sqrt(NPeak))^2);

mewy=(1/sum(cw))*sum(cw.*dispy);
sigy=(1/sum(cw))*sum(cw.*((dispy-mewy).^2));

DISPARITY_Y = t_inv * sqrt(mewy^2 + (sigy/sqrt(NPeak))^2);

if isnan(DISPARITY_X) || isnan(DISPARITY_Y) 
    DISPARITY_X = nan;
    DISPARITY_Y = nan;
end


end