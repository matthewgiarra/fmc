%#codegen
function [peak_shift_rows, peak_shift_cols, corr_max_val, corr_peak_diameter, corr_max_val_delete] = subpixel(SPATIALCORRELATION, CORRELATION_WINDOW, Method, Peakswitch, COMPILED)
% This poorly commented function was taken from PRANA

% Default to using compiled codes.
if nargin < 5
    COMPILED = 1;
end

% This is redundant and can be deleted.
% Make sure deleting this won't break other function calls furst.
corr_max_val_delete = max(SPATIALCORRELATION(:));
    
% Determine size of spatial correlation
[correlationHeight, correlationWidth] = size(SPATIALCORRELATION);

%intialize indices
cc_x = -correlationWidth/2 : correlationWidth/2 - 1;
cc_y = -correlationHeight/2:correlationHeight/2-1;

%find maximum correlation value
[corr_max_val, corr_max_val_index] = max(SPATIALCORRELATION(:));

% Use 4 standard deviations for the peak sizing (e^-2)
sigma = 4;

%if correlation empty
if corr_max_val == 0
    if Peakswitch
        peak_shift_cols=zeros(1,3);
        peak_shift_rows=zeros(1,3);
        corr_max_val=zeros(1,3);
        corr_peak_diameter=zeros(1,3);
    else
        peak_shift_cols=0; peak_shift_rows=0; corr_max_val=0; corr_peak_diameter=0; 
    end
else
    if Peakswitch
        peak_shift_cols=zeros(1,3);
        peak_shift_rows=zeros(1,3);
        corr_peak_diameter=zeros(1,3);
        
        % Choose whether or not to use compiled codes.
        if COMPILED
            %Locate peaks using imregionalmax
            peakmat = locate_correlation_peaks_mex(SPATIALCORRELATION);
        else
            peakmat = locate_correlation_peaks(SPATIALCORRELATION);
        end
        
        for i=2:3
            peakmat(peakmat==corr_max_val(i-1))=0;
            [corr_max_val(i),corr_max_val_index(i)]=max(peakmat(:));
        end
        j=length(corr_max_val);
    else
        peak_shift_cols=zeros(1,1);
        peak_shift_rows=zeros(1,1);
        corr_peak_diameter=zeros(1,1);
        j=1;    
    end
    
    for i=1:j
        method=Method;
        
        %find x and y indices
        shift_locy = 1+mod(corr_max_val_index(i)-1,correlationHeight);
        shift_locx = ceil(corr_max_val_index(i)/correlationHeight);

        shift_errx=[];
        shift_erry=[];
%         find subpixel displacement in x
        if shift_locx == 1
%             keyboard
            %boundary condition 1
            shift_errx =  SPATIALCORRELATION( shift_locy , shift_locx+1 )/corr_max_val(i); method=1;
        elseif shift_locx == correlationWidth
            %boundary condition 2
            shift_errx = -SPATIALCORRELATION( shift_locy , shift_locx-1 )/corr_max_val(i); method=1;
        elseif SPATIALCORRELATION( shift_locy , shift_locx+1 ) == 0
            %endpoint discontinuity 1
            shift_errx = -SPATIALCORRELATION( shift_locy , shift_locx-1 )/corr_max_val(i); method=1;
        elseif SPATIALCORRELATION( shift_locy , shift_locx-1 ) == 0
            %endpoint discontinuity 2
            shift_errx =  SPATIALCORRELATION( shift_locy , shift_locx+1 )/corr_max_val(i); method=1;
        end
        if shift_locy == 1
%             keyboard
            %boundary condition 1
            shift_erry = -SPATIALCORRELATION( shift_locy+1 , shift_locx )/corr_max_val(i); method=1;
        elseif shift_locy == correlationHeight
            %boundary condition 2
            shift_erry =  SPATIALCORRELATION( shift_locy-1 , shift_locx )/corr_max_val(i); method=1;
        elseif SPATIALCORRELATION( shift_locy+1 , shift_locx ) == 0
            %endpoint discontinuity 1
            shift_erry =  SPATIALCORRELATION( shift_locy-1 , shift_locx )/corr_max_val(i); method=1;
        elseif SPATIALCORRELATION( shift_locy-1 , shift_locx ) == 0
            %endpoint discontinuity 2
            shift_erry = -SPATIALCORRELATION( shift_locy+1 , shift_locx )/corr_max_val(i); method=1;
        end

        
        if method==2
            
            %%%%%%%%%%%%%%%%%%%%
            % 4-Point Gaussian %
            %%%%%%%%%%%%%%%%%%%%
            
            %Since the case where M is located at a border will default to
            %the 3-point gaussian and we don't have to deal with
            %saturation, just use 4 points in a tetris block formation:
            %
            %             *
            %            ***
            
            points=[shift_locy   shift_locx   SPATIALCORRELATION(shift_locy  ,shift_locx  );...
                    shift_locy-1 shift_locx   SPATIALCORRELATION(shift_locy-1,shift_locx  );...
                    shift_locy   shift_locx-1 SPATIALCORRELATION(shift_locy  ,shift_locx-1);...
                    shift_locy   shift_locx+1 SPATIALCORRELATION(shift_locy  ,shift_locx+1)];
                
            [~,IsortI] = sort(points(:,3),'descend');
            points = points(IsortI,:);

            x1=points(1,2); x2=points(2,2); x3=points(3,2); x4=points(4,2);
            y1=points(1,1); y2=points(2,1); y3=points(3,1); y4=points(4,1);
            a1=points(1,3); a2=points(2,3); a3=points(3,3); a4=points(4,3);

            alpha(1) = (x4^2)*(y2 - y3) + (x3^2)*(y4 - y2) + ((x2^2) + (y2 - y3)*(y2 - y4))*(y3 - y4);
            alpha(2) = (x4^2)*(y3 - y1) + (x3^2)*(y1 - y4) - ((x1^2) + (y1 - y3)*(y1 - y4))*(y3 - y4);
            alpha(3) = (x4^2)*(y1 - y2) + (x2^2)*(y4 - y1) + ((x1^2) + (y1 - y2)*(y1 - y4))*(y2 - y4);
            alpha(4) = (x3^2)*(y2 - y1) + (x2^2)*(y1 - y3) - ((x1^2) + (y1 - y2)*(y1 - y3))*(y2 - y3);

            gamma(1) = (-x3^2)*x4 + (x2^2)*(x4 - x3) + x4*((y2^2) - (y3^2)) + x3*((x4^2) - (y2^2) + (y4^2)) + x2*(( x3^2) - (x4^2) + (y3^2) - (y4^2));
            gamma(2) = ( x3^2)*x4 + (x1^2)*(x3 - x4) + x4*((y3^2) - (y1^2)) - x3*((x4^2) - (y1^2) + (y4^2)) + x1*((-x3^2) + (x4^2) - (y3^2) + (y4^2));
            gamma(3) = (-x2^2)*x4 + (x1^2)*(x4 - x2) + x4*((y1^2) - (y2^2)) + x2*((x4^2) - (y1^2) + (y4^2)) + x1*(( x2^2) - (x4^2) + (y2^2) - (y4^2));
            gamma(4) = ( x2^2)*x3 + (x1^2)*(x2 - x3) + x3*((y2^2) - (y1^2)) - x2*((x3^2) - (y1^2) + (y3^2)) + x1*((-x2^2) + (x3^2) - (y2^2) + (y3^2));

            delta(1) = x4*(y2 - y3) + x2*(y3 - y4) + x3*(y4 - y2);
            delta(2) = x4*(y3 - y1) + x3*(y1 - y4) + x1*(y4 - y3);
            delta(3) = x4*(y1 - y2) + x1*(y2 - y4) + x2*(y4 - y1);
            delta(4) = x3*(y2 - y1) + x2*(y1 - y3) + x1*(y3 - y2);

            deno = 2*(log(a1)*delta(1) + log(a2)*delta(2) + log(a3)*delta(3) + log(a4)*delta(4));

            x_centroid = (log(a1)*alpha(1) + log(a2)*alpha(2) + log(a3)*alpha(3) + log(a4)*alpha(4))/deno;
            y_centroid = (log(a1)*gamma(1) + log(a2)*gamma(2) + log(a3)*gamma(3) + log(a4)*gamma(4))/deno;
            shift_errx=x_centroid-shift_locx;
            shift_erry=y_centroid-shift_locy;
            
            betas = abs((log(a2)-log(a1))/((x2-x_centroid)^2+(y2-y_centroid)^2-(x1-x_centroid)^2-(y1-y_centroid)^2));
            corr_peak_diameter(i)=sqrt(sigma^2/(2*betas));
            
        elseif any(method==[3 4])
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            % Gaussian Least Squares %
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %Find a suitable window around the peak (5x5 preferred)
            x_min=shift_locx-2; x_max=shift_locx+2;
            y_min=shift_locy-2; y_max=shift_locy+2;
            if x_min<1
                x_min=1;
            end
            if x_max>correlationWidth
                x_max=correlationWidth;
            end
            if y_min<1
                y_min=1;
            end
            if y_max>correlationHeight
                y_max=correlationHeight;
            end
            points=SPATIALCORRELATION(y_min:y_max,x_min:x_max).*CORRELATION_WINDOW(y_min:y_max,x_min:x_max);
            
            %Options for the lsqnonlin solver
            options=optimset('MaxIter',1200,'MaxFunEvals',5000,'TolX',5e-6,'TolFun',5e-6,...
                'LargeScale','off','Display','off','DiffMinChange',1e-7,'DiffMaxChange',1,...
                'Algorithm','levenberg-marquardt');
            
            %Initial values for the solver
%             x0=[M(i) 1 shift_locx shift_locy];
              x0 = [corr_max_val(i), 1, 1, shift_locx, shift_locy, 0];

            [xloc, yloc] = meshgrid(x_min:x_max,y_min:y_max);

            %Run solver; default to 3-point gauss if it fails
            xvars=lsqnonlin(@leastsquares2D,x0,[],[],options,points(:),[yloc(:),xloc(:)],method);
            shift_errx=xvars(4)-shift_locx;
            shift_erry=xvars(5)-shift_locy;
            corr_peak_diameter(i) = sqrt(sigma^2/(2*abs(xvars(2))));

        end
        
        if method==1
%             error

            %%%%%%%%%%%%%%%%%%%%
            % 3-Point Gaussian %
            %%%%%%%%%%%%%%%%%%%%
            
            if isempty(shift_errx)
                %gaussian fit
                % Log of correlation @ minus 1 
                lCm1 = log(SPATIALCORRELATION( shift_locy , shift_locx-1 )* CORRELATION_WINDOW( shift_locy , shift_locx-1 ));
                lC00 = log(SPATIALCORRELATION( shift_locy , shift_locx   )* CORRELATION_WINDOW( shift_locy , shift_locx   ));
                % 
                lCp1 = log(SPATIALCORRELATION( shift_locy , shift_locx+1 )* CORRELATION_WINDOW( shift_locy , shift_locx+1 ));
                
                if (2*(lCm1+lCp1-2*lC00)) == 0
                    shift_errx = 0;
                    Dx = nan;
                else
                    shift_errx = (lCm1-lCp1)/(2*(lCm1+lCp1-2*lC00));
                    betax = abs(lCm1-lC00)/((-1-shift_errx)^2-(shift_errx)^2);
                    Dx = sigma./sqrt((2*betax));
                end
            else
                Dx = nan;
            end
            
            if isempty(shift_erry)
                lCm1 = log(SPATIALCORRELATION( shift_locy-1 , shift_locx )* CORRELATION_WINDOW( shift_locy-1 , shift_locx ));
                lC00 = log(SPATIALCORRELATION( shift_locy   , shift_locx )* CORRELATION_WINDOW( shift_locy   , shift_locx ));
                lCp1 = log(SPATIALCORRELATION( shift_locy+1 , shift_locx )* CORRELATION_WINDOW( shift_locy+1 , shift_locx ));
                if (2*(lCm1+lCp1-2*lC00)) == 0
                    shift_erry = 0;
                    Dy = nan;
                else
                    shift_erry = (lCm1-lCp1)/(2*(lCm1+lCp1-2*lC00));
                    betay = abs(lCm1-lC00)/((-1-shift_erry)^2-(shift_erry)^2);
                    Dy = sigma./sqrt((2*betay));
                end
            else
                Dy = nan;
            end
            
            corr_peak_diameter(i) = nanmean([Dx Dy]);
        end
        
        peak_shift_cols(i)=cc_x(shift_locx)+shift_errx;
        peak_shift_rows(i)=cc_y(shift_locy)+shift_erry;
        
        if isinf(peak_shift_cols(i)) || isinf(peak_shift_rows(i))
            peak_shift_cols(i)=0; peak_shift_rows(i)=0;
        end
    end
end

end