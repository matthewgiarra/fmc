function [UVAL, VVAL, ISOUTLIER] = validateField_prana(X, Y, U, V)
% Size of the grid
[gridHeight, gridWidth] = size(X);

% Hard code a bunch of validation parameters
Eval = zeros(size(U));
ThreshSwitch = 0;
UODswitch = 1;
BootSwitch = 0;
extraPeakSwitch = 0;
Uthresh = [-inf inf];
Vthresh = [-inf inf];
UODthresh = [3 2];
UODwinsize = [3 3; 3 3];
BootPer = 0;
BootIter = 0;
BootKMax = 0;
c = zeros(length(X(:)), 1);
d = zeros(length(X(:)), 1);

v = -V;
y = flipud(Y);

[UvalRaw, VvalRaw, IsOutlierRaw]=VAL(X(:), y(:), U(:), v(:), Eval(:), c(:), d(:),ThreshSwitch,UODswitch,BootSwitch,extraPeakSwitch,...
                        Uthresh,Vthresh,UODwinsize,UODthresh,BootPer,BootIter,BootKMax);


% % Vector validation
% [UvalRaw, VvalRaw, RvalRaw, ISOUTLIER] = VAL_2d3c(X(:), Y(:), U(:), V(:), W(:), Eval(:), C(:), D(:), ThreshSwitch,UODswitch,...
%     BootSwitch, extraPeakSwitch, Uthresh, Vthresh, UODwinsize, UODthresh, BootPer,BootIter,BootKMax);

% Reshape matrices.
UVAL = flipud(reshape(UvalRaw, [gridHeight, gridWidth]));
VVAL = -1 * flipud(reshape(VvalRaw, [gridHeight, gridWidth]));
ISOUTLIER = flipud(reshape(IsOutlierRaw, [gridHeight, gridWidth]));
   
end

function [Uval,Vval,Evalval,Cval,Dval,DXval,DYval,ALPHAval]=VAL(X, Y, U, V, Eval,C,D,Threshswitch,UODswitch,Bootswitch,extrapeaks,Uthresh,Vthresh,UODwinsize,UODthresh,Bootper,Bootiter,Bootkmax,DX,DY,ALPHA)
% --- Validation Subfunction ---
if extrapeaks
    j=3;
else
    j=1;
end

imClass = 'double';

if exist('DX','var')
    [~,~,~,~,DX,DY,ALPHA]=matrixform(X,Y,U,V,DX,DY,ALPHA);
    [X,Y,U,V,Eval,C,D]=matrixform(X,Y,U,V,Eval,C,D);
else
    [X,Y,U,V,Eval,C,D]=matrixform(X,Y,U,V,Eval,C,D);
    DX    = zeros(size(D),imClass);
    DY    = zeros(size(D),imClass);
    ALPHA = zeros(size(D),imClass);
end



Uval=U(:,:,1);Vval=V(:,:,1);Evalval=Eval(:,:,1);
if ~isempty(C)
    Cval=C(:,:,1);
    Dval=D(:,:,1);
else
    Cval=[];
    Dval= [];
end

    DXval    = DX(:,:,1);
    DYval    = DY(:,:,1);
    ALPHAval = ALPHA(:,:,1);
    
    
    
S=size(X);

if Threshswitch || UODswitch
    for i=1:j
        %Thresholding
        if Threshswitch
            [Uval,Vval,Evalval] = Thresh(Uval,Vval,Uthresh,Vthresh,Evalval);
        end

        %Univeral Outlier Detection
        if UODswitch
            t=permute(UODwinsize,[2 3 1]);
            t=t(:,t(1,:)~=0);
            [Uval,Vval,Evalval] = UOD(Uval,Vval,t',UODthresh,Evalval);
        end
%         disp([num2str(sum(sum(Evalval>0))),' bad vectors'])
        %Try additional peaks where validation failed
        if i<j
            Utemp=U(:,:,i+1);Vtemp=V(:,:,i+1);Evaltemp=Eval(:,:,i+1);Ctemp=C(:,:,i+1);Dtemp=D(:,:,i+1);
            DXtemp=DX(:,:,i+1);DYtemp=DY(:,:,i+1);ALPHAtemp=ALPHA(:,:,i+1);
            Uval(Evalval>0)=Utemp(Evalval>0);
            Vval(Evalval>0)=Vtemp(Evalval>0);
            Cval(Evalval>0)=Ctemp(Evalval>0);
            Dval(Evalval>0)=Dtemp(Evalval>0); 
            Evalval(Evalval>0)=Evaltemp(Evalval>0);
            
            DXval(Evalval>0)=DXtemp(Evalval>0); 
            DYval(Evalval>0)=DYtemp(Evalval>0); 
            ALPHAval(Evalval>0)=ALPHAtemp(Evalval>0); 
        end
    end
end

%limit replacement search radius to largest region tested during UOD
maxSearch = floor( (max(UODwinsize(:))-1)/2 );

%replacement
for i=1:S(1)
    for j=1:S(2)
        if Evalval(i,j)>0
            %initialize replacement search size
            q=0;
            s=0;

            %get replacement block with at least 8 valid points, but stop if region grows bigger than UOD test region
            while s==0
                q=q+1;
                Imin = max([i-q 1   ]);
                Imax = min([i+q S(1)]);
                Jmin = max([j-q 1   ]);
                Jmax = min([j+q S(2)]);
                Iind = Imin:Imax;
                Jind = Jmin:Jmax;
                Ublock = Uval(Iind,Jind);
                if q >= maxSearch || length(Ublock(~isnan(Ublock)))>=8 
                    Xblock = X(Iind,Jind)-X(i,j);
                    Yblock = Y(Iind,Jind)-Y(i,j);
                    Vblock = Vval(Iind,Jind);
                    s=1;
                end
            end
            
            %distance from erroneous vector
            Dblock = (Xblock.^2+Yblock.^2).^-0.5;
            Dblock(isnan(Ublock))=nan;
            Dblock(isinf(Dblock))=nan;

            %validated vector
            Uval(i,j) = nansum(nansum(Dblock.*Ublock))/nansum(nansum(Dblock));
            Vval(i,j) = nansum(nansum(Dblock.*Vblock))/nansum(nansum(Dblock));       
        end
    end
end

%clean up any remaining NaNs
Uval(isnan(Uval)) = 0;
Vval(isnan(Vval)) = 0;

%Bootstrapping
if Bootswitch
    [Uval,Vval,Evalval] = bootstrapping(X,Y,Uval,Vval,Bootper,Bootiter,Bootkmax,Evalval);
end

%convert back to vector
[Uval,Vval,Evalval,Cval,Dval]=vectorform(X,Y,Uval,Vval,Evalval,Cval,Dval);
[DXval,DYval,ALPHAval]=vectorform(X,Y,DXval,DYval,ALPHAval,[],[]);

end



function [X,Y,U,V,Eval,C,D]=matrixform(x,y,u,v,eval,c,d)
% --- Vector to Matrix Subfunction ---

imClass = 'double';

%find unique x and y grid points
a=sort(unique(x));
b=sort(unique(y));
N=length(x);

%initialize matrices
U=nan(length(b),length(a),size(u,2));
V=nan(length(b),length(a),size(v,2));
Eval=-1*ones(length(b),length(a),size(eval,2),imClass);

%generate grid matrix
[X,Y]=meshgrid(a,b);

%generate variable matrices (nans where no data available)
for i=1:size(U,3)
    for n=1:N
        I=find(b==y(n));
        J=find(a==x(n));
        U(I,J,i) = u(n,i);
        V(I,J,i) = v(n,i);
        Eval(I,J,i) = eval(n);
    end
end
if ~isempty(c)
    C=nan(length(b),length(a),size(c,2));
    for i=1:size(c,2)
        for n=1:N
            I= b==y(n);
            J= a==x(n);
            C(I,J,i)=c(n,i);
        end
    end
else
    C=[];
end
if ~isempty(d)
    D=nan(length(b),length(a),size(d,2));
    for i=1:size(d,2)
        for n=1:N
            I= b==y(n);
            J= a==x(n);
            D(I,J,i)=d(n,i);
        end
    end
else
    D=[];
end

end


function [u,v,eval,c,d]=vectorform(x,y,U,V,Eval,C,D)
% --- Matrix to Vector Subfunction ---
x=x(:);y=y(:);
%find unique x and y grid points
a=sort(unique(x));
b=sort(unique(y));
N=length(x(:));

imClass = 'double';

%initialize vectors
S=size(x(:));
u    = zeros(S,imClass);
v    = zeros(S,imClass);
eval = zeros(S,imClass);
if ~isempty(C)
    c = zeros(S,imClass);
    d = zeros(S,imClass);
else
    c = [];
    d = [];
end

%generate data vectors where data is available
for n=1:N
    I=find(b==y(n));
    J=find(a==x(n));
    u(n)    = U(I,J);
    v(n)    = V(I,J);
    eval(n) = Eval(I,J);
    if ~isempty(C)
        c(n)    = C(I,J);
        d(n)    = D(I,J);
    end
end

end

function [U,V,Eval] = UOD(U,V,t,tol,Eval)
% --- Universal Outlier Detection Validation Subfunction ---

%number of validation passes
pass = length(tol);

S=size(U);

%outlier searching
for k=1:pass
    
    q0 = floor((t(k,:)-1)/2);
    
    for i=1:S(1)
        for j=1:S(2)
            if Eval(i,j)==0           
                %get evaluation block with at least 8 valid points
                s=0;
                q = q0;
                while s==0
                    Imin = max([i-q(2) 1   ]);
                    Imax = min([i+q(2) S(1)]);
                    Jmin = max([j-q(1) 1   ]);
                    Jmax = min([j+q(1) S(2)]);
                    Iind = Imin:Imax;
                    Jind = Jmin:Jmax;
                    Ublock = U(Iind,Jind);
                    if length(Ublock(~isnan(Ublock)))>=8 || any(q >= 2*q0)
%                         Xblock = X(Iind,Jind)-X(i,j);
%                         Yblock = Y(Iind,Jind)-Y(i,j);
                        Vblock = V(Iind,Jind);
                        s=1;
                    else
                        q=q+1;
                    end
                end

%                 %distance from vector location
%                 Dblock = (Xblock.^2+Yblock.^2).^0.5;
%                 Dblock(isnan(Ublock))=nan;

                %universal outlier detection
                Ipos = find(Iind==i);
                Jpos = find(Jind==j);
                [Ru]=UOD_sub(Ublock,Ipos,Jpos);
                [Rv]=UOD_sub(Vblock,Ipos,Jpos);

                if Ru > tol(k) || Rv > tol(k)
                    %UOD threshold condition
                    U(i,j)=nan;
                    V(i,j)=nan;
                    Eval(i,j)=k;
                end

            end

        end
    end
end

end


function [R]=UOD_sub(W,p,q)
% --- Universal Outlier Detection Algorithm ---

%minimum variance assumption
e=0.1; 

%remove value from query point
x=W(p,q);
W(p,q)=nan;

%remove any erroneous points
P=W(:);
Ps = sort(P);
Psfull = Ps(~isnan(Ps));
N=length(Psfull);

if N<=floor(length(W)/3)
    %return negative threshold value if no valid vectors in block
    R = inf;
else
    %return the median deviation normalized to the MAD
    if mod(N,2)==0
        
        % Median value of surrounding vectors
        M = (Psfull(N/2)+Psfull(N/2+1))/2;
        
        % Sorted difference between the non-NAN values and the median of
        % the surrounding vectors
        MADfull = sort(abs(Psfull-M));
        
        % Median of the medians of the surrounding vectors
        Q = (MADfull(N/2)+MADfull(N/2+1))/2;
        
        % Ratio of the magnitude of the (test vector - the median of
        % surrounding vectors) to the (median of the medians plus the expected variance)
        R = abs(x-M)/(Q+e);
    else
        M = Psfull((N+1)/2);
        MADfull = sort(abs(Psfull-M));
        Q = MADfull((N+1)/2);
        R = abs(x-M)/(Q+e);
    end
end

end