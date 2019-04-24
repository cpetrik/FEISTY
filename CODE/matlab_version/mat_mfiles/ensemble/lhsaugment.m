function xF = lhsaugment(x1,nPoi)
%function xF = lhsaugment(x1,nPoi)
%Function to augment a given latin hypercube x1 by a number of points,
%nPoi. Only the length is changed, i.e. points are added to the length.
%The original points are left unctouched and appear first in the output
%xF. Thus the size of xF is [size(x1,1)+nPoi size(x1,2)].
x2=lhsdesign(nPoi,size(x1,2));
nPoi=size(x2,1);
oPoi=size(x1,1);
tPoi=nPoi+oPoi;
fInt=1/tPoi;
for i=1:tPoi
    cBound(i,:)=[(i-1)*fInt i*fInt];
end
xF=zeros(tPoi,size(x1,2));
bX1=zeros(size(x1));
bX2=zeros(size(x2));
bF=zeros(tPoi,size(x1,2));
iF=zeros(1,size(x1,2));
iMove=0;
for i=1:oPoi
    for j=1:size(cBound,1)
        for l=1:size(x1,2)
            if (x1(i,l)>cBound(j,1))&&(x1(i,l)<=cBound(j,2))&&(bF(j,l)==0)
                iF(1,l)=iF(1,l)+1;
                xF(iF(1,l),l)=x1(i,l);
                bX1(i,l)=1;
                bF(j,l)=1;
            elseif (x1(i,l)>cBound(j,1))&&(x1(i,l)<=cBound(j,2))&&(bF(j,l)~=0)
                iMin=size(cBound,1);
                pMin=size(cBound,1);
                for m=j:-1:1
                    if (bF(m,l)==0)
                        iMin=m;
                        pMin=j-m;
                        break
                    end
                end
                for m=j:size(cBound,1)
                    if (bF(m,l)==0)&&(m-j<pMin)
                        iMin=m;
                        pMin=j+m;
                        break
                    end
                end
                iF(1,l)=iF(1,l)+1;
                xF(iF(1,l),l)=x1(i,l);
                bX1(i,l)=1;
                bF(iMin,l)=1;
            end
        end
    end
end
for i=1:nPoi
    for j=1:size(cBound,1)
        for l=1:size(x2,2)
            if (x2(i,l)>cBound(j,1))&&(x2(i,l)<=cBound(j,2))&&(bF(j,l)==0)
                iF(1,l)=iF(1,l)+1;
                xF(iF(1,l),l)=x2(i,l);
                bX2(i,l)=1;
                bF(j,l)=1;
            elseif (x2(i,l)>cBound(j,1))&&(x2(i,l)<=cBound(j,2))&&(bF(j,l)~=0)
                iMin=size(cBound,1);
                pMin=size(cBound,1);
                for m=j:-1:1
                    if (bF(m,l)==0)
                        iMin=m;
                        pMin=j-m;
                        break
                    end
                end
                for m=j:size(cBound,1)
                    if (bF(m,l)==0)&&(m-j<pMin)
                        iMin=m;
                        pMin=j+m;
                        break
                    end
                end
                iF(1,l)=iF(1,l)+1;
                xF(iF(1,l),l)=(x2(i,l)-(floor(x2(i,l)/fInt)*fInt))+((iMin-1)*fInt);
                bX2(i,l)=1;
                bF(iMin,l)=1;
                if l==1
                iMove=iMove+1;
                end
            end
        end
    end    
end