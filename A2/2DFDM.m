function [ M ] = 2DFDM( xbc1,xbc2,ybc1,ybc2,L,W,size )
M(:,0)=xbc1;
M(:,L)=xbc2;
M(1,:)=ybc1;
M(W,:)=ybc2;
for y=2:W-size
    if y==1
       
    for x=2:L-size
        
    end
end


end

