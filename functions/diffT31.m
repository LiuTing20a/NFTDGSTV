function diffT_a = diffT3(a,b,sizeD)

tenX = reshape(a, sizeD);
tenY = reshape(b, sizeD);

dfx     = diff(tenX, 1, 1);
dfy     = diff(tenY, 1, 2);


dfxT   = zeros(sizeD);
dfyT   = zeros(sizeD);

dfxT(1,:,:) = tenX(end, :, :) - tenX(1, :, :); %
dfxT(2:end,:,:) = -dfx;
dfyT(:,1,:)     =  tenY(:,end,:) - tenY(:,1,:);
dfyT(:,2:end,:) = -dfy;

diffT_a = dfxT + dfyT;
end