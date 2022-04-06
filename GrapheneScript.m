myInd = [];
for j = 1:1296
    for i = 1:50
        if (.01 <= MyEig(i,j)) && (MyEig(i,j) <= .04)
           myInd = [myInd; i j];
        end
    end
end
[numpoints,~] = size(myInd);
for i = 1:numpoints
     scatter3(kPoints(1,myInd(i,2)), kPoints(2,myInd(i,2)), MyEig(myInd(i,1),myInd(i,2))); hold on
end
