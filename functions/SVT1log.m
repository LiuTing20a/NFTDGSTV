
function [ myM ] = SVT1log(myM, mu,speedup)
if(speedup==0)
%         [U,S,V]=svd(myM,'econ');
        [m, n] = size(myM);
if m <= n
    AAT = myM*myM';
    [S, Sigma, ~] = svd(AAT);
    Sigma         = sqrt(diag(Sigma));    
    tol           = eps;
    temp          = (Sigma-tol).^2-4*(mu-tol*Sigma); % log(eps) = 36.0437
    ind           = find (temp>0);
    n             = length(ind);
    SigmaNew      = (Sigma(1:n)-tol+sqrt(temp(1:n)))/2;
    SigmaNew      = SigmaNew ./ Sigma(1:n);
  myM = S(:, 1:n) * diag(SigmaNew) * S(:, 1:n)' * myM;
    return;
end
if m > n
    myM = SVTlog(myM', mu);
   myM = myM';
    return;
end    
    else
%         [U,S,V]=FastSVD(myM,100);
        [m, n] = size(myM);
if m <= n
    AAT = myM*myM';
    [S, Sigma, ~] = FastSVD(AAT,100);
    Sigma         = sqrt(diag(Sigma));    
    tol           = eps;
    temp          = (Sigma-tol).^2-4*(mu-tol*Sigma); % log(eps) = 36.0437
    ind           = find (temp>0);
    n             = length(ind);
    SigmaNew      = (Sigma(1:n)-tol+sqrt(temp(1:n)))/2;
    SigmaNew      = SigmaNew ./ Sigma(1:n);
  myM = S(:, 1:n) * diag(SigmaNew) * S(:, 1:n)' * myM;
    return;
end
if m > n
    myM = SVTlog(myM', mu);
   myM = myM';
    return;
end
%     [row,col]=size(S);
%     bound=min([row col]);
%     for i=1:bound
%         if(S(i,i)>mu)
%             S(i,i)=S(i,i)-mu;
%         else
%             S(i,i)=0;
%         end
%     end
%     myM=U*S*V'; 
%  myM = S(:, 1:n) * diag(SigmaNew) * S(:, 1:n)' * myM;
end

