function [tensor_X,  tensor_S ,change] = NFTDGSTV(tenD, lambda)

% Solve the WSNM problem 
% ---------------------------------------------
% Input:
%       tenD       -    n1*n2*n3 tensor
%       lambda  -    >0, parameter
%       mu-      regelarization parameter
%       p    -    key parameter of Schattern p norm
%
% Output:
%       tenB       -    n1*n2*n3 tensor
%       tenT       -    n1*n2*n3 tensor
%       change     -    change of objective function value

max_iter = 100;
lambda1 =0.001;%0.001
[row, col, channel] = size(tenD);
tsize  = [row, col, channel];
sizeD = size(tenD);
p = sizeD(3);
change=zeros(1,max_iter);
lambda_2   =0.01;
lambda_3=lambda;
N          = 3;
betat =0.002;
gamma =0.001;
deta = 0.002;
alpha      = [1/N, 1/N, 1/N];
myInitial_v=0.0002;% tuning
factor = 1.15;
M = sizeD(1);
N2 = sizeD(2);
Eny_x_fft   = (abs(psf2otf([+1; -1], [M,N2,p]))).^2  ;
Eny_y_fft   = (abs(psf2otf([+1, -1], [M,N2,p]))).^2  ;
Eny_fft  =  Eny_x_fft + Eny_y_fft;
R1 = zeros(sizeD);      % R1 : auxiliary variable for GS_x
R2 = zeros(sizeD);      % R2 : auxiliary variable for GS_y
W3 = R1;                % W3 : multiplier for DQ_x-R1
W4 = R2;  
W1 = zeros(size(tenD));  % W1 : multiplier for Y - X - S          
W2 = W1; 

    U=cell(N,1);
    V=cell(N,1);
    myV=cell(N,1);
    
    X=cell(N,1);
    W=cell(N,1);
    Gamma=cell(N,1);
    epsilon=1.0e-6;%1.0e-6
    for i=1:N
        q=1;
        for j=[1:i-1,i+1:N]
            q=q*tsize(j);
        end
        rand('seed',1)
        U{i}=rand(tsize(i),tsize(i));   
        rand('seed',2)
        V{i}=rand(tsize(i),tsize(i));   
    end
    
    for i=1:N
       Gamma{i}=sign(U{i})/max([norm(U{i}), norm(U{i},Inf), epsilon]);      
    end 
    tensor_X=tenzeros(tsize);
    tensor_G=tensor_X;
    tensor_S=tensor_X;          
    tensor_W=tenzeros(tsize);
    tensor_K=tenzeros(tsize);  
    Q = tensor_X;  
  %% Parameter settings
    iteration=1;
    rho_1=myInitial_v;
    rho_2=myInitial_v;
    rho_3=myInitial_v;%%
    rho_4=myInitial_v;
    rho_5= 0.005; %tuning 0.01
  
    while(true)
        tensor_X_pre=tensor_X;
        for n=1:N
            X{n}=double(tenmat(tensor_X,n));
            W{n}=double(tenmat(tensor_W,n));   %% XµÄ³Ë×Ó    
        end
    %% update U 
        for n=1:N
              U{n}=SVT1log(V{n} + Gamma{n}/rho_3,alpha(n)/rho_3,0);  
        end
    %% update V 
       for n=1:N
            tmp=ttm(tensor_G,V,-n);
            tmp=tenmat(tmp,n);
            tmp=double(tmp);
            V{n}=(- Gamma{n} + rho_3 * U{n} + W{n} * tmp' + rho_4 * X{n}* tmp')/((rho_3*eye(tsize(n)) + rho_4*(tmp*tmp'))); 
        end
    %% update X 
    tmp = - tensor_W + rho_4* ttm(tensor_G,V,1:N) + rho_5*tenD - rho_5*tensor_S + tensor_K+ betat*Q-W2;
    tensor_X = tmp/(rho_4 + rho_5+betat);
    %% update G 
    for n=1:N
        myV{n}=V{n}';
    end
    myG = optimize_Z(myV,double(tensor_X),double(tensor_W),rho_4,lambda_2);
    tensor_G = tensor(myG);        
    %% Update sparse target tensor S
    tensor_S = prox_l1(double(-tensor_X+tenD+tensor_K/rho_5),lambda_3/rho_5);%
    %% - Q subproblem update     
    diffT_p = diffT31(betat * R1 - W3, betat * R2 - W3, sizeD);
    temp1 = reshape(diffT_p + betat*tensor_X + W2, sizeD);
    z = real(ifftn(fftn(double(temp1)) ./ (betat*Eny_fft + betat)));
    Q = reshape(z,sizeD);
    %% - R1 and R2 subproblem update
    [diff_Qx, diff_Qy] = diff31(Q, sizeD); 
    R1 = Thres_21(diff_Qx+ W3/betat, lambda1/betat);  
    R2 = Thres_21(diff_Qy+ W4/betat, lambda1/betat); 
    %% update multiplers
     for n=1:N
        Gamma{n}=Gamma{n}+ rho_3*(V{n}-U{n});
     end
        tensor_W = tensor_W + rho_4*(tensor_X - ttm(tensor_G, V,1:N));
        tensor_K = tensor_K + rho_5*(tenD - tensor_X - tensor_S);
        W2 = W2+betat*(tensor_X-Q);   
        W3 = W3+betat*(diff_Qx-R1);    
        W4 = W4+betat*(diff_Qy-R2);
        diff=norm(tensor_X-tensor_X_pre)/norm(tensor_X);
        change(iteration)=(diff);
        gamma = 1.1*gamma;  
        deta  = 1.1*deta;  
        betat = 1.25*betat;
        rho_1=rho_1*factor;
        rho_2=rho_2*factor;
        rho_3=rho_3*factor;
        rho_4=rho_4*factor;
        rho_5=rho_5*factor;
        fprintf('iter=%d,diff=%f\n',iteration,diff);
        if(iteration>80)||diff<1e-5
            break;
        end
    end          
end