function [B,Mtt]=KernelDeriv_single(X,Y,K,SGX,SGY,EPS,Xdata,obsY)

    [N,M]=size(X);  % N: data size, M: dim of X.
    Xdatam=repmat(Xdata,N,1);
    Ydatam=repmat(obsY,N,1);
    I=eye(N);

    sx2=2*SGX*SGX;
    sy2=2*SGY*SGY;

    % Gram matrix of X
    ab=X*X';
    aa=diag(ab);
    D=repmat(aa,1,N);
    xx=max(D + D' - 2*ab, zeros(N,N));
    Kx=exp(-xx./sx2);  

    % Gram matrix of Y
    ab=Y*Y';
    aa=diag(ab);
    D=repmat(aa,1,N);
    yy=max(D + D' - 2*ab, zeros(N,N));
    Ky=exp(-yy./sy2);  

    % Derivative of k(X_i, x) w.r.t. x
%     Dx=reshape(repmat(X,N,1),N,N,M);
%     Xij=Dx-permute(Dx,[2 1 3]);
%     Xij=Xij./SGX/SGX;
%     H=Xij.*repmat(Kx,[1 1 M]);
%     Hp=permute(H,[1,3,2]);
    Xij=Xdatam-X;
    diff=Xij.^2;
    Hp=Xij./(SGX*SGX).*exp(-diff./sx2);
    
    F=((Kx+N*EPS.*I)\Ky)/(Kx+N*EPS.*I);
    Mt=zeros(M,M,N);
    Mtt=zeros(M,M);
   % parfor i=1:N
   %     Mt(:,:,i)=Hp(:,:,i)'*F*Hp(:,:,i)*W(i);
%         fprintf('%d/n',i);
  %  end
  Mtt(:,:)=Hp'*F*Hp;
  %  R=sum(Mt,3);
   % [B,L]=eigs(R,K);
   [B,L]=eigs(Mtt,K);

end
