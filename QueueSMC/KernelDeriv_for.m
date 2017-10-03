function [B,t,V]=KernelDeriv_for(X,Y,K,SGX,SGY,EPS,W)

    [N,M]=size(X);  % N: data size, M: dim of X.

    I=eye(N);

    sx2=2*SGX*SGX;
    sy2=2*SGY*SGY;

    % Gram matrix of X
    ab=X*X';
    aa=diag(ab);
    D=repmat(aa,1,N);
    xx=max(D + D' - 2*ab, zeros(N,N));
    Kx=exp(-xx./sx2);  
%     Wm=reshape(repmat(W',M*N,1),[N,M,N]);
    % Gram matrix of Y
    ab=Y*Y';
    aa=diag(ab);
    D=repmat(aa,1,N);
    yy=max(D + D' - 2*ab, zeros(N,N));
    Ky=exp(-yy./sy2);  

    % Derivative of k(X_i, x) w.r.t. x
    Dx=reshape(repmat(X,N,1),N,N,M);
    Xij=Dx-permute(Dx,[2 1 3]);
    Xij=Xij./SGX/SGX;
    H=Xij.*repmat(Kx,[1 1 M]);
    Hp=permute(H,[1,3,2]);
    
    F=((Kx+N*EPS.*I)\Ky)/(Kx+N*EPS.*I);
    Mt=zeros(M,M,N);
    for i=1:N
        Mt(:,:,i)=Hp(:,:,i)'*F*Hp(:,:,i)*W(i);
%         fprintf('%d/n',i);
    end
    R=sum(Mt,3);
    [V,L]=eig(R);
    [e,idx]=sort(diag(L),'descend');
    B=V(:,idx(1:K));
    t=sum(e(idx(1:K)));

end
