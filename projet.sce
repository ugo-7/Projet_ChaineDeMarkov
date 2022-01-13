N=50
p=19/20

// Espaces d'états E et O

E=[0:N;N:-1:0];

function res=nbcase(n)
    if(n==1)
        res=3;
    else
        res=nbcase(n-1)+n+1;
    end
endfunction

O=zeros(nbcase(N), 2);
Obis=zeros((N+1)^2,2);
for(i=0:N)
    Obis(((N+1)*i+1):((N+1)*i+N+1),1)=linspace(i,i,N+1);
    Obis(((N+1)*i+1):((N+1)*i+N+1),2)=[0:N];
end
k=1
for(i=1:size(Obis)(1))
    if(Obis(i,1)+Obis(i,2)<=N)
        O(k,1)=Obis(i,1);
        O(k,2)=Obis(i,2);
        k=k+1;
    end
end


// Matrice P

P=(1/N)*(diag(1:N,-1)+diag(N:-1:1,1));

// Matrice e

function res=bin(k, n)
    if(k<=n)
        res=prod(1:n)/(prod(1:k)*prod(1:(n-k)));
    else
        res=0;
    end
endfunction

e=zeros(size(O)(1), size(E)(2));
for(i=1:size(O)(1))
    for(j=1:size(E)(2))
        xA=O(i,1);
        xB=O(i,2);
        yA=E(1,j);
        yB=E(2,j);
        e(i,j)=bin(xA,yA)*p^xA*(1-p)^(yA-xA)*bin(xB,yB)*p^xB*(1-p)^(yB-xB);
    end
end

// Matrice Pi 

pi = linspace(1/(N+1),1/(N+1),N+1);

// Fonction alpha

function ret=alpha(x, L)
    ret=0
    for y=1:N+1
        if L==1 then
            res=pi(y)*e(x(1),y)
        else 
            res=0
            for i=1:N+1
                res=res+(delta(x, i, L-1)*P(i,y)) 
            end
            res=res*e(x(L),y)
        end
        ret=ret+res
    end
endfunction

// Fonction delta

function ret=delta(x, y, L)
    if L==1 then
        res=pi(y)*e(x(1),y)
    else 
        for z=1:N+1
            u(z)=delta(x, z, L-1)*P(z,y)
        end
        res=max(u)*e(x(L),y)
    end
    ret=res
endfunction

// Fonction Psi

function res=psi(x, y, L)
    for z=1:N+1
        u(z)=delta(x, z, L-1)*P(z,y)
    end
    [a, res]=max(u)
endfunction

// Fonctions Yet

function res=yLet(x, L)
    if L==length(x) then
        for y=1:N+1
            u(y)=delta(x, y, L)
        end
        [a, b]=max(u)
    else
        b = psi(x,yLet(x, L+1), L+1)
    end
    res=b
endfunction

function res=yet(x)
    res=1:length(x)
    for k=1:length(x)
        res(k) = yLet(x, k)
    end
endfunction

// Application

size(O)
a=998
O(a,1)
O(a,2)
x1=[1000,1020,998]
yet(x1)
alpha(x1,3)
