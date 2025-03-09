% The Thomas Algorithm is specifically used for determining the solutions
% of a tridiagonal matrix by use of gauss elimination
m=60;
n=60;

left_matrix=zeros(m,n);
right_matrix=zeros(m,1);

%Generating the coefficients for the tridiagonal matrix
for i=1:m
    r0=rand(1,1);
    r1=rand(1,1);
    r2=rand(1,1);

    left_matrix(i,i)=r1;
    if(i-1>=1) 
        left_matrix(i,i-1)=r0;
    end

    if(i+1<=m)
        left_matrix(i,i+1)=r2;
    end
end

%Generating right hand side constant matrix
for i=1:m
    r0=rand(1,1);
    right_matrix(i,1)=r0;
end

% Acquiring the solution on basis of the inverse method
soln=left_matrix\right_matrix


%Applying thomas algorithm

newleft_matrix=zeros(m,m);
newright_matrix=zeros(m,1);

%generating the modified left coefficient matrix and the right constant
%matrix simultaneously
for i=1:m
    newleft_matrix(i,i)=1;
    ri=right_matrix(i,1);
    if(i==1)
        b=left_matrix(i,i);
        c=left_matrix(i,i+1);
        gamma_n=c/b;
        newleft_matrix(i,i+1)=gamma_n;
        newleft_matrix(i,i)=1;

        rho_n=ri/b;
    elseif(i==m)
         a=left_matrix(i,i-1);
        b=left_matrix(i,i);
        newleft_matrix(i,i)=1;

        rho_n=(ri-a*newright_matrix(i-1,1))/(b-a*newleft_matrix(i-1,i));
    else
        a=left_matrix(i,i-1);
        b=left_matrix(i,i);
        c=left_matrix(i,i+1);
        gamma_n=c/(b-a*(newleft_matrix(i-1,i)));
        newleft_matrix(i,i)=1;
        newleft_matrix(i,i+1)=gamma_n;

        rho_n=(ri-a*newright_matrix(i-1,1))/(b-a*newleft_matrix(i-1,i));
    end
    newright_matrix(i,1)=rho_n;
end


for i=m:-1:1
    if(i==m) soln(i,1)=newright_matrix(i,1);
    else
        soln(i,1)=newright_matrix(i,1)-newleft_matrix(i,i+1)*soln(i+1,1);
    end
end
soln %Numeric solution
