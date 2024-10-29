function M=costrrmat(X)


M=zeros(4,4);
Rx=[1 0 0;...
    0 cos(X(1)) -sin(X(1));...
    0 sin(X(1)) cos(X(1))];

Ry=[cos(X(2)) 0 sin(X(2));...
    0 1 0;...
    -sin(X(2)) 0 cos(X(2))];

Rz=[cos(X(3)) -sin(X(3)) 0;...
    sin(X(3)) cos(X(3)) 0;...
    0 0 1];

M(1:3,1:3)=Rz*Ry*Rx;
%M(1:3,1:3)=M(1:3,1:3)';
M(1,4)=X(4);
M(2,4)=X(5);
M(3,4)=X(6);
M(4,4)=1;
