function [R] =  getRotMatrixfrom2Vec(v1,v2)

v = cross(v1,v2);
c = dot(v1,v2);

v_sk = [0 -v(3) v(2);
        v(3) 0 -v(1);
        -v(2) v(1) 0];

R_ = eye(3,3) + v_sk + v_sk^2*(1/(1+c));

if not(isnan(R_(1,1)))
    r = rodrigues(R_);

   
    R = rodrigues(r);
else
    R = R_;
end
end