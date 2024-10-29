function plot3c(x,y,z,varargin)  
    x = [x(:) ; x(1)];   
    y = [y(:) ; y(1)];
    z = [z(:) ; z(1)];
    plot3(x,y,z,varargin{:})  
end