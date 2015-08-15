function T=get_transfer_matrix1(n,w,d,c)
I=sqrt(-1);

T=[cos(n*w*d/c)+I*sin(n*w*d/c) 0;0 cos(n*w*d/c)-I*sin(n*w*d/c)];
end