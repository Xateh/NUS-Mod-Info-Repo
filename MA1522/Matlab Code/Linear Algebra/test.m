%{
syms a b;

A = [a,2,a,a+b;
    a,2,a,a;
    3,3,-b,3;
    a+1,3,a+1,a+1];
B = [a-b;
    a-b;
    -b;
    a-b+1];

display(solution_space(A, B));
%}

syms a

A = [a a*(a-1) 0; 0 a-1 (a-1)*(a+1); 0 0 a+1];

display(solution_space(A));