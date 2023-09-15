M = [1,2,3;4,5,6;7,8,9;10,11,12];

[dMdx, dMdy] = gradient(M, 1/3, 1)

U = zeros(4,5);

U(2,:)