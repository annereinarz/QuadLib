function  test_quad2PhyP(  )

[X,Y]=meshgrid(0:0.5:1, 0:0.5:1);
X=reshape(X, 9,1);
Y=reshape(Y, 9,1);
X=[X; zeros(18,1)];
Y=[Y; (0:1/17:1)'];
x=[Y,X];
plot(x(:,1),x(:,2),'rs','MarkerSize',10)

vertexlist=[1 0 -1;
            1 2  0];
 v0 = vertexlist(:, 1);
 A(:,1:2) = bsxfun(@minus,vertexlist(:, 2:3 ),v0); 
x1=bsxfun(@plus, v0, A*x')';



vertexlist=[1 0 4;
            1 2 0.5];
v0 = vertexlist(:, 1);
A(:,1:2) = bsxfun(@minus,vertexlist(:, 2:3 ),v0);
x2=bsxfun(@plus, v0, A*x')';


plot(x1(:,1),x1(:,2),'rs','MarkerSize',10);
hold on;
plot(x2(:,1),x2(:,2),'bx','MarkerSize',10);


end

