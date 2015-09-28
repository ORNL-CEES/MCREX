% I need to build an undirected graph
S = A + A';
S=sparse(S);

G=spones(S);

% number of vertices
nv=size(G,1);

% number of edges
D=diag(G);

% construction of the matrix without loops
 Gnl=sparse(G-diag(D));
 ne=nnz(Gnl)/2 ;

fileID = fopen('graph.txt','w');
fprintf(fileID,'%u %u\n', nv, ne);
for i=1:size(Gnl,1)
    
    aux=find(Gnl(i,:));
    fprintf(fileID,'%u ',aux(1:end-1));
    fprintf(fileID,'%u',aux(end));
    fprintf(fileID,'\n');
    
end
fclose(fileID);