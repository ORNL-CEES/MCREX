function [X, VAR, NWALKS]=MC_forward_adapt(A, b, P, cdf, stat)

if stat.adapt_cutoff==1 && stat.adapt_walks==1
    
    n_walks=stat.nwalks;
    max_step=stat.max_step;
    walkcut=stat.walkcut;
    var_cut=stat.varcut;

    X=zeros(size(b));
    NWALKS=zeros(size(b));
    VAR=zeros(size(b));

   for k=1:size(b,1)
       display(strcat('k= ',num2str(k)))
       count=0;
       x=[];
       ratio=10^6;
       while count < (10^(7))*size(b,1) && ratio > var_cut
            for walk=1:n_walks
                estim=0;
                previous=k;
                W=1;
                Wf=walkcut*W;
                current=k;
                estim=estim+W*b(current);
                i=1;
                while i<=max_step && W>Wf
                aux=rand;
                cdfrow_ind=find(cdf(previous,:));
                if ~isempty(cdfrow_ind)
                    current=find(cdf(previous,cdfrow_ind)>aux, 1 );
                    current=cdfrow_ind(current);
                    W=W*A(previous,current)/P(previous,current);
                else 
                    W=0;
                end

                    if W==0
                        break;
                    end
                    estim=estim+W*b(current);
                    i=i+1;
                    previous=current;
                end
                x=[x; estim];
            end
            dvar=sqrt(var(x)/size(x,1));
            if abs(mean(x))~=0
               ratio=dvar/abs(mean(x));
            else
               ratio=0;
            end       
            count=count+n_walks;
            display(strcat('ratio =', num2str(ratio)));
       end
       X(k)=mean(x);
       NWALKS(k)=count;
       VAR(k)=dvar;
   end
   
elseif stat.adapt_cutoff==0 && stat.adapt_walks==1
    n_walks=stat.nwalks;
    max_step=stat.max_step;
    var_cut=stat.varcut;

    X=zeros(size(b));
    NWALKS=zeros(size(b));
    VAR=zeros(size(b));

   for k=1:size(b,1)
       count=0;
       x=[];
       ratio=1;
       while count < 1000*size(b,1) && ratio > var_cut
            for walk=1:n_walks
                estim=0;
                previous=k;
                W=1;
                current=k;
                estim=estim+W*b(current);
                i=1;
                while i<=max_step 
                aux=rand;
                cdfrow_ind=find(cdf(previous,:));
                if ~isempty(cdfrow_ind)
                    current=find(cdf(previous,cdfrow_ind)>aux, 1 );
                    current=cdfrow_ind(current);
                    W=W*A(previous,current)/P(previous,current);
                else 
                    W=0;
                end
                    if W==0
                        break;
                    end
                    estim=estim+W*b(current);
                    i=i+1;
                    previous=current;
                end
                x=[x; estim];
            end
            dvar=sqrt(var(x)/size(x,1));
            if abs(mean(x))~=0
               ratio=dvar/abs(mean(x));
            else
               ratio=0;
            end       
           
            count=count+n_walks;
       end
       X(k)=mean(x);
       NWALKS(k)=count;
       VAR(k)=dvar;
   end
     
elseif stat.adapt_cutoff==1 && stat.adapt_walks==0   
    
    n_walks=stat.nwalks;
    max_step=stat.max_step;
    walkcut=stat.walkcut;
 
    X=zeros(size(b));
    NWALKS=zeros(size(b));
    VAR=zeros(size(b));

   for k=1:size(b,1)
       count=0;
       x=zeros(n_walks,1);

       for walk=1:n_walks
           estim=0;
           previous=k;
           W=1;
           Wf=walkcut*W;
           current=k;
           estim=estim+W*b(current);
           i=1;
           while i<=max_step && W>Wf
                aux=rand;
                cdfrow_ind=find(cdf(previous,:));
                if ~isempty(cdfrow_ind)
                    current=find(cdf(previous,cdfrow_ind)>aux, 1 );
                    current=cdfrow_ind(current);
                    W=W*A(previous,current)/P(previous,current);
                else 
                    W=0;
                end
               if W==0
                   break;
               end
               estim=estim+W*b(current);
               i=i+1;
               previous=current;
           end
           x(walk)=estim;          
       end
       count=count+n_walks;
       X(k)=mean(x);
       NWALKS(k)=count;
       dvar=sqrt(var(x)/size(x,1));
       VAR(k)=dvar;
   end

else
    n_walks=stat.nwalks;
    max_step=stat.max_step;
 
    X=zeros(size(b));
    NWALKS=zeros(size(b));
    VAR=zeros(size(b));

   for k=1:size(b,1)
       count=0;
       x=zseros(n_walks);

       for walk=1:n_walks
           estim=0;
           previous=k;
           W=1;
           current=k;
           estim=estim+W*b(current);
           i=1;
           while i<=max_step 
                aux=rand;
                cdfrow_ind=find(cdf(previous,:));
                if ~isempty(cdfrow_ind)
                    current=find(cdf(previous,cdfrow_ind)>aux, 1 );
                    current=cdfrow_ind(current);
                    W=W*A(previous,current)/P(previous,current);
                else 
                    W=0;
                end
               if W==0
                   break;
               end
               estim=estim+W*b(current);
               i=i+1;
               previous=current;
           end
           x(walk)=estim;
            
       end
       count=count+n_walks;       
       X(k)=mean(x);
       dvar=sqrt(var(x)/size(x,1));
       NWALKS(k)=count;
       VAR(k)=dvar;
   end

   
end