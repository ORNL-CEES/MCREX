function [X, VAR, NWALKS, reject]=MC_forward_adapt2(A, b, P, cdf, stat)

reject=zeros(size(b));

if stat.adapt_cutoff==1 && stat.adapt_walks==1
    
    n_walks=stat.nwalks;
    max_step=stat.max_step;
    walkcut=stat.walkcut;
    var_cut=stat.varcut;
    var_diff=stat.vardiff;

    X=zeros(size(b));
    NWALKS=zeros(size(b));
    VAR=cell(size(b));

   for k=1:size(b,1)
       count=0;
       x=[];
       ratio=1;
       
       for j=1:2
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
                dvar=sqrt(var(x)/size(x,1));
                if abs(mean(x))~=0
                    ratio=dvar/abs(mean(x));
                else
                    ratio=0;
                end
            end
            count=count+n_walks;  
            VAR{k}=[VAR{k} dvar];
       end
       
       while   ratio > var_cut || ( ratio > 0 && abs((VAR{k}(end)-VAR{k}(end-1)))/abs(mean(x)) > var_diff ) 
       %while   ( ratio > 0 && abs((VAR{k}(end)-VAR{k}(end-1)))/(VAR{k}(end)*sqrt(length(x))) > var_diff ) 
               
               %display(strcat('ratio =', num2str(ratio), '   ', num2str(count)))
               %display(strcat('ratio =', num2str(abs((VAR{k}(end)-VAR{k}(end-1)))/abs(mean(x)))))
               if ratio < var_cut 
                   reject(k)=reject(k)+1;
               end
                x2=zeros(n_walks,1);
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
               end
               x2(walk)=estim;
               x=[x; x2];
               dvar=sqrt(var(x)/size(x,1));
               if abs(mean(x))~=0
                  ratio=dvar/abs(mean(x));
               else
                   ratio=0;
               end            
              count=count+n_walks;
              VAR{k}=[VAR{k} dvar];
       end
       
       X(k)=mean(x);
       NWALKS(k)=count;
   end

   
elseif stat.adapt_cutoff==0 && stat.adapt_walks==1
    n_walks=stat.nwalks;
    max_step=stat.max_step;
    var_cut=stat.varcut;
    var_diff=stat.vardiff;

    X=zeros(size(b));
    NWALKS=zeros(size(b));
    VAR=zeros(size(b));

   for k=1:size(b,1)
       count=0;
       x=[];
       ratio=1;
       
       for j=1:2
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
                    current=min(find(cdf(previous,cdfrow_ind)>aux));
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
                dvar=sqrt(var(x)/size(x,1));
                if abs(mean(x))~=0
                    ratio=dvar/abs(mean(x));
                else
                    ratio=0;
                end

            end
            count=count+n_walks;
            VAR{k}=[VAR{k} dvar];
       end
           
       
       while   ratio > var_cut || ( ratio > 0 && abs((VAR{k}(end)-VAR{k}(end-1)))/abs(mean(x)) > var_diff ) 
                               
                if ratio < var_cut 
                   reject(k)=reject(k)+1;
                end          

                x2=zeros(n_walks,1);
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
                        current=min(find(cdf(previous,cdfrow_ind)>aux));
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
                    x2(walk)=estim;
                    x=[x; x2];
                    dvar=sqrt(var(x)/size(x,1));
                    if abs(mean(x))~=0
                        ratio=dvar/abs(mean(x));
                    else
                        ratio=0;
                    end

                end
                count=count+n_walks;
                VAR{k}=[VAR{k} dvar];
           
                X(k)=mean(x);
                NWALKS(k)=count;
       end
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
                    current=min(find(cdf(previous,cdfrow_ind)>aux));
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
           dvar=sqrt(var(x)/size(x,1));
       end
       count=count+n_walks;

       X(k)=mean(x);
       NWALKS(k)=count;
       VAR{k}=dvar;
   end

else
    n_walks=stat.nwalks;
    max_step=stat.max_step;
 
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
           current=k;
           estim=estim+W*b(current);
           i=1;
           while i<=max_step 
                aux=rand;
                cdfrow_ind=find(cdf(previous,:));
                if ~isempty(cdfrow_ind)
                    current=min(find(cdf(previous,cdfrow_ind)>aux));
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
           dvar=sqrt(var(x)/size(x,1));
       end
       count=count+n_walks;

       X(k)=mean(x);
       NWALKS(k)=count;
       VAR{k}=dvar;
   end

   
end