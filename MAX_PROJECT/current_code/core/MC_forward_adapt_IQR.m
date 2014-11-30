function [X, VAR, NWALKS]=MC_forward_adapt_IQR(A, b, P, cdf, stat)


if stat.adapt_cutoff==1 && stat.adapt_walks==1
    
    n_walks=stat.nwalks;
    max_step=stat.max_step;
    walkcut=stat.walkcut;
    var_diff=stat.vardiff;

    X=zeros(size(b));
    NWALKS=zeros(size(b));
   % IQR=cell(size(b,1),1);
   VAR=cell(size(b,1),1);

   for k=1:size(b,1)
       count=0;
       x=[];
       y=[]; 
       
       
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
                if sum(abs(P(previous,:)))>0
                    current=min(find(cdf(previous,:)>aux));
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
%                 if (std(x)~=0)
%                     %Iqr=iqr(x)/std(x);
%                     VAR{k}=[VAR{k} var(x)];
%                 else
%                     %Iqr=0;
%                     VAR{k}=[VAR{k} 1];
%                 end
            end
            count=count+n_walks;  
            %IQR{k}=[IQR{k} Iqr];
       if (std(x)~=0)
            %Iqr=iqr(x)/std(x);
            VAR{k}=[VAR{k} var(x)];
        else
            %Iqr=0;
            VAR{k}=[VAR{k} 1];
       end            
 
      end
       
       
      % while  abs(IQR{k}(end)-IQR{k}(end-1)) /sqrt(length(x)) > var_diff 
      while( abs(var(x)/var(x(1:(end-n_walks))) -1 ) > var_diff ) 
         %display(strcat(num2str(k), '_', num2str(VAR{k}(end)/VAR{k}(end-1)))); 
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
                if sum(abs(P(previous,:)))>0
                    current=min(find(cdf(previous,:)>aux));
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
%                 if (std(x)~=0)
%                     %Iqr=iqr(x)/std(x);
%                     VAR{k}=[VAR{k} var(x)];
%                 else
%                     %Iqr=0;
%                     VAR{k}=[VAR{k} 1];
%                 end
           end
         count=count+n_walks;  
         %IQR{k}=[IQR{k} Iqr];
        
      end
    if (std(x)~=0)
        %Iqr=iqr(x)/std(x);
        VAR{k}=[VAR{k} var(x)];
    else
        %Iqr=0;
        VAR{k}=[VAR{k} 1];
    end
      
   X(k)=mean(x);
   NWALKS(k)=count;
   end

   
elseif stat.adapt_cutoff==0 && stat_adapt_walks==1
    n_walks=stat.nwalks;
    max_step=stat.max_step;
    var_diff=stat.vardiff;

    X=zeros(size(b));
    NWALKS=zeros(size(b));
    IQR=cell(size(b,1),1);

   parfor k=1:size(b,1)
       count=0;
       x=[];
      
       
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
                if sum(abs(P(previous,:)))>0
                    current=min(find(cdf(previous,:)>aux));
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
                if (std(x)~=0)
                    Iqr=iqr(x)/std(x);
                else
                    Iqr=0;
                end
            end
            count=count+n_walks;  
            IQR{k}=[IQR{k} Iqr];
       end
           
       
       while   abs(IQR{k}(end)-IQR{k}(end-1)) > var_diff  
                               
             for walk=1:n_walks
                    estim=0;
                    previous=k;
                    W=1;
                    current=k;
                    estim=estim+W*b(current);
                    i=1;
                    while i<=max_step 
                    aux=rand;
                    if sum(abs(P(previous,:)))>0
                        current=min(find(cdf(previous,:)>aux));
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
                if (std(x)~=0)
                    Iqr=iqr(x)/std(x);
                else
                    Iqr=0;
                end
            end
            count=count+n_walks;  
            IQR{k}=[IQR{k} Iqr];
       end
        
   end
  
     
elseif stat.adapt_cutoff==1 && stat.adapt_walks==0   
    
    n_walks=stat.nwalks;
    max_step=stat.max_step;
    walkcut=stat.walkcut;
 
    X=zeros(size(b));
    NWALKS=zeros(size(b));
    IQR=cell(size(b,1),1);

   parfor k=1:size(b,1)
       count=0;
       x=[];

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
                if sum(abs(P(previous,:)))>0
                    current=min(find(cdf(previous,:)>aux));
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
           Iqr=iqr(x)/std(x);
       end
       count=count+n_walks;  
       IQR{k}=[IQR{k} Iqr];
   end
else
    n_walks=stat.nwalks;
    max_step=stat.max_step;
 
    X=zeros(size(b));
    NWALKS=zeros(size(b));
    IQR=cell(size(b,1),1);

   parfor k=1:size(b,1)
       count=0;
       x=[];

       for walk=1:n_walks
           estim=0;
           previous=k;
           W=1;
           current=k;
           estim=estim+W*b(current);
           i=1;
           while i<=max_step 
                aux=rand;
                if sum(abs(P(previous,:)))>0
                    current=min(find(cdf(previous,:)>aux));
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
                 Iqr=iqr(x)/std(x);
       end
    count=count+n_walks;  
    IQR{k}=[IQR{k} Iqr];
  end
   
end