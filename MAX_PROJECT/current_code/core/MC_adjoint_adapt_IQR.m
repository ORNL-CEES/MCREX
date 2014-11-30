function [x, IQR, Z, tally, res]=MC_adjoint_adapt_IQR(A, b, P, cdf, Pb, cdfb, stat)

Z=[];

x=[];
IQR=[];
res=[];

if stat.adapt_cutoff==1 && stat.adapt_walks==1
    
    n_walks=stat.nwalks;
    max_step=stat.max_step;
    walkcut=stat.walkcut;
    var_diff=stat.iqr;

    X=zeros(n_walks,size(b,1));
    tally=zeros(size(b));
    
    for j=1:2 
        for walk=1:n_walks
            aux=rand;

            %it detects what is the inital status of the chain
            previous=min(find(cdfb>aux));
            W=norm(b,1)*sign(b(previous));
            Wf=W*walkcut;
            tally(previous)=tally(previous)+1;
        %    W=b(previous)*length(find(P(previous,:)));
            X(walk,previous)=X(walk,previous)+W;
            i=1;
            while i<=max_step && W>Wf
                  aux=rand;
                  if sum(abs(P(previous,:)))>0
                    current=min(find(cdf(previous,:)>aux));
                    W=W*A(current,previous)/P(previous,current);
                  else
                     W=0;
                  end
                  if W==0
                      % if I reach a null state, I will never move away from
                      % there and the related permutation will not give any
                      % contribution to the estimation of the updating vector,
                      % thus this permutation can be neglected
                     break;
                  end

                  % I insert a contribution on the estimator related to the
                  % current state
                  X(walk,current)=X(walk,current)+W;
                  i=i+1;
                  previous=current;
            end
        end

        %computation of the expected value for the updating vector
        Z=[Z; X];
        x=mean(Z,1)';    
        res=[res norm( (b-(eye(size(A))-A)*x),2)/norm(b,2)];
        if isempty(find(std(Z)==0))
            IQR=[IQR max(iqr(Z)./std(Z))];
        else
            Iqr=iqr(Z)./std(Z);
            findout=find(std(Z)==0);
            Iqr(findout)=0;
            IQR=[IQR max(Iqr)*length(b)];
%             if(max(Iqr)==0 && ~isempty(find(Z)))
%                 IQR(end)=sum(IQR(1:end))*n_walks;
%             end
        end    
    end 
    
    
    while   abs(IQR(:,end) - IQR(:,end-1))> var_diff  
        
        %display(strcat('ratio= ', num2str(ratio)))
        
        X=zeros(n_walks,size(b,1));
        for walk=1:n_walks
                aux=rand;

                %it detects what is the inital status of the chain
                previous=min(find(cdfb>aux));
                W=norm(b,1)*sign(b(previous));
                Wf=W*walkcut;
                tally(previous)=tally(previous)+1;
            %    W=b(previous)*length(find(P(previous,:)));
                X(walk,previous)=X(walk,previous)+W;
                i=1;
                    while i<=max_step && W>Wf
                          aux=rand;
                          if sum(abs(P(previous,:)))>0
                            current=min(find(cdf(previous,:)>aux));
                            W=W*A(current,previous)/P(previous,current);
                          else
                             W=0;
                          end
                          if W==0
                              % if I reach a null state, I will never move away from
                              % there and the related permutation will not give any
                              % contribution to the estimation of the updating vector,
                              % thus this permutation can be neglected
                             break;
                          end

                          % I insert a contribution on the estimator related to the
                          % current state
                          X(walk,current)=X(walk,current)+W;
                          i=i+1;
                          previous=current;
                    end
         end

        %computation of the expected value for the updating vector
        Z=[Z; X];
        x=mean(Z,1)';    
        res=[res norm( (b-(eye(size(A))-A)*x),2)/norm(b,2)];
        if isempty(find(std(Z)==0))
            IQR=[IQR max(iqr(Z)./std(Z))];
        else
            Iqr=iqr(Z)./std(Z);
            findout=find(std(Z)==0);
            Iqr(findout)=0;
            IQR=[IQR max(Iqr)];
%             if(max(Iqr)==0 && ~isempty(find(Z)))
%                 IQR(end)=sum(IQR(1:end))*n_walks;
%             end
        end

    end
    
elseif stat.adapt_cutoff==0 && stat.adapt_walks==1
    
    n_walks=stat.nwalks;
    max_step=stat.max_step;
    var_diff=stat.iqr;
    
    X=zeros(n_walks,size(b,1));
    tally=zeros(size(b));
    
    for j=1:2
        for walk=1:n_walks
            aux=rand;

            %it detects what is the inital status of the chain
            previous=min(find(cdfb>aux));
            W=norm(b,1)*sign(b(previous));
            tally(previous)=tally(previous)+1;
        %    W=b(previous)*length(find(P(previous,:)));
            X(walk,previous)=X(walk,previous)+W;
            i=1;
            while i<=max_step 
                  aux=rand;
                  if sum(abs(P(previous,:)))>0
                    current=min(find(cdf(previous,:)>aux));
                    W=W*A(current,previous)/P(previous,current);
                  else
                     W=0;
                  end
                  if W==0
                      % if I reach a null state, I will never move away from
                      % there and the related permutation will not give any
                      % contribution to the estimation of the updating vector,
                      % thus this permutation can be neglected
                     break;
                  end

                  % I insert a contribution on the estimator related to the
                  % current state
                  X(walk,current)=X(walk,current)+W;
                  i=i+1;
                  previous=current;
            end
        end

        %computation of the expected value for the updating vector
        Z=[Z; X];
        x=mean(Z,1)';    
        res=[res norm( (b-(eye(size(A))-A)*x),2)/norm(b,2)];
        if isempty(find(std(Z)==0))
            IQR=[IQR max(iqr(Z)./std(Z))];
        else
            Iqr=iqr(Z)./std(Z);
            findout=find(std(Z)==0);
            Iqr(findout)=0;
            IQR=[IQR max(Iqr)];
        end    

    end

    while abs(IQR(:,end) - IQR(:,end-1))> var_diff        
        
        X=zeros(n_walks,size(b,1));
        for walk=1:n_walks
                aux=rand;

                %it detects what is the inital status of the chain
                previous=min(find(cdfb>aux));
                W=norm(b,1)*sign(b(previous));
                tally(previous)=tally(previous)+1;
            %    W=b(previous)*length(find(P(previous,:)));
                X(walk,previous)=X(walk,previous)+W;
                i=1;
                    while i<=max_step 
                          aux=rand;
                          if sum(abs(P(previous,:)))>0
                            current=min(find(cdf(previous,:)>aux));
                            W=W*A(current,previous)/P(previous,current);
                          else
                             W=0;
                          end
                          if W==0
                              % if I reach a null state, I will never move away from
                              % there and the related permutation will not give any
                              % contribution to the estimation of the updating vector,
                              % thus this permutation can be neglected
                             break;
                          end

                          % I insert a contribution on the estimator related to the
                          % current state
                          X(walk,current)=X(walk,current)+W;
                          i=i+1;
                          previous=current;
                    end
         end

        %computation of the expected value for the updating vector
        Z=[Z; X];
        x=mean(Z,1)';    
        res=[res norm( (b-(eye(size(A))-A)*x),2)/norm(b,2)];
        if isempty(find(std(Z)==0))
            IQR=[IQR max(iqr(Z)./std(Z))];
        else
            Iqr=iqr(Z)./std(Z);
            findout=find(std(Z)==0);
            Iqr(findout)=0;
            IQR=[IQR max(Iqr)];
        end    
        
    end

elseif  stat.adapt_cutoff==1 && stat.adapt_walks==0   
    
    n_walks=stat.nwalks;
    max_step=stat.max_step;
    walkcut=stat.walkcut;
    
    X=zeros(n_walks,size(b,1));
    tally=zeros(size(b));
    
    for walk=1:n_walks
        aux=rand;

        %it detects what is the inital status of the chain
        previous=min(find(cdfb>aux));
        W=norm(b,1)*sign(b(previous));
        tally(previous)=tally(previous)+1;
        Wf=W*walkcut;
    %    W=b(previous)*length(find(P(previous,:)));
        X(walk,previous)=X(walk,previous)+W;
        i=1;
        while i<=max_step && W>Wf
              aux=rand;
              if sum(abs(P(previous,:)))>0
                current=min(find(cdf(previous,:)>aux));
                W=W*A(current,previous)/P(previous,current);
              else
                 W=0;
              end
              if W==0
                  % if I reach a null state, I will never move away from
                  % there and the related permutation will not give any
                  % contribution to the estimation of the updating vector,
                  % thus this permutation can be neglected
                 break;
              end

              % I insert a contribution on the estimator related to the
              % current state
              X(walk,current)=X(walk,current)+W;
              i=i+1;
              previous=current;
        end
    end

    %computation of the expected value for the updating vector
    Z=[Z; X];
    x=mean(Z,1)';
    if isempty(find(std(Z)==0))
        IQR=[IQR max(iqr(Z)./std(Z))];
    else
        Iqr=iqr(Z)./std(Z);
        findout=find(std(Z)==0);
        Iqr(findout)=0;
        IQR=[IQR max(Iqr)];
    end    

else
    n_walks=stat.nwalks;
    max_step=stat.max_step;
    
    X=zeros(n_walks,size(b,1));
    tally=zeros(size(b));
    
    for walk=1:n_walks
        aux=rand;

        %it detects what is the inital status of the chain
        previous=min(find(cdfb>aux));
        W=norm(b,1)*sign(b(previous));
        tally(previous)=tally(previous)+1;
    %    W=b(previous)*length(find(P(previous,:)));
        X(walk,previous)=X(walk,previous)+W;
        i=1;
        while i<=max_step 
              aux=rand;
              if sum(abs(P(previous,:)))>0
                current=min(find(cdf(previous,:)>aux));
                W=W*A(current,previous)/P(previous,current);
              else
                 W=0;
              end
              if W==0
                  % if I reach a null state, I will never move away from
                  % there and the related permutation will not give any
                  % contribution to the estimation of the updating vector,
                  % thus this permutation can be neglected
                 break;
              end

              % I insert a contribution on the estimator related to the
              % current state
              X(walk,current)=X(walk,current)+W;
              i=i+1;
              previous=current;
        end
    end

    %computation of the expected value for the updating vector
    Z=[Z; X];
    x=mean(Z,1)';
    if isempty(find(std(Z)==0))
        IQR=[IQR max(iqr(Z)./std(Z))];
    else
        Iqr=iqr(Z)./std(Z);
        findout=find(std(Z)==0);
        Iqr(findout)=0;
        IQR=[IQR max(Iqr)];
    end    


end

end
    
