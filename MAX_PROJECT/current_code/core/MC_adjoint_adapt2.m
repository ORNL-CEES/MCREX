function [x, y, Z, tally, res, reject]=MC_adjoint_adapt2(A, b, P, cdf, Pb, cdfb, stat)

Z=[];

x=[];
y=[];
res=[];
reject=0;

if stat.adapt_cutoff==1 && stat.adapt_walks==1
    
    n_walks=stat.nwalks;
    max_step=stat.max_step;
    walkcut=stat.walkcut;
    var_cut=stat.varcut;
    var_diff=stat.vardiff;

    X=zeros(n_walks,size(b,1));
    tally=zeros(size(b));
    
    for j=1:2 
        for walk=1:n_walks
            aux=rand;

            %it detects what is the inital status of the chain
            previous=min(find(cdfb>aux));
            W=sum(abs(b))*sign(b(previous));
            Wf=W*walkcut;
            tally(previous)=tally(previous)+1;
        %    W=b(previous)*length(find(P(previous,:)));
            X(walk,previous)=X(walk,previous)+W;
            i=1;
            while i<=max_step && W>Wf
                  aux=rand;
                  if sum(abs(P(previous,:)))==1
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
        Y=Z.^2;
        y=[ y sqrt((mean(Y,1)'-(x.^2))./(size(Z,1))) ];     
        res=[res norm( (b-(eye(size(A))-A)*x),2)/norm(b,2)];
    end 
    
    
    while norm(y(:,end), 1)/norm(x, 1) > var_cut  || norm(y(:,end) - y(:,end-1), 1)/norm(y(:,1)) > var_diff
        
        if norm(y(end), 1)/norm(x(end), 1) < var_cut
            reject=reject+1;
        end        
        
        X=zeros(n_walks,size(b,1));
        for walk=1:n_walks
                aux=rand;

                %it detects what is the inital status of the chain
                previous=min(find(cdfb>aux));
                W=sum(abs(b))*sign(b(previous));
                Wf=W*walkcut;
                tally(previous)=tally(previous)+1;
            %    W=b(previous)*length(find(P(previous,:)));
                X(walk,previous)=X(walk,previous)+W;
                i=1;
                    while i<=max_step && W>Wf
                          aux=rand;
                          if sum(abs(P(previous,:)))==1
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
        Y=Z.^2;
        y=[ y sqrt((mean(Y,1)'-(x.^2))./(size(Z,1))) ]; 
        res=[res norm( (b-(eye(size(A))-A)*x),2)/norm(b,2)];
    end
    
elseif stat.adapt_cutoff==0 && stat.adapt_walks==1
    
    n_walks=stat.nwalks;
    max_step=stat.max_step;
    var_cut=stat.varcut;
    var_diff=stat.vardiff;
    
    X=zeros(n_walks,size(b,1));
    tally=zeros(size(b));
    
    for j=1:2
        for walk=1:n_walks
            aux=rand;

            %it detects what is the inital status of the chain
            previous=min(find(cdfb>aux));
            W=sum(abs(b))*sign(b(previous));
            tally(previous)=tally(previous)+1;
        %    W=b(previous)*length(find(P(previous,:)));
            X(walk,previous)=X(walk,previous)+W;
            i=1;
            while i<=max_step 
                  aux=rand;
                  if sum(abs(P(previous,:)))==1
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
        Y=Z.^2;
        y=[y sqrt((mean(Y,1)'-(x.^2))./(size(Z,1)))]; 
        res=[res norm( (b-(eye(size(A))-A)*x),2)/norm(b,2)];
    end

    while norm(y(:,end), 1)/norm(x, 1) > var_cut  || norm(y(:,end) - y(:,end-1), 1)/norm(y(:,end)) > var_diff
        
        if norm(y(end), 1)/norm(x(end), 1) < var_cut
            reject=reject+1;
        end                
        
        X=zeros(n_walks,size(b,1));
        for walk=1:n_walks
                aux=rand;

                %it detects what is the inital status of the chain
                previous=min(find(cdfb>aux));
                W=sum(abs(b))*sign(b(previous));
                tally(previous)=tally(previous)+1;
            %    W=b(previous)*length(find(P(previous,:)));
                X(walk,previous)=X(walk,previous)+W;
                i=1;
                    while i<=max_step 
                          aux=rand;
                          if sum(abs(P(previous,:)))==1
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
        Y=Z.^2;
        y=[ y sqrt((mean(Y,1)'-(x.^2))./(size(Z,1))) ];   
        res=[res norm( (b-(eye(size(A))-A)*x),2)/norm(b,2)];
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
        W=sum(abs(b))*sign(b(previous));
        tally(previous)=tally(previous)+1;
        Wf=W*walkcut;
    %    W=b(previous)*length(find(P(previous,:)));
        X(walk,previous)=X(walk,previous)+W;
        i=1;
        while i<=max_step && W>Wf
              aux=rand;
              if sum(abs(P(previous,:)))==1
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
    Y=Z.^2;
    y=sqrt((mean(Y,1)'-(x.^2))./(size(Z,1))); 
    
else
    n_walks=stat.nwalks;
    max_step=stat.max_step;
    
    X=zeros(n_walks,size(b,1));
    tally=zeros(size(b));
    
    for walk=1:n_walks
        aux=rand;

        %it detects what is the inital status of the chain
        previous=min(find(cdfb>aux));
        W=sum(abs(b))*sign(b(previous));
        tally(previous)=tally(previous)+1;
    %    W=b(previous)*length(find(P(previous,:)));
        X(walk,previous)=X(walk,previous)+W;
        i=1;
        while i<=max_step 
              aux=rand;
              if sum(abs(P(previous,:)))==1
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
    Y=Z.^2;
    y=sqrt((mean(Y,1)'-(x.^2))./(size(Z,1))); 

end

end
    
