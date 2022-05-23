function [D,Idx] = pdist2ties(X,Y,Distance,additionalArg,K,includeTies)
%#codegen
%Helper codegen function for calling pdist2 with passed arguments and
%   handling ties if includeTies is true.

%   Copyright 2017-2019 The MathWorks, Inc.

Kinit = K;
coder.internal.prefer_const(includeTies,Distance,K,additionalArg);
[nx,~] = size(X);
[ny,~] = size(Y);
if includeTies
    if (K < 3)
        extra = 2;     
    elseif K < 5
        extra = K;
    else
        extra = 5;
    end
    K = K + extra;
    
    if isempty(additionalArg)
        [distT, idxT] = pdist2(X,Y, Distance, 'smallest',K);
    else
        [distT, idxT] = pdist2(X,Y, Distance, additionalArg, 'smallest',K);
    end
    
    %  Degenerate case, just return an empty cell of the proper size.
    if (nx == 0 || ny == 0)
        D = repmat({zeros(1,0)},coder.internal.indexInt(ny),1);
        Idx = repmat({zeros(1,0)},coder.internal.indexInt(ny),1);
        return;
    end

    D = coder.nullcopy(cell(coder.internal.indexInt(ny),1));
    Idx = coder.nullcopy(cell(coder.internal.indexInt(ny),1));    

    if Kinit >= nx 
        for ii = 1 : coder.internal.indexInt(ny)
            Idx{ii} = idxT(:,ii)';
            D{ii} = distT(:,ii)';
        end
        return; 
    end
    
    %Expect that most query points don't have extra tie neighbors
    %handle no ties first
    notDone = true(coder.internal.indexInt(ny),1);
    allDone = true;
    % first pass, assign no ties
    if coder.internal.isConst(ny)
        dounroll = ny < 10;
    else
        dounroll = false;
    end
    for ii = coder.unroll(1:coder.internal.indexInt(ny), dounroll) 
        notDone(ii) = ~(distT(coder.internal.indexInt(Kinit)+1,ii) > distT(coder.internal.indexInt(Kinit),ii) ||...
                        isnan(distT(coder.internal.indexInt(Kinit)+1,ii)));
        allDone = allDone && ~notDone(ii); 
        if ~notDone(ii)
            Idx{ii} = idxT(1:coder.internal.indexInt(Kinit),ii)';
            D{ii} = distT(1:coder.internal.indexInt(Kinit),ii)'; 
        end
    end   
                  
    if allDone
       return;
    end
    
    %second pass
    if K >= nx
        for ii = 1:coder.internal.indexInt(ny)
            if notDone(ii)
                tempIdx = coder.internal.indexInt(-1);
                for jj = coder.internal.indexInt(Kinit+1): coder.internal.indexInt(K)
                    if distT(jj,ii) > distT(coder.internal.indexInt(Kinit),ii) 
                        tempIdx = jj;
                        break;
                    end
                end
                if tempIdx == -1
                    for jj = coder.internal.indexInt(Kinit+1): coder.internal.indexInt(K)
                        if isnan(distT(jj,ii)) 
                            tempIdx = jj;
                            break;
                        end
                    end
                    if tempIdx == -1
                       N = coder.internal.indexInt(coder.internal.indexInt(K)-coder.internal.indexInt(Kinit));
                    else
                       N = coder.internal.indexInt(tempIdx-1); 
                    end    
                else
                    N = coder.internal.indexInt(tempIdx-1);
                end
                Idx{ii} = idxT(1:N,ii)';
                D{ii} = distT(1:N,ii)'; 
            end
        end
    else
        for ii = 1:coder.internal.indexInt(ny)
            if notDone(ii)
                tempIdx = coder.internal.indexInt(-1);
                for jj = coder.internal.indexInt(Kinit+1): coder.internal.indexInt(K)
                    if distT(jj,ii) > distT(coder.internal.indexInt(Kinit),ii) 
                        tempIdx = jj;
                        break;
                    end
                end
                if tempIdx == -1
                    if distT(coder.internal.indexInt(K),ii) == coder.internal.inf
                        r = distT(coder.internal.indexInt(K),ii);
                    else
                        r = distT(coder.internal.indexInt(K),ii)+10*eps(distT(coder.internal.indexInt(K),ii));
                    end
                    if ~isempty(additionalArg)
                        [dTemp, iTemp] = pdist2(X,Y(ii,:),Distance,additionalArg,'Radius',r);
                    else
                        [dTemp, iTemp] = pdist2(X,Y(ii,:),Distance,'Radius',r);
                    end
                    Idx{ii} = iTemp{:};
                    D{ii} = dTemp{:};  
                else
                    N = coder.internal.indexInt(tempIdx-1); 
                    Idx{ii} = idxT(1:N,ii)';
                    D{ii} = distT(1:N,ii)';                     
                end
            end
        end        
    end
else  
  if ~isempty(additionalArg)
     [Dout, IdxOut] = pdist2(X,Y,Distance,additionalArg,'Smallest',K);
  else
     [Dout, IdxOut] = pdist2(X,Y,Distance,'Smallest',K);  
  end
  D = Dout';
  Idx = IdxOut';
end
end