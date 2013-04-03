function result = GARP(pin,qin,fuzz);

% GARP(p,q) returns a either 1 for Passes or 0 for Fails 
m = eye(size(pin,2));
for i = 1:size(pin,2);
    for j = 1:size(pin,2);
         m(i,j)=((pin(:,i)'*qin(:,i))>=(pin(:,i)'*qin(:,j) - fuzz));
    end;        
%   for i = 1:size(pin,2);
%         m(i,:)=(pin(:,i)'*qin(:,i)>=(pin(:,i)'*qin));
%   end;   
end;   
for i = 1:size(pin,2);                              
    for j=1:size(pin,2);                          
        m(i,:) = max([m(i,:); m(i,j)*m(j,:)]);
    end;                                    
end;                                        
result=1;
for i = 1:size(pin,2);
    for j = 1:size(pin,2);
        if m(i,j)==1; 
            if (pin(:,j)'*qin(:,j))>(pin(:,j)'*qin(:,i) + fuzz); 
                result=0; break;
            end;
        end;
    end;
    if result==0; break; end;
end;
    