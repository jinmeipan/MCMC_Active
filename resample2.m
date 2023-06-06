%% MCMC layer re-organizing
function profile=resample2(profile_in,n_out);

if(0)
profile_in=[1	267.7878966	0	275.7533998	0.15 0	0.278088604	1.5
2	266.8502055	0	245.2168551	0.15	0   0.278088604	1.5
3	265.5202771	0	199.02267	0.15	0   0.278088604	1.5
4	264.8803371	0	223.8651685	0.07	0   0.278088604	1.5
5	264.35	0	165.2743363	0.13	    0   0.09	0.25
6	265.1296117	0	109.0873786	0.19	0   0.09	0.25];

n_out=2;
end


n=size(profile_in,1); %current layers
profile=profile_in;   %the output

%            1             2     3       4         5        6     7   8
%profile=[ [1:1:nlayer]', m_T,  Zeros, m_dens, thickness, zeros, pex, davg];

dz0=profile(:,5); %units are cm, order is top-first to corespond with other memls_inputs


if n<n_out,
    
    dz1=dz0;
    while n<n_out,
        %Divide layers. The idea is to pick the largest layer, and split it
        %into two layers of equal size    
        %[~,ilargest]=max(dz1);
        ilargest=find( dz1 == max(dz1), 1, 'last');       
        
        [profile,n]=calcf_divide(profile,ilargest);
        dz1=profile(:,5);
    end

elseif n>n_out,
    
    dz1=dz0;
    while n>n_out,
%         %Combine layers. The idea is to pick the smallest layer, and
%         %combine it with the smaller adjacent layer.
%         %1) find the smallest layer        
%         [~,ismallest]=min(dz1);
%         %2) find the smallest adjacent layer to the smallest layer
%         iadj=[ismallest+1 ismallest-1];
%         jgood= iadj>0&iadj<=length(dz1);
%         iadj=iadj(jgood);
%         [~,iadjsmall]=min(dz1(iadj));
%         icombine=iadj(iadjsmall);
%         
%         [profile,n]=calcf_combine(profile,ismallest,icombine);
%         

        %Revised Sep21. Combine layers. The idea is to pick the two layers
        %with most similarity
        similarity=zeros(n-1,1);
        for i=1:n-1
            
            if(size(profile_in,2)>7)
                useful = [4,7,8];  %4=density,y=pex,8=dmax
            else
                useful = [4,7];
            end
            
            %idx
            temp = profile(:,useful);
            temp = sum(temp);
            idx = find(isnan(temp)==0);
            
            if(length(idx)==0)
                similarity(i) = nan;
            else
                %density and dmax
                prop1 = profile(i,useful(idx));
                prop2 = profile(i+1,useful(idx));
            
                %weights = [200, 1.0];
                %weighted by themselves;
                weights = mean([prop1;prop2]);
                weighted = (prop1 - prop2).^2./weights.^2;
                
                similarity(i) = sum(weighted);
            end
        end
        
%         disp('profile/similarity')
%         disp(num2str([profile,[similarity;nan]]))
        
        i1=find( similarity == min(similarity), 1, 'first');         
        %[~,i1]=min(similarity);
        i2=i1+1;
        
        %if all are nan values, use the smallest layer to combine
        if sum(isnan(similarity)) == n-1
            %[~,i1]=min(dz1);
            i1 = find( dz1 == min(dz1), 1, 'first');
            if(i1==n)
                i1=i1-1;
            end
            
            if(i1~=1)
                if(dz1(i1-1) < dz1(i1+1))
                    i1=i1-1;
                end
            end
            
            i2=i1+1;
        end
        
%         disp('chosen i')
%         disp(num2str(i1))

        
        
        [profile,n]=calcf_combine2(profile,i1,i2);
        dz1=profile(:,5);
    end        
else
    'nothing to do';
end


%rename layers
%profile(:,1)=[1:1:size(profile,1)]';


end


% 
% function [profile_out,n]=calcf_combine(profile_in,ismallest,icombine);
%    n=size(profile_in,1);
% %profile=[ [1:1:nlayer]', m_T,  Zeros, m_dens, thickness, pex, davg];
%     
%     dz1=profile_in(ismallest,5);
%     dz2=profile_in(icombine,5);
%     
%     dens1=profile_in(ismallest,4);
%     dens2=profile_in(icombine,4);
%     %Comb=(profile_in(ismallest,:)*dz1 + profile_in(icombine,:)*dz2)/(dz1+dz2);
%     Comb=(profile_in(ismallest,:)*dz1*dens1 + profile_in(icombine,:)*dz2*dens2) ...
%         /(dz1*dens1+dz2*dens2); %changed to mass weighted average
%     Comb(:,5)=dz1+dz2;
%     
%     min_i=min(ismallest,icombine);
%     max_i=max(ismallest,icombine);
%     
%     if ismallest==1 | icombine==1
%         profile_out=[Comb; profile_in(3:end,:)];
%     elseif ismallest==n | icombine==n
%         profile_out=[profile_in(1:end-2,:); Comb];
%     else
%         profile_out=[profile_in(1:min_i-1,:); Comb; profile_in(max_i+1:end,:)];
%     end
%     
%     n=size(profile_out,1);
%     profile_out(:,1)=[1:1:n]';
% end


function [profile_out,n]=calcf_combine2(profile_in,i1,i2)

    n=size(profile_in,1);
    
    %combine i1 and i2
    dz1=profile_in(i1,5);
    dz2=profile_in(i2,5);
    dens1=profile_in(i1,4);
    dens2=profile_in(i2,4);
    mass=dz1*dens1+dz2*dens2;
    
    if(isnan(mass)==1)
        dens1=1;
        dens2=1;
        mass=dz1+dz2; %use depth-weighted average
        
        if(isnan(mass)==1)
            %use average directly, if thickness is also missing
            dz1=1;
            dz2=1;
            mass=2;
        end
    end
            
    comb=(profile_in(i1,:)*dz1*dens1 + profile_in(i2,:)*dz2*dens2)/mass;
    
    %profile=[ [1:1:nlayer]', m_T,  Zeros, m_dens, thickness, pex, davg];
    comb(1,5)=dz1+dz2;
    
    if(n==2)
        profile_out=comb;
        n=n-1;
        profile_out(:,1)=[1:1:n]';
        return
    end
        
    if i1==1
        profile_out = [ comb; profile_in(i2+1:end,:)];
    else if i2==n
            profile_out = [ profile_in(1:i1-1,:);  comb];
        else
            profile_out = [ profile_in(1:i1-1,:);  comb; profile_in(i2+1:end,:)];
        end
    end
    n=n-1;
    profile_out(:,1)=[1:1:n]';
            
end



function [profile_out,n]=calcf_divide(profile_in,isplit)
    n=size(profile_in,1);
    
%     iall=1:n;
%            1             2     3       4         5        6     7
%profile=[ [1:1:nlayer]', m_T,  Zeros, m_dens, thickness, pex, davg]; 
    profile_in(isplit,5)=profile_in(isplit,5)/2;
    
    if isplit==1
        profile_out=[profile_in(1,:) ; profile_in(1,:); profile_in(2:end,:)];
    elseif isplit==n
        profile_out=[profile_in(1:end-1,:); profile_in(end,:) ; profile_in(end,:); ];
    else
        profile_out=[profile_in(1:isplit-1,:); ...
            profile_in(isplit,:); profile_in(isplit,:); ...
            profile_in(isplit+1:end,:)];
    end
    
    n=n+1;
    profile_out(:,1)=[1:1:n]';
end