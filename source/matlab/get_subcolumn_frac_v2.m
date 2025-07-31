function [unique_col_frac ucol_num ucol_num_same] = get_subcolumn_frac_v2(nboxes, nlev, ncol, cc, overlap)

% this subroutinue is to assign cloud or clear to each sub-column of each 
% box by threecloud overlap schemes
%  overlap =1 maximum overlap
%  overlap =2 random overlap
%  overlap =3 maximum-random overlap
% nlev is the number of layers of cloud profiles
% ncol is the number of sub-column
% cc is the cloud fraction profile
%   output
%   subcol_frac  nboxes * ncol * nlev , the cloud assignments
%   unique_col_frac   nboxes * ucol_num * nlev, unique cloud assignments
%   ucol_num(nboxes)  the number of unique cloud assignments
%   ucol_num_same(nboxes,ucol_num)  the number of same cloud assignments for each unique cloud assignment
%   sum(ucol_num_same(nboxes,i), i=1:ucol_num) = ncol

    % add a clear layer over the top of atmosphere
   idd = find(cc<0.001);
   
    cc(idd) =0.0;
   
   tcc(:,2:nlev+1) = cc;
   tcc(:,1) = 0.0;

  for icol=1:ncol
      
       boxpos(:,icol) = (icol - 0.5)/ncol;
       
  end
  % initialize the sub-column with 0 (clear sky)
 
  subcol_frac = zeros(nboxes,ncol,nlev);
 
  cloud_threshold = zeros(nboxes,ncol);
  for ilev =1:nlev
     
       if ilev==1
           % maximum overlap
           if overlap==1
              
               cloud_threshold = boxpos;
             
           % random overlap or maximum and random overlap
           else
               cloud_threshold =rand(nboxes,ncol);
           end
       end
      ran =rand(nboxes,ncol);
      for icol=1:ncol
         % Maximum overlap
         if overlap ==1
              
               tmin(1:nboxes,icol)=0.0;
               flag(1:nboxes,icol) = 1;
              
         % Random overlap
         elseif overlap ==2
              
               tmin(1:nboxes,icol)=0.0;
               flag(1:nboxes,icol) = 0;
               
         % Maximum and Random overlap
         elseif overlap ==3
            for ibox=1:nboxes
                 tmin(ibox,icol)=min(tcc(ibox,ilev), tcc(ibox,ilev+1));
                 if cloud_threshold(ibox,icol)<min(tcc(ibox,ilev), tcc(ibox,ilev+1)) & cloud_threshold(ibox,icol)>0
                     flag(ibox,icol) = 1;
                 else
                     flag(ibox,icol) =0;
                 end
            end
         end
  
         %ran=get_ran(nboxes, seed);
  
          
         cloud_threshold(:,icol) = cloud_threshold(:,icol) .*flag(:,icol) + (1-flag(:,icol)) .* (tmin(:,icol) + (1 - tmin(:,icol)) .* ran(:,icol));
          
      end
      
      % cloud assignment
     for icol =1:ncol
          for ibox=1:nboxes
            if cc(ibox,ilev) > cloud_threshold(ibox,icol)
               subcol_frac(ibox, icol, ilev) =1;
            else
               subcol_frac(ibox, icol, ilev) =0;
            end
          end
      end
      
  end
  unique_col_frac= zeros(nboxes,ncol,nlev);
  for ibox=1:nboxes
      subcol = reshape(subcol_frac(ibox,:,:),ncol,nlev);
      [fracnew,pos1,pos2] = unique(subcol, 'rows');
      
      ucol_num(ibox) = length(pos1); % the number of unique cloud profiles
      
      unique_col_frac(ibox,1:ucol_num(ibox),:) = fracnew;
      
      % the number of same cloud profiles for each unique cloud profile
      for i=1:ucol_num(ibox)
          ucol_num_same(ibox,i) = length(find (ismember(subcol, fracnew(i,:),'rows') ==1 ));
      end
      %%%%%%%%%%%%%%%%%%%%%5
  %    frac =zeros(nlev,ncol+1)+NaN;  
 %  newlen = length(fracnew(:,1));
 
 %  icol=0
 %  for i=1:length(pos1)
 %      for j=1:ucol_num_same(ibox,i)
  %            icol=icol+1;
  %            for ilev =1:nlev
  %                   if fracnew(i,ilev)==1  
              
   %                      frac(ilev, icol) =fracnew(i,ilev);
  %                   end
   %           end
  %     end
  % end
  % for ilev=1:nlev
 %   for icol=1:ncol
 %       subcol_tmp(ilev,icol) = subcol(icol,ilev);
 %   end
%end
%       figure;
%   pcolor(1:ncol,1:nlev,subcol_tmp);
 %  set(gca,'Ydir','reverse');
%   setFontSize;
%   xlabel('Column Number');
%   ylabel('Model Level Number');
 %  Title('Maximum-Random Overlap');
 %  figure;
%   pcolor(1:ncol+1,1:nlev,frac);
%   set(gca,'Ydir','reverse');
%   setFontSize;
 %  xlabel('Column Number');
 %  ylabel('Model Level Number');
  % Title('Maximum-Random Overlap');
  % Title('Random Overlap');
  end
  
    
