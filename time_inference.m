
% mat - matrix with CyTOF measures (cells are in the rowa and markers in the columns)
% marker_names - vector with the marker names corresponding to the columns in mat

load('mouse_data.mat');

mat=sessionData; 
marker_names=gates{3};
marker_names=marker_names';
define_dt = 1; % time in hours post IDU injection (1,4,6 or 12 hours)

% remove non T-cells: 
cofactor=5;  mat=asinh(mat/cofactor);
non_T_cells=mat(:,38); % 161Dy_CD11b_CD11c_CD19 markers
k2=find(non_T_cells<=5); 
mat=mat(k2,:); 
mat=cofactor*sinh(mat);

% remove outliers (according to specific markers)  
markers= [3,18,23,27,28,29,30,33,37,41,43,48,49,31,51,52]; 
remove_cells=remove_outliers(mat(:,markers) , percentile); %%% romove cells above 99%
mat(remove_cells,:)=[]; 
mat2=mat;

% normalize 
mat=bsxfun(@rdivide, mat, mean(mat));

% Insert your method for pseudotime inference here  
% Trajectory (output of the pseudotime inference) : a vector with the segment index assigned to each cell (if necessary group cells into segments)  
% ell : the spatial coordinated in CD4/CD8 space of the centers of each segment  
% Segment_index : if the trajectory bifurcates this variable containes the index of each portion of the trajectory (e.g. index 1 for segments before the bifurcation, and indices 2 and 3 for segments after the bifurcation)   
% gamma1_thresh and gamma2_thresh are the respective thesholds for division and apoptosis (obtained using the approach shown in Figure3)


% calculate the mean marker values in each segment 
gamma1=mat(:,12); % IDU 
gamma2=mat(:,20); % Caspase3
GM1=zeros(1,max(Trajectory)); %% average gamma1 per segment 
GM2=zeros(1,max(Trajectory)); %% average gamma2 per segment
mean_marker = zeros( length(marker_names) , max(Trajectory) );

 for J=1:max(Trajectory)
     
     cells_in_current_segment = find(Trajectory==J);
     mean_marker(:,J)         = mean(mat(cells_in_current_segment,:)); 
    
     GM1(J)=length(cells_in_current_segment(ismember(cells_in_current_segment,gamma1_thresh)))/length(cells_in_current_segment); % percentage of cells within the segment that are dividing
     if ~isempty(cells_in_current_segment(ismember(cells_in_current_segment,gamma1_thresh))) 
       GM2(J)=length(cells_in_current_segment(ismember(cells_in_current_segment,gamma2_thresh)))/length(cells_in_current_segment(ismember(cells_in_current_segment,gamma1_thresh))); % percentage of apoptosis cells is calculated out of the dividing population 
     else
       GM2(J)=0;   
     end

 end


% time in each segment (breaking the calculation into 2 steps - first calculate time in the CD4/CD8 space) 

choose_segment=1; 
t_numeric_cal_DN_DP  = return_t_segment(1:max(Trajectory),Segment_index ,Trajectory,GM1,GM2,mat,choose_segment,define_dt);

choose_segment=2; 
t_numeric_cal_CD4    = return_t_segment(1:max(Trajectory),Segment_index ,Trajectory,GM1,GM2,mat,choose_segment,define_dt);

choose_segment=3; 
t_numeric_cal_CD8    = return_t_segment(1:max(Trajectory),Segment_index ,Trajectory,GM1,GM2,mat,choose_segment,define_dt);



% sub-divide segments using higher dimentions 

choose_segment=1; 
t_numeric_cal_DN_DP_NEW  = return_t_segment_ND(mat,mat2,1:max(Trajectory),Segment_index,Trajectory,ell, choose_segment,t_numeric_cal_DN_DP  ,define_dt , gamma1_thresh, gamma2_thresh );

choose_segment=2; 
t_numeric_cal_CD4_NEW    = return_t_segment_ND(mat,mat2,1:max(Trajectory),Segment_index,Trajectory,ell, choose_segment,t_numeric_cal_DN_DP  ,define_dt , gamma1_thresh, gamma2_thresh );

choose_segment=3; 
t_numeric_cal_CD8_NEW    = return_t_segment_ND(mat,mat2,1:max(Trajectory),Segment_index,Trajectory,ell, choose_segment,t_numeric_cal_DN_DP  ,define_dt , gamma1_thresh, gamma2_thresh );


time_CD4_trajectory_NEW=[t_numeric_cal_DN_DP_NEW  t_numeric_cal_DN_DP_NEW(end)+t_numeric_cal_CD4_NEW];
time_CD8_trajectory_NEW=[t_numeric_cal_DN_DP_NEW  t_numeric_cal_DN_DP_NEW(end)+t_numeric_cal_CD8_NEW];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Functions: 

function ignore_cells=remove_outliers(DATA_mat, percentile)
 P= prctile(DATA_mat, percentile); 
 DATA_mat=bsxfun(@rdivide, DATA_mat, P);
 [ignore_cells,~]=find(DATA_mat>=1);
 ignore_cells=unique(ignore_cells);
 return
end



function t_numeric_cal = return_t_segment(ind,Segment_index,Trajectory,GM1,GM2,mat, choose_segment,define_dt)

ind=ind(Segment_index==choose_segment); 
n=length(ind); 
[k1,k2] = size(mat);

% probability densities
 ff=zeros(1,n);
 for II=1:n
    D=find(A==ind(II));
    ff(II)=length(D)/k1;  % proportion of cells corrisponding to ell(II)  
 end
FF=cumsum(ff);
 
dt=define_dt/24; 
GM1=GM1/dt;
GM2=GM2/dt;

if choose_segment==1  %%% DN+DP segment ends now at T=20
  [t_numeric_cal, t_numeric_initial_current ,GM2] = return_t_numeric(GM1,GM2,ff,n); 
  if min(diff(t_numeric_cal))<0 %%% if "return_t_numeric" failed
  [t_numeric_cal, t_numeric_initial_current] = return_t_numeric2(GM1,GM2,ff,n); 
  end     
else  %%% DN+DP segment ends now at T=5
  [t_numeric_cal, t_numeric_initial_current] = return_t_numeric_SP(GM1,GM2,ff,n); 
  if min(diff(t_numeric_cal))<0 %%% if "return_t_numeric" failed
  [t_numeric_cal, t_numeric_initial_current] = return_t_numeric2_SP(GM1,GM2,ff,n); 
  end  
end

return

end


function [t_numeric_cal, t_numeric_initial_current , GM2] = return_t_numeric(GM1,GM2,ff,n) 

T=20; %%% time cells spend in the thymus 
t2 = 0.05;    

t_numeric_cal=zeros(1,n);     
t_numeric_cal(2)= t2;    
for i=3: (n+1)
    
    param = 0.1; 
    first_comp(i) = ff(i-2)/(t_numeric_cal(i-1)-t_numeric_cal(i-2));
    second_comp(i) = (GM1(i-2)-GM2(i-2))*ff(i-2); 
     
    if i>42
         GM2(i-2) = GM1(i-2) - (param-first_comp(i))/ff(i-2);
    end
     
    t_numeric_cal(i)=t_numeric_cal(i-1)+ ff(i-1)/(  ff(i-2)/(t_numeric_cal(i-1)-t_numeric_cal(i-2))+(GM1(i-2)-GM2(i-2))*ff(i-2)  );
     
end
t_numeric_cal=(t_numeric_cal(1:end-1)+t_numeric_cal(2:end))/2; 
   
t_numeric_cal = (t_numeric_cal / max(t_numeric_cal))*T;

t_numeric_initial_current = t2; 

return

end



function [t_numeric_cal, t_numeric_initial_current] = return_t_numeric2(GM1,GM2,ff,n) 

T=10; %%% time cells spend in the thymus 
deviation_tolerance=0.1; %%% allow a tolerance so that T can equal T plus/minus deviation_tolerance
limit1=T-deviation_tolerance; limit2=T+deviation_tolerance;
t_initial_lim1=0.0244;  
t_initial_lim2=0.12;
t_initial_lim=linspace(t_initial_lim1,t_initial_lim2,10000); 
t_numeric_cal=zeros(1,n);
j=1;

while j<=10000 && t_initial_lim1<t_initial_lim2 &&  all(t_numeric_cal(end)<=limit1 || t_numeric_cal(end)>=limit2)
    
t_numeric_cal=zeros(1,n);     
t_numeric_cal(2)= t_initial_lim(j); 
   for i=3: (n+1)
     t_numeric_cal(i)=t_numeric_cal(i-1)+ ff(i-1)/(  ff(i-2)/(t_numeric_cal(i-1)-t_numeric_cal(i-2))+(GM1(i-2)-GM2(i-2))*ff(i-2)  );
   end
   t_numeric_cal=(t_numeric_cal(1:end-1)+t_numeric_cal(2:end))/2; 
t_numeric_initial_current= t_initial_lim(j);   
j=j+1; 
end 

return

end




function [t_numeric_cal, t_numeric_initial_current] = return_t_numeric_SP(GM1,GM2,ff,n) 

T=5; %%% time cells spend in the thymus 
t2 = 0.05;   

t_numeric_cal=zeros(1,n);     
t_numeric_cal(2)= t2;    
  
for i=3: (n+1)
     t_numeric_cal(i)=t_numeric_cal(i-1)+ ff(i-1)/(  ff(i-2)/(t_numeric_cal(i-1)-t_numeric_cal(i-2))+(GM1(i-2)-GM2(i-2))*ff(i-2)  );
end

t_numeric_cal=(t_numeric_cal(1:end-1)+t_numeric_cal(2:end))/2;
   
t_numeric_cal = (t_numeric_cal / max(t_numeric_cal))*T;

t_numeric_initial_current =t2;   

return

end



function [t_numeric_cal, t_numeric_initial_current] = return_t_numeric2_SP(GM1,GM2,ff,n) 

T=4; %%% time cells spend in the thymus 
deviation_tolerance=0.1; %%% allow a tolerance so that T can equal T plus/minus deviation_tolerance
limit1=T-deviation_tolerance; limit2=T+deviation_tolerance;
t_initial_lim1=0.0244;  
t_initial_lim2=0.12;
t_initial_lim=linspace(t_initial_lim1,t_initial_lim2,10000); 
t_numeric_cal=zeros(1,n);
j=1; 

while j<=10000 && t_initial_lim1<t_initial_lim2 &&  all(t_numeric_cal(end)<=limit1 || t_numeric_cal(end)>=limit2)
    
t_numeric_cal=zeros(1,n);     
t_numeric_cal(2)= t_initial_lim(j); 
   for i=3: (n+1)
     t_numeric_cal(i)=t_numeric_cal(i-1)+ ff(i-1)/(  ff(i-2)/(t_numeric_cal(i-1)-t_numeric_cal(i-2))+(GM1(i-2)-GM2(i-2))*ff(i-2)  );
   end
   t_numeric_cal=(t_numeric_cal(1:end-1)+t_numeric_cal(2:end))/2; 
t_numeric_initial_current= t_initial_lim(j);   
j=j+1; 
end 

return

end


function inds_between_points=return_inds_between(ell1,ell2,markerA,markerB,ind_cells_current)

dist=zeros(length(markerA),2); %%% will contain distance between each cell in pool and both points (ell1,ell2) 

for I=1:size(dist,1)

    dist(I,1)=sqrt( (ell1(1,1)-markerA(I)).^2 +  (ell1(1,2)-markerB(I)).^2 ); % dist between point and ell1
    dist(I,2)=sqrt( (ell2(1,1)-markerA(I)).^2 +  (ell2(1,2)-markerB(I)).^2 ); % dist between point and ell2

end
    
inds_between_points= (dist(:,1)-dist(:,2))<0 ; %% cells closer to ell1
inds_between_points= ind_cells_current(inds_between_points);

return

end



function [ell_cell , ell_cell_ind]= cell_closest_to_ell2(mat,ell)

ell_cell=zeros(size(ell,1),  size(mat,2));
ell_cell_ind=zeros(size(ell,1),1);
markerA=mat(:,1); 
markerB=mat(:,2); 

D=sqrt( (ell(1,1) - ell(2,1))^2 + (ell(1,2) - ell(2,2))^2  ); %%% distance between 2 points in ell

%%% distances from first point
dist1=sqrt( (ell(1,1)-markerA).^2 +  (ell(1,2)-markerB).^2 );
dist2=sqrt( (ell(2,1)-markerA).^2 +  (ell(2,2)-markerB).^2 );

indxs=[]; 
p=0.25*D; %% radius around first point 

while isempty(indxs) && p < 0.5*D %%% so cells closest to 2 points on ell don't overlap 
 indxs = find(dist1<=p); 
 p=p*1.1;    
end
    
if isempty(indxs) %%% choose cells further away out of radius 
    p=0.5*D;
    k=find(dist2>p); %% cells out of second radius
    [n,m]=sort(dist1(k)); 
    indxs=k(m); 
    if length(indxs)>=5
      indxs=indxs(1:5); 
    end
end

cell_reference=[ell(1,1) ell(1,2)  mean(mat(indxs,3:end) , 1)]; 
mat_reference=mat(indxs,:); 

R=sum(bsxfun(@minus, mat_reference, cell_reference).^2,2);
ind=find(R==min(R));
ell_cell(1,:)=mat_reference(ind(1),:);
ell_cell_ind(1)=indxs(ind(1));

indxs1=indxs;

indxs=[]; 
p=0.25*D; %% radius around first point 

while isempty(indxs) && p < 0.5*D %%% so cells closest to 2 points on ell don't overlap 
 indxs = find(dist2<=p); 
 p=p*1.1;    
end

if isempty(indxs) %%% choose cells further away out of radius 
    p=0.5*D;
    k=find(dist1>p); %% cells out of second radius
    [n,m]=sort(dist2(k)); 
    indxs=k(m); 
    if length(indxs)>=5
      indxs=indxs(1:5); 
    end
end

cell_reference=[ell(2,1) ell(2,2)  mean(mat(indxs,3:end) , 1)]; 
mat_reference=mat(indxs,:); 

R=sum(bsxfun(@minus, mat_reference, cell_reference).^2,2);
ind=find(R==min(R));
ell_cell(2,:)=mat_reference(ind(1),:);
ell_cell_ind(2)=indxs(ind(1));

indxs2=indxs;
return

end




function trajectory = calculate_trajectory(G,anchors) 

trajectory = anchors(1,:); 
cur = trajectory;  

indxs = repmat({[-1 0 1]},1,2);
Y0=cell(1,2);
[Y0{:}]=ndgrid(indxs{:});
Y1=cell2mat(Y0); 
vicinity_mat2=reshape(Y1,[3^2,2]); 
vicinity_mat=zeros(3^ndims(G),ndims(G));
vicinity_mat(:,1:2)=repmat(vicinity_mat2,3^(ndims(G)-2) ,1); 
div_factor=3^ndims(G); 

for ii=ndims(G):-1:3
    
    rep_factor=length(vicinity_mat)/div_factor; %%% number of repeats of -1,0,1
    temp= reshape(repmat([-1 0 1], div_factor/3,1),1, div_factor );
    vicinity_mat(:,ii)= repmat(temp,1,rep_factor);
    div_factor=div_factor/3; 
end

prev_pos=[];   no_progression=[]; 
for I=1:size(anchors,1)-1 
 
 next_anchor = anchors(I+1,:);
 vicinity_curr = bsxfun(@plus, vicinity_mat, cur); 
 
 while any(abs(cur - next_anchor) > 1)
  
     vicinity_curr = bsxfun(@plus, vicinity_mat, cur); 
     [n,m]=find(vicinity_curr==0); 
     n1=unique(n);     
     q=size(G); 
     vicinity_temp=bsxfun(@minus, vicinity_curr, q);
     [n,m]=find(vicinity_temp>0);  
     n2=unique(n);
     ignore_position =[n1;n2]; 

     K=num2cell(vicinity_curr); 
     nghbrs=zeros(1,length(vicinity_curr));
      for j=1:length(vicinity_curr)
         if ~isempty(find(ignore_position==j, 1))
            nghbrs(j)=-inf;
         else
            nghbrs(j)=G(K{j,:});
         end
      end
     nghbrs=reshape(nghbrs,3*ones(1,ndims(G))); 
      
     cur_cell=num2cell(cur); cur_val=G(cur_cell{:}); 
     
     if ~isempty(prev_pos)
       [nghbrs,no_progression] = prevent_progression(G,prev_pos,cur,vicinity_prev,vicinity_curr, nghbrs);  
     end 
     
     [~,next_step] = max(reshape(nghbrs, 1, []));
     next_step1=cell(1,ndims(G));
     [next_step1{:}] = ind2sub(size(nghbrs), next_step);
     next_step=cell2mat(next_step1)-2;
     
     step_size = sort(nghbrs(:) - cur_val); 
     a=num2cell(cur + sign(next_anchor - cur)); 
     step_size_moment =G(a{:})-cur_val; 
     
     if  ~isempty(no_progression) || ~any(next_step)  || (step_size(end)-step_size_moment < 5000) 
         next_step = sign(next_anchor - cur);
     end
     
    prev_pos=cur; 
    vicinity_prev=vicinity_curr;
 
    cur = cur + next_step;
    trajectory(end+1,:) = cur;
 end
 prev_pos=cur; vicinity_prev=vicinity_curr; 
cur=next_anchor;
trajectory(end+1,:) = cur;  

end
return 
end



function trajectory_val=return_trajec_val(DATA_mat,trajectory,loc)

trajectory_val=zeros(size(trajectory));
A=max(loc); %%%% for each of the 'm' markers (theoreticly each marker can be divided to a grid of different size)
[~,m]=size(DATA_mat);
for I=1:m
   
   B= linspace(min(DATA_mat(:,I)),max(DATA_mat(:,I)),A(I)); 
   trajectory_val(:,I)=B(trajectory(:,I));  
   
end
return
end



function A_new=re_calculate_Trajectory(ell_new,size_segments,segmented_original_cells,ell_index_original,mat,markers_segment)

A_new=zeros(size(segmented_original_cells));
add_to_A=0;

for I=1: length(size_segments)
   
    ell_cur=ell_new( 1:size_segments(I)  ,:);
    curr_cells_in_segment=segmented_original_cells(ell_index_original==I); 
    markers_expand=markers_segment(I); 
    DATA_mat=mat(curr_cells_in_segment,markers_expand); 
    
    %%%% assign each point in current pool to new small ell
    A_curr=zeros(1,size(DATA_mat,1)); 
    for I2=1:length(A_curr)   
      distance_I=0;
       for J=1:size(DATA_mat,2)
         distance_I=distance_I+ (DATA_mat(I2,J)-ell_cur(:,J)).^2; 
       end
      distance_I=sqrt(distance_I); 
      A_temp=find(distance_I==min(distance_I));
      A_curr(I2)=A_temp(1); 
    end
    
     A_new(curr_cells_in_segment) = A_curr + add_to_A; 
     add_to_A = add_to_A+ size_segments(I)-1;
     
    ell_new(1:size_segments(I) - 1,:)=[];
    
end

return

end



 function t_numeric_cal  = return_t_segment_ND(mat,mat2,ind,Segment_index,Trajectory,ell, choose_segment,t_numeric_cal_original ,define_dt , gamma1_thresh, gamma2_thresh )

ind=ind(Segment_index==choose_segment); 
ind_cells_curr=ismember(Trajectory,ind); %%% current cells assigned to chosen ell 
Trajectory=Trajectory(ismember(Trajectory,ind));
[length_original,k2] = size(mat);

markers = mat(ind_cells_curr,:); 

grid_vec=diff(t_numeric_cal_original); grid_vec=grid_vec/max(grid_vec); 
grid_max=10; 
grid_vec=grid_vec*grid_max; grid_vec=ceil(grid_vec); 

markers_test=[3,18,23  ,27,28,29,30,33,37,  41,43,48   ,51,52]; %% without 24 and 49 which will be included anyway

markers_segment=zeros(max(Trajectory)-1 , 8); % for using different markers each segment  
B=ones(size(Trajectory));

ell_new=[];
ell_index_original=[];
segmented_original_cells=[]; 

Trajectory_new=zeros(size(Trajectory)); 
add_to_Trajectory=0; 
size_segments=zeros(1,max(Trajectory)-1); 

for I=1:max(Trajectory)-1    
    
    ind_cells_current=find( (Trajectory==ind(I) | Trajectory==ind(I+1)) & B==1 ); %%% Pool of cells closest to 2 points 

   if I <= (max(Trajectory)-2)
    inds_between_points=return_inds_between(ell(I,:),ell(I+2,:),markerA(ind_cells_current),markerB(ind_cells_current),ind_cells_current); 
    B(inds_between_points)=0; 
   else
    inds_between_points=find(B==1); %% remaining cells    
   end
    
    %%% choose markest for current segment with maximum variance 
    mat_cur=mat(inds_between_points,markers_test); 
    var_cur=std(mat_cur)./mean(mat_cur); 
    [~ , marker_test_inds]=sort(var_cur); 
    markers_expand=markers_test( marker_test_inds(end-5:end) ); % markers with max variance  
    markers_segment(I,:)=[49,24,   markers_expand]; % add CD4 and CD8 
    markers_expand=[49,24,   markers_expand]; 
    
    %%% Expand ell
    DATA_mat=mat(inds_between_points,markers_expand); 
    [~,m]=size(DATA_mat);
    grid_size=grid_vec(I);  %10    grid_vec(I)
    J1= num2cell(grid_size*ones(1,m));
    J1{1}=20; J1{2}=20; %% for CD4 and CD8 
    [G edges mid loc] = histcn(DATA_mat,J1{:});
    

    [ell_cell_curr , ell_cell_ind_curr]= cell_closest_to_ell2(mat(inds_between_points,markers_expand), [ell(I,:) ; ell(I+1,:)] );

    start_loc=loc(ell_cell_ind_curr(1),:); %%% location of ell(I) on grid (starting point) 
    end_loc=loc(ell_cell_ind_curr(2),:);  %%% location of ell(I+1) on grid (ending point) 
    anchors = [start_loc ; end_loc];   
    
    trajectory_curr = calculate_trajectory(G,anchors);
    trajectory_val_curr=return_trajec_val(DATA_mat,trajectory_curr,loc); %%% return marker values in each trajectory point 
    
    trajectory_curr=zeros(1,size(DATA_mat,1)); 
    for I2=1:length(trajectory_curr)   
      distance_I=0;
       for J=1:size(DATA_mat,2)
         distance_I=distance_I+ (DATA_mat(I2,J)-trajectory_val_curr(:,J)).^2; 
       end
      distance_I=sqrt(distance_I); 
      trajectory_temp=find(distance_I==min(distance_I));
      trajectory_curr(I2)=trajectory_temp(1); 
    end
    
    segmented_original_cells=[segmented_original_cells  inds_between_points]; 
    ell_index_original=[ell_index_original ; I*ones(length(inds_between_points),1)]; 
    
    Trajectory_new(inds_between_points) = trajectory_curr + add_to_Trajectory; 
    add_to_Trajectory = add_to_Trajectory+ size(trajectory_val_curr,1)-1;
        
    ell_new=[ell_new(1:end-1,:) ; trajectory_val_curr];  
    size_segments(I)=size(trajectory_val_curr,1); 
end

Trajectory_new=re_calculate_Trajectory(ell_new,size_segments,segmented_original_cells,ell_index_original,mat,markers_segment);

[n,m]=size(ell_new);
GM1=zeros(1,n); %% average gamma1
GM2=zeros(1,n); %% average gamma2

 for J=1:n   
     K=find(Trajectory_new==ind(J));
     
     GM1(J)=length(K(ismember(K,gamma1_thresh)))/length(K);
     if ~isempty(K(ismember(K,gamma1_thresh))) 
       GM2(J)=length(K(ismember(K,gamma2_thresh)))/length(K(ismember(K,gamma1_thresh))); 
     else
       GM2(J)=0;   
     end
 end    
 
 ff=zeros(1,n);
 for II=1:n
    D=find(Trajectory_new==ind(II));
    ff(II)=length(D)/length_original;   
 end
 FF=cumsum(ff);    
 dt=define_dt/24; 
 GM1=GM1/dt;
 GM2=GM2/dt;
 

if choose_segment==1  %%% DN+DP segment ends at T=10
  [t_numeric_cal, t_numeric_initial_current , GM2] = return_t_numeric(GM1,GM2,ff,n); 
  if min(diff(t_numeric_cal))<0 %%% if "return_t_numeric" failed
  [t_numeric_cal, t_numeric_initial_current] = return_t_numeric2(GM1,GM2,ff,n); 
  end     
else  %%% DN+DP segment ends at T=4
  [t_numeric_cal, t_numeric_initial_current] = return_t_numeric_SP(GM1,GM2,ff,n); 
  if min(diff(t_numeric_cal))<0 %%% if "return_t_numeric" failed
  [t_numeric_cal, t_numeric_initial_current] = return_t_numeric2_SP(GM1,GM2,ff,n); 
  end  
end

return 

end





