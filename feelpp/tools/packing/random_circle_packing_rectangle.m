function [C,R]=random_circle_packing_rectangle(ab,R_min,R_max,cnst,vis)
% Random circle packing inside a rectangle.
%
% INPUT:
%   - ab    : 1-by-2 vector specify rectangle dimensions; ab=[width height].
%             ab=[500 500] is the default setting. Coordinates of the 
%             lower-left and upper-right corners of the rectangle are 
%             (0,0) and (ab(1),ab(2)), respectively.
%   - R_min : minimum circle radius. R_min=min(ab)/100 is the default setting.
%   - R_max : maximum circle radius. R_max=min(ab)/20 is the default setting.
%   - cnst :  set cnst=true to ensure all circles fit into the rectangle.
%             cnst=false is the default setting, meaning that only circle
%             centroids will be constrained to the boundary and interior of
%             the rectangle.
%   - vis   : set vis=false to suppress visualization. vis=true is the
%             default setting.
%
% OUTPUT:
%   - C     : Q-by-2 array of sphere centroids
%   - r     : Q-by-1 array of sphere radii
%

% Default settings
if nargin<1 || isempty(ab)
    ab=[500 500];
elseif numel(ab)~=2 || ~isnumeric(ab) || ~ismatrix(ab) || sum(ab(:)<1E-6 | isinf(ab(:)))>0
    error('Invalid entry for 1st input argument (ab)')
else
    ab=ab(:)';
end

if nargin<2 || isempty(R_min)
    R_min=min(ab)/100;
elseif numel(R_min)~=1 || ~isnumeric(R_min) || R_min>(min(ab)/4-2E-12) || R_min<min(ab)/1E3
    error('Invalid entry for 2nd input argument (D_min)')
end

if nargin<3 || isempty(R_max)
    R_max=min(ab)/20;
elseif numel(R_max)~=1 || ~isnumeric(R_max) || R_max<R_min || R_max>(min(ab)/4-2E-12)
    error('Invalid entry for 3rd input argument (D_max)')
end

if nargin<4 || isempty(cnst)
    cnst=false;
elseif numel(cnst)~=1 || ~islogical(cnst)
    error('Invalid entry for 4th input argument (cnst)')
end

if nargin<5 || isempty(vis)
    vis=true;
elseif numel(vis)~=1 || ~islogical(vis)
    error('Invalid entry for 5th input argument (vis)')
end

% Grid used to keep track of unoccupied space inside the rectangle
dx=max(min(ab)/2E3,R_min/50);
x=0:dx:ab(1);
y=0:dx:ab(2);
[x,y]=meshgrid(x,y);
G=[x(:) y(:)];
clear x y

% Start by placing circles along the edges if cnst=false
dR=R_max-R_min;
Cc=bsxfun(@times,ab,[0 0; 1 0; 1 1; 0 1]); % corner vertices
if ~cnst

      [Xa,Xb]=deal(Cc,circshift(Cc,[-1 0]));    
      Rc=dR*rand(4,1)+R_min;
      [Rc_a,Rc_b]=deal(Rc,circshift(Rc,[-1 0]));

      [C,R]=deal(cell(4,1));
      for i=1:4       
          [Ci,Ri]=SampleLineSegment(Xa(i,:),Xb(i,:),Rc_a(i),Rc_b(i),R_min,R_max);
          Ci(end,:)=[];
          C{i}=Ci;        
          Ri(end)=[];
          R{i}=Ri;        
      end
      C=cell2mat(C);
      R=cell2mat(R);

      % Update grid 
      for i=1:size(C,1), G=update_grid(C(i,:),R(i),G,R_min); end

else

      % Remove all grid points less than R_min units from the boundary
      G_max=G+R_min+1E-12;
      G_min=G-R_min-1E-12;
      chk_in=bsxfun(@le,G_max,ab) & bsxfun(@ge,G_min,[0 0]);
      chk_in=sum(chk_in,2)==2;
      G(~chk_in,:)=[];
      clear G_max G_min chk_in

      C=[]; R=[];
  end

% Begin visualization
if vis
    hf=figure('color','w');
    axis equal
    set(gca,'XLim',[0 ab(1)]+ab(1)*[-1/20 1/20],'YLim',[0 ab(2)]+ab(2)*[-1/20 1/20],'box','on')
    hold on
    drawnow

      Hg=plot(G(:,1),G(:,2),'.k','MarkerSize',2);

      t=linspace(0,2*pi,1E2)';
      P=[cos(t) sin(t)]; % unit circle

      if ~isempty(C)
          for i=1:size(C,1)
              Pm=bsxfun(@plus,R(i)*P,C(i,:));
              h=fill(Pm(:,1),Pm(:,2),'r');
              set(h,'FaceAlpha',0.25)

          end
          drawnow
      end    

end

f=5:5:100;
fprintf('Progress : ')
fprintf(' %-3u ',f)
fprintf(' (%%complete)\n')
fprintf('            ')
Ng=size(G,1);
f=size(G,1)-round((f/100)*Ng);

% Use rejection sampling to populate interior of the rectangle
flag=true;
n=0; cnt=0; m=0;
while ~isempty(G) && cnt<1E4 && (~vis || (vis && ishandle(hf))) 

      n=n+1;

      % New circle
      if flag && (cnt>500 || size(G,1)<0.95*Ng)
          flag=false;
          Rg=R_max*ones(size(G,1),1);
      end

      i=[];
      if cnt<=500 && flag
          X_new=ab.*rand(1,2);   % centroid
      else
          i=randi(size(G,1));
          X_new=G(i,:)+(dx/2)*(2*rand(1,2)-1); 
          X_new=min(max(X_new,[0 0]),ab);
          if cnt>1E3
              Rg(:)=max(0.95*Rg,R_min);
          end
      end
      if isempty(i)
          R_new=dR*rand(1)+R_min; % radius
      else
          R_new=(Rg(i)-R_min)*rand(1)+R_min;
      end

      % Check if the circle fits inside the rectangle when cnst=true
      if cnst
          X_new_max=X_new+R_new+1E-12;
          X_new_min=X_new-R_new-1E-12;
          chk_in=X_new_max<=ab & X_new_min>=[0 0];
          if sum(chk_in)<2
              cnt=cnt+1;
              continue
          end
      end

      % Does the new circle intersect with any other circles?
      if ~isempty(C)
          d_in=sqrt(sum(bsxfun(@minus,C,X_new).^2,2));
          id=d_in<(R+R_new);
          if sum(id)>0
              cnt=cnt+1;
              if ~isempty(i)
                  Rg(i)=min(0.95*Rg(i),min(0.99*(R_new+dx/2),R_max));
                  Rg(i)=max(Rg(i),R_min);
              end
              continue
          end
      end

      % Accept new circle
      cnt=0;
      m=m+1;
      C=cat(1,C,X_new);
      R=cat(1,R,R_new);
      [G,id]=update_grid(X_new,R_new,G,R_min);
      if ~flag, Rg(id)=[]; end

      % Visualization
      if vis && ishandle(hf)
          Pm=bsxfun(@plus,R_new*P,X_new);
          h=fill(Pm(:,1),Pm(:,2),'r');
          set(h,'FaceAlpha',0.25)
          if mod(m,50)==0            
              delete(Hg)
              Hg=plot(G(:,1),G(:,2),'.k','MarkerSize',2);
              drawnow
          end
      end

      % Progress update   
      if size(G,1)<=f(1)
          f(1)=[];
          fprintf('*    ')
      end

end
fprintf('\n')

% Show rectangle
if vis && ishandle(hf)
    Cc=[Cc;Cc(1,:)];
    plot(Cc(:,1),Cc(:,2),'--k','LineWidth',2)    
    delete(Hg)
end

if nargout<1, clear C R; end

function [G,id]=update_grid(X_new,R_new,G,R_min)
% Remove grid points within R_new+R_min units of new circle

D=sum(bsxfun(@minus,G,X_new).^2,2);
id=D<(R_new+R_min+1E-12)^2;
G(id,:)=[];

function [C,R]=SampleLineSegment(Xa,Xb,Ra,Rb,R_min,R_max)
% Place circles along line segment between points Xa and Xb

r=Xb-Xa;
L=norm(r);
r=r/norm(L);

dR=R_max-R_min;
C=Xa; R=Ra;
while true    
    R_new=dR*rand(1)+R_min;
    C_new=C(end,:)+r*(R(end)+R_new+R_max*rand(1));        
    D=L - norm(C_new + r*(R_new+Rb) - Xa); % will there be enough space left for the end point with radius Rb?     
    if D<2*(R_min+1E-12)    
        C=cat(1,C,Xb);
        R=cat(1,R,Rb);
        break
    else
        C=cat(1,C,C_new);
        R=cat(1,R,R_new);
    end     
end