function replD=MDrepl1(missD,MDmark,varargin)

%One-dimensional interpolation. Using INTERP1 function.
%
%    replD=MDrepl1(missD, MDmark)
%
%     replD  =  data with replaced values
%     missD  =  data containing missing data
%     MDmark =  mark used to symbolise the missing values
%
%replD = MDrepl1(missD,MDmark,'method') specifies alternate methods.
%The default is linear interpolation.  Available methods are:
%
%  'nearest' = nearest neighbour interpolation
%  'linear'  = linear interpolation
%  'spline'  = cubic spline interpolation
%  'cubic'   = cubic interpolation
%
%All the interpolation methods require that X be monotonic. X can be
%non-uniformly spaced. For more details look HELP INTERP1. 

%Heikki Junninen 22.09.2000
%Updated 


%Check if missing data mark is entered correctly. If not display error.
if ~any(isnan(missD));
   if MDmark~=missD
      error('Make sure your MDmark is correct!')
   end
      
   Inan=find(missD==MDmark);
   missD(Inan)=NaN;
end

[l m]=size(missD);
ind=(1:l)';

%Extract data without MD 
for j=2:m+1
   Icoml=find(~isnan(missD(:,(j-1))));
   dd=[Icoml missD(Icoml,(j-1))];
   
   %Make interpolation 
  replD(:,j)=interp1(dd(:,1),dd(:,2),ind(:,1),varargin{:});
end  
replD(:,1)=[];

