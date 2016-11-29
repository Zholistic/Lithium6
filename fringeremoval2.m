function [ odimages,optrefimages,avgimage,timer ] = fringeremoval2( absimages,refimages,bgmask,returnimgs )
% FRINGEREMOVAL2 - Fringe noise removal from absorption images using 'mix-and-match' method 
% Creates an optimal reference image for each absorption image in a set as a linear
% combination of reference images with coefficients chosen to minimize least-squares residuals. 
% The coefficients are obtained by solving a linear set of equations using matrix inverse by LU decomposition. 
%
% Algorithm inspired by J. Kronj\"ager PhD thesis (Hamburg-2007).
%
% Syntax:  [odimages,optrefimages,avgimage] = fringeremoval(absimages,refimages,bgmask);
%
% Inputs:
%    absimages - cell array or three-dimensional matrix of raw absorption image data
%    refimages - cell array or three-dimensional matrix of reference images for background division
%    bgmask - Binary image specifying background region used in removal, 1=bg,0=data 
%    returnimgs (optional) - vector of indices referring to the absorption images to process and return. All reference images will be used as a basis for fringe-removal.
%
% Outputs:
%    odimages - Cell array or three-dimensional matrix (depending on input) of optimised optical density images: -log(abs-dark)/(ref-dark)
%    optrefimages - Cell array or three-dimensional matrix of optimal reference images
%    avgimage - Two-dimensional containing average optical density image computed over full set

% Example: 
%    
% Other m-files required: none
%
%
% Author: Shannon Whitlock
% email: S.M.Whitlock@uva.nl
% May 2009; Last revision: 05-May-2009

%------------- BEGIN CODE --------------
timer=[];

tic
if not(exist('bgmask','var'))
   bgmask=ones(size(absimages(:,:,1)));end     %check if background region specified
k=find(bgmask(:)==1);                      %index k specifying background region

cellmode=iscell(absimages);            % check the input data format (cell array of image matrices or three-dimensional array)

if cellmode                            % check if input is a cell array
    nimgs=length(absimages);           % number of images contained in cell array
    nimgsR=length(refimages);
    
    xdim=length(absimages{1}(1,:));
    ydim=length(absimages{1}(:,1));    % image dimensions used to reshape data for output

  
    R=single(reshape(cat(3,refimages{:}),xdim*ydim,nimgsR));
    A=single(reshape(cat(3,absimages{:}),xdim*ydim,nimgs)); %flatten images
else
     nimgs=size(absimages,3);           % number of images contained in cell array
     nimgsR=size(refimages,3);           % number of images contained in cell array
     
     xdim=size(absimages,2);
     ydim=size(absimages,1);            % image dimensions used to reshape data for output

     R=single(reshape(refimages,xdim*ydim,nimgsR));
     A=single(reshape(absimages,xdim*ydim,nimgs)); %flatten images 
end

if not(exist('returnimgs','var'))
    returnimgs=1:nimgs;
end

% ensure there are no duplicate reference images (R must be linearly independent)
%R=unique(R','rows')'; % comment this line if you run out of memory
[L,U,p] = lu(R(k,:)'*R(k,:),'vector'); % LU decomposition
%L = chol (R(k,:)'*R(k,:), 'lower') ;  % Cholesky decomposition; may perform
                                       % better in some circumstances
avgimage=zeros(size(bgmask));
odimages=cell(nimgs,1);
optrefimages=cell(nimgs,1);

timer(1)=toc;
for j=returnimgs
    tic
    b=R(k,:)'*A(k,j);     %b=sum(R(k,:).*repmat(A(k,j),[1,length(R(1,:))]))';

    % Obtain coefficients c which minimise least-square residuals  
    lower.LT = true; upper.UT = true;
    c = linsolve(U,linsolve(L,b(p,:),lower),upper); % using LU decomposition
    %c = linsolve(L,linsolve(L,b,lower),upper); % using Cholesky decomposition

    if cellmode 
        % compute optimised reference image
        O=R*c;               % flat optimised reference image
     
        odimages{j}=reshape(-log(A(:,j)./O),[ydim,xdim]);     % optical density image
        optrefimages{j}=reshape(O,[ydim,xdim]);              % optimised reference image
        avgimage=avgimage+odimages{j};                        % add this odimage to averaged image
    else
      
        O(:,j)=R*c; %O(:,j)=(sum(repmat(c',[xdim*ydim,1]).*R,2));
    end
    timer=[timer,toc];
end

if cellmode
    avgimage=avgimage/nimgs;
else
    odimages=reshape(-log(A(:,returnimgs)./O),[ydim,xdim,length(returnimgs)]);
    optrefimages=reshape(O,[ydim,xdim,length(returnimgs)]);
    avgimage=mean(odimages,3);
end
 

%------------- END CODE --------------