function [evt,xy,bound,interior,eboundt] = mesh_ref(MMele,MMedge,evt,xy,bound,evtY,xyY,boundY)  
%MESH_REF mesh refinement based on longest-edge bisection (LEB) algorithm
%       BMK Minor fix to prevent a single element refinement causing an
%       error.
%
% [evt,xy,bound,interior,eboundt] = mesh_ref(MMelem,MMedge,evt,xy,bound,evtY,xyY,boundY)
%
% input:
%         MMelem     set of overall marked elements
%         MMedge     set of overall marked edges to be refined
%            evt     element mapping matrix (before refinement)
%             xy     vertex coordinate matrix (before refinement)
%          bound     boundary vertices vector (before refinement)
%           evtY     element mapping matrix for midpoints
%            xyY     vertex coordinate vector for midpoints
%         boundY     boundary midpoints vector
%
% output: 
%            evt     element-mapping matrix
%             xy     vertex coordinate matrix
%          bound     boundary vertices vector
%       interior     interior vertices vector
%        eboundt     element-boundary-mapping matrix
%
% The function performs a mesh refinement based on the Longest Edge Bisection 
% (LEB) as a variant of the Newest Vertex Bisection (NVB) algorithm. The 
% longest edges of the marked elements are marked first (reference edges). 
% Hanging nodes are then avoided by completion steps.
%
% Function(s) called: bisection
%                     
%   TIFISS function: LR; 22 June 2018
% Copyright (c) 2018 A. Bespalov, L. Rocchi                                 
     
% Overall number of edges
  nedg = size(xyY,1);
    
% Global numbers of the new nodes (midpoints) to be inserted in the mesh
  newNodes = size(xy,1) + (1:length(MMedge));
  
% Midpoints' global numbers: 
% The following vector 'markEdge' is sparse nedg-by-1 and it contains 
% either 0 in i-th position (if the i-th edge has not been marked)
% or the global number of the midpoint that would be inserted on that edge
  markEdge = sparse(MMedge,1,newNodes,nedg,1);
  
% -----------------------------------------------------------------
% New vertex coordinate matrix xy
% -----------------------------------------------------------------
% Appending the the midpoints' coordinates to current xy
  xy = [xy; xyY(MMedge,:)];
  
% -----------------------------------------------------------------
% New element-mapping matrix 
% -----------------------------------------------------------------
% Refinement of marked elements/edges 
  [evt] = bisection_bk(MMele,MMedge,markEdge,evt,evtY);
   
% -----------------------------------------------------------------
% New interior/boundary nodes and boundary mapping matrix eboundt
% -----------------------------------------------------------------
   
% New boundary nodes
  bound = [bound; nonzeros( markEdge(boundY) )];
  
% New interior vertices vector
  totalnodes = 1:size(xy,1);
  interior   = totalnodes( ~ismember(totalnodes,bound) )';
  
% New element boundary mapping matrix eboundt:
% This is done by updating (locally inside this function) the detail space Y
% for the new refined mesh. To this end, the function P1GRID_DETAIL_SPACE needs 
% to run just up to step 3 (see inside): one way to do this and save some time 
% is to give an optional input parameter such that, using nargin, we can decide 
% if let run the last step 4 or not...
  [evtY,~,boundY,~] = p1grid_detail_space(xy,evt);  
  [belem,bedge]     = find( ismember(evtY,boundY) );
  eboundt           = [belem, bedge];
  
end  % end function
