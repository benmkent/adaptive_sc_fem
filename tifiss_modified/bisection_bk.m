function [evt] = bisection(MMelem,MMedge,markEdge,evt,evtY)
%BISECTION performs bisections of marked elements/edges 
%
% [evt] = bisection(MMelem,MMedge,markEdge,evt,evtY)
%
% input:
%         MMelem      vector of (overall) marked elements
%         MMedge      vector of (overall) marked edges to be refined
%       markEdge      marked edge-midpoint-position vector
%            evt      element mapping matrix
%           evtY      element mapping matrix for midpoints
%
% output: 
%            evt      element mapping matrix (after refinement)
%
% See also MESH_REF
%                     
%   TIFISS function: LR; 22 June 2018
% Copyright (c) 2018 A. Bespalov, L. Rocchi  

% Extract marked elements according to the refinement type 1/2/3
  refxelem = sum( ismember( evtY(MMelem,:) , MMedge) , 2);
  Mone     = sort( MMelem( refxelem == 1 ) );                                     % Elements marked for 1 bisection
  Mtwo     = MMelem( refxelem == 2 );                                             % Elements marked for 2 bisections
  Mthree   = sort( MMelem( refxelem == 3 ) );                                     % Elements marked for 3 bisections
        
% -----------------------------------------------------------------------------
% STEP 1: refine elements marked for 1 bisection (green refinement)
% -----------------------------------------------------------------------------

  vtx      = evt(Mone,:);                                                         % Nodes of elements in Mone
  newnodes = full( markEdge( evtY(Mone,2) ) );                                    % Midpoints of longest edges of Mone

% Refine elements
  firstChild  = [Mone, vtx(:,2), newnodes, vtx(:,1)];                             % Top child elements
  secondChild = [Mone, vtx(:,3), newnodes, vtx(:,2)];                             % Bottom child elements

% Save elements
  newrows_one = [firstChild, secondChild];       
  newrows_one = reshape(newrows_one',4,2*length(Mone))';

% Transform the new rows per each child in a cell array: it will be length(Mone) 
% long and each cell is a 2-by-4 matrix (where 2 is due to the 2 children)
  cell_one = mat2cell(newrows_one, 2*ones(1,size(newrows_one,1)/2), 4);
  
  
% -----------------------------------------------------------------------------
% STEP 2: refine elements marked for 2 bisections (blue refinement)
% -----------------------------------------------------------------------------

% Separate the elements marked for blue-left from those marked for blue-right 
% refinement. NOTE that: 
% - blue-left refinement means that (also) the 1st element's edge is bisected
% - blue-right refinement means that (also) the 3rd element's edge is bisected
  checkleft    = repmat( [1 1 0] , [length(Mtwo),1] );
  checkright   = repmat( [0 1 1] , [length(Mtwo),1] );
  markEdgesxel = ismember( evtY(Mtwo,:) , MMedge);
  Mtwoleft     = Mtwo( sum( abs(markEdgesxel - checkleft)  , 2) == 0 );           % Marked elements for blue-left refinement
  Mtworight    = Mtwo( sum( abs(markEdgesxel - checkright) , 2) == 0 );           % Marked elements for blue-right refinement
  Mtwoleft     = sort( Mtwoleft  );
  Mtworight    = sort( Mtworight );
  
% Allocate an empty vector for all children space
  [firstChildLeft,  secondChildLeft,  thirdChildLeft ,...
   firstChildRight, secondChildRight, thirdChildRight] = deal([]); 
  
  if ~isempty(Mtwoleft) 
      vtxLeft     = evt(Mtwoleft,:);                                              % Nodes of elements in Mtwoleft 
      leNodeLeft  = full( markEdge( evtY(Mtwoleft,2) ) );                         % Midpoints of longest edges of Mtwoleft's elements
      mrkNodeLeft = full( markEdge( evtY(Mtwoleft,1) ) );                         % Midpoints of 1st edges of Mtwoleft's elements
      % Left-Blue refinement
      firstChildLeft  = [Mtwoleft, vtxLeft(:,2),  leNodeLeft, vtxLeft(:,1)];      % Top child elements
      secondChildLeft = [Mtwoleft, vtxLeft(:,2), mrkNodeLeft,   leNodeLeft];      % Middle child elements
      thirdChildLeft  = [Mtwoleft,   leNodeLeft, mrkNodeLeft, vtxLeft(:,3)];      % Bottom child elements
  end
  
  if ~isempty(Mtworight)
      vtxRight     = evt(Mtworight,:);                                            % Nodes of elements in Mtworight
      leNodeRight  = full( markEdge( evtY(Mtworight,2) ) );                       % Midpoints of longest edges of Mtworigth's elements
      mrkNodeRight = full( markEdge( evtY(Mtworight,3) ) );                       % Midpoints of 3rd edges of Mtworight's elements
      % Right-Blue refinement
      firstChildRight  = [Mtworight, vtxRight(:,1), mrkNodeRight,   leNodeRight]; % Top child elements
      secondChildRight = [Mtworight,   leNodeRight, mrkNodeRight, vtxRight(:,2)]; % Middle child elements
      thirdChildRight  = [Mtworight, vtxRight(:,3),  leNodeRight, vtxRight(:,2)]; % Bottom child elements
  end   
  
% Save elements  
  newrows_two_left  = [firstChildLeft,  secondChildLeft,  thirdChildLeft];
  newrows_two_left  = reshape(newrows_two_left',4,3*length(Mtwoleft))';
%
  newrows_two_right = [firstChildRight, secondChildRight, thirdChildRight];
  newrows_two_right = reshape(newrows_two_right',4,3*length(Mtworight))';
        
% Transform the new rows per each child in cells array: they will be length(Mtwoleft) 
% and length(Mtworight) long and each cell would be a 3-by-4 matrix 
% (where 3 is due to the 3 children)
  cell_two_left  = mat2cell(newrows_two_left,  3*ones(1,size(newrows_two_left,1)/3),  4);
  cell_two_right = mat2cell(newrows_two_right, 3*ones(1,size(newrows_two_right,1)/3), 4);
  

% -----------------------------------------------------------------------------             
% STEP 3: refine elements marked for 3 bisections (bisec3 refinement)
% -----------------------------------------------------------------------------

  vtx           = evt(Mthree,:);                                                  % Nodes of elements in Mthree
  leNode        = full( markEdge( evtY(Mthree,2) ) );                             % Midpoints of longest edges of Mthree
  mrkNodeLeft   = full( markEdge( evtY(Mthree,1) ) );                             % Midpoints of 1st edges of Mthree
  mrkNodeRight  = full( markEdge( evtY(Mthree,3) ) );                             % Midpoints of 3rd edges of Mthree
  
% Refine elements
  firstChild    = [Mthree, vtx(:,1), mrkNodeRight,   leNode];                     % Top child elements
  secondChild   = [Mthree,   leNode, mrkNodeRight, vtx(:,2)];                     % First middle child elements
  thirdChild    = [Mthree, vtx(:,2),  mrkNodeLeft,   leNode];                     % Second middle child elements
  fourtChild    = [Mthree,   leNode,  mrkNodeLeft, vtx(:,3)];                     % Bottom child elements

% Save elements
  newrows_three = [firstChild, secondChild, thirdChild, fourtChild];
  newrows_three = reshape(newrows_three',4,4*length(Mthree))';
  
% Transform the new rows per each child in a cell array: it will be length(Mthree) 
% long and each cell would be a 4-by-4 matrix (where 4 is due to the 4 children)
  cell_three = mat2cell(newrows_three, 4*ones(1,size(newrows_three,1)/4), 4);
  

% -----------------------------------------------------------------------------
% STEP 4: create the new element mapping matrix
% -----------------------------------------------------------------------------    

% New elements are inserted in the current evt in the following way.
% Suppose the 23rd element was marked for 2 refinements. This will produce 3 
% children. Then, the 23rd row in evt should be deleted and 3 new rows 
% corresponding to the 3 children will take its place. The 3 children will 
% stay on the 23rd, 24th, and 25th row. 
% Accordingly, the old 24th row will become the 26th row.
%
% The insertion explained above is done by transforming the evt matrix in a cell array

% Create a cell containing the rows of evt (element)
  evt = [(1:size(evt,1))', evt];
  evtcell = mat2cell(evt, ones(1,size(evt,1)), 4);

% Update elements with the new rows of all children that have to be inserted
  evtcell( Mone )      = cell_one(:);
  evtcell( Mtwoleft )  = cell_two_left(:);
  evtcell( Mtworight ) = cell_two_right(:);
  evtcell( Mthree )    = cell_three(:);
 
% Convert the cell to matrix
  evt = cell2mat( evtcell );  
  evt = evt(:,[2 3 4]);
  
end % end function