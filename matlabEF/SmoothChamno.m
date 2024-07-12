function [chamno2] = SmoothChamno(chamno1,modl1,mail1,intgno1,xcrit1);
% Smooth element field defined at nodes with discontinuity criteria
%
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 11 / 04 / 2004
%
% Inputs
%   chamno1	Field to be smoothed
%   modl1	Its model
%   mail1	Its mesh
%   intgno1	Its location points (each node of each element)
%   xcrit1	Relative proximity criteria
% Outputs
%   chamno2	Smoothed field
% Comments
%   See also ChamnoToChpo

% Assemble chamno1 into U1a
% """""""""""""""""""""""""
[listComp1,listUnit1] = ListCompCham2(chamno1);
Ncomp1 = length(listComp1);
Nptg1 = Nptg(mail1,intgno1);
numerptg1 = [1:Nptg1];
mapComp1 = MapCompCham(chamno1,mail1,intgno1,numerptg1,listComp1);
U1 = ChamToVect3(chamno1,mail1,intgno1,numerptg1,mapComp1,listComp1);
U1a = U1(mapComp1); % U1a(ptg,component)

UU = U1a;

% Local numbering of nodes
% """"""""""""""""""""""""
nmail1 = ChangeMesh2(mail1,'POI1');
numer1 = sort(nmail1{1}.MAIL');
clear nmail1;
numer_inv1 = InverseList(numer1,max(numer1));
lmail1 = RenumMesh(mail1,numer_inv1); % Locally numbered mesh

% Build inverse graph (elements connected to nodes)
% """""""""""""""""""
% lnod1(lind1(i):lind1(i+1)-1) : nodes connected to element i
[lnod1,lind1] = MeshToGraph(lmail1);
% lnod2(lind2(i):lind2(i+1)-1) : elements connected to node i
% lord2(lind2(i):lind2(i+1)-1) : positions of node i in elements
[lnod2,lind2,lord2] = InverseGraph2(lnod1,lind1);
clear numer_inv1 lmail1;

% Correspondance with the ptg number in each element
% lnod3(lind2(i):lind2(i+1)-1) : ptg global numbers of node i
lnod3 = zeros(1,length(lnod2));
% Lnptg1(izo1) = number of ptg in zone izo1
Lnptg1 = Nptg(mail1,intgno1,[1:length(mail1)]);
Lnptg2 = cumsum(Lnptg1);
% Lnbel1(izo1) = number of elements in zone izo1
Lnbel1 = Nelements2(mail1,[1:length(mail1)]);
Lnbel2 = cumsum(Lnbel1);
% Loop on nodes
Nbnot1 = length(lind2) - 1;
for ino1 = 1:Nbnot1
  lelem1 = lnod2(lind2(ino1):lind2(ino1+1)-1); % connected elements
  lptg1  = lord2(lind2(ino1):lind2(ino1+1)-1); % position in elements
  for i = 1:length(lelem1)
%   node at local position ipos1 in element of global number ielem1
    ielem1 = lelem1(i);
    ipos1  = lptg1(i);
%   zone izo1 and local element number iel1
    junk = find(Lnbel2>=ielem1);
    izo1 = junk(1);
    if izo1 == 1; iel1 = ielem1;
    else;         iel1 = ielem1 - Lnbel2(izo1-1);
    end
%   corresponding global ptg number iptg1
    if izo1 == 1; iptg1 = (iel1-1)*(Lnptg1(izo1)/Lnbel2(izo1))+ipos1;
    else;         iptg1 = Lnptg2(izo1-1)+(iel1-1)*(Lnptg1(izo1)/Lnbel2(izo1))+ipos1;
    end
%
    lnod3(lind2(ino1)+i-1:lind2(ino1)+i-1) = iptg1;
    clear iptg1;
  end
  clear lelem1 lptg1;
end

% For node ino1,
% the connected elements are lnod2(lind2(ino1):lind2(ino1+1)-1)
% the corresponding local positions in elements are
%   lord2(lind2(ino1):lind2(ino1+1)-1)
% the corresponding global ptg numbers are
%   lnod3(lind2(ino1):lind2(ino1+1)-1)
% the corresponding values are 
%   U1a(lnod3(lind2(ino1):lind2(ino1+1)-1),:)

% Smooth values and backstore them in U1a
% """""""""""""""""""""""""""""""""""""""
Nbnot1 = length(lind2) - 1;
for ino1 = 1:Nbnot1
  xval1 = U1a(lnod3(lind2(ino1):lind2(ino1+1)-1),:);
  nval1 = size(xval1,1);
%
% Group the values: in a same group, they must not differ
% from more than xcrit1 from each other
  xcross = zeros(nval1,nval1);
  for ival1 = 1:nval1
    for jval1 = ival1+1:nval1
      xcross(ival1,jval1) = norm(xval1(ival1,:) - xval1(jval1,:)) / ...
                            norm(xval1(ival1,:) + xval1(jval1,:));
      xcross(jval1,ival1) = xcross(ival1,jval1);
    end
  end
% Exclusion-from-the-group procedure
  color = zeros(1,nval1);
  ncolor = 0; lgroup = [1:nval1]; j = lgroup;
  while ~isempty(j)
    ncolor = ncolor + 1; lgroup = lgroup(j); color(lgroup) = ncolor;
    [junk,j] = find(xcross(lgroup(1),lgroup) > xcrit1);
  end
  if ncolor > 1
    disp(['Node ' int2str(numer1(ino1)) ' smoothing ' int2str(ncolor) ...
          ' groups from ' int2str(nval1) ' elements'])
%xcross
%ii = lnod2(lind2(ino1):lind2(ino1+1)-1)
%for icolor = 1:ncolor
%ii(find(color==icolor))
%end
  end
  for icolor = 1:ncolor
    lgroup = find(color==icolor);
    ngroup = length(lgroup);
    xaverage1 = sum(xval1(lgroup,:),1)/ngroup;
    for igroup = 1:ngroup
      xval1(lgroup(igroup),:) = xaverage1;
    end
  end
%
% Storage
  U1a(lnod3(lind2(ino1):lind2(ino1+1)-1),:) = xval1;
  clear xval1 xcross;
end

% Disassemble U1a into chamno2
% """"""""""""""""""""""""""""
U1(mapComp1) = U1a;
chamno2 = VectToCham2(U1,mail1,intgno1,numerptg1,mapComp1,listComp1);
clear lnod1 lind1 lnod2 lind2 lord2 lnod3;
clear U1a U1 numerptg1 mapComp1;
