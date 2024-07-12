function [mail2,varargout] = ChangeMesh2(mail1,mot,varargin)
% Change element type within a mesh
%
% DUREISSEIX David    L.M.T. STRUCTURES et SYSTEMES  le 29 / 07 / 2002
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 11 / 12 / 2002
%  Ajout passage TRI6 vers TRI3
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 15 / 05 / 2005
%  Ajout passage QUA4 vers RAC2
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 14 / 11 / 2006
%  Ajout mots cles QUAD et LINE
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 24 / 11 / 2007
%  Ajout mots cles TET5b
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 09 / 08 / 2008
%  Ajout mots cles TET9b
%
% Inputs
%  mail1		Mesh to be changed
%  mot			Desired type of element (POI1,TRI3,JOIN)
%                       or change to be applied (LINE,QUAD)
% Optional inputs
%  xcoor1(nbnot,idim)	Coordinates of nodes
%                       (needed if QUAD,TET5b,TET9b)
% Output
%  mail2		New mesh
% Optional output
%  xcoor2(nbnot_add,idim)	New nodes (needed if QUAD,TET5n,TET9b)
%
% If QUAD option is used, the coordinates have to be given,
% and new nodes are created (beware that there could be duplicated
% nodes)
%

nin1 = nargin - 2;
nou1 = nargout - 1;

clear mail2;
GlobalVar;

switch mot
  case 'POI1',
    nbzone1 = length(mail1);
    maile2 = zeros(0,1);
    for izo = 1:nbzone1
      maile1 = mail1{izo}.MAIL;
      [nbel1,nbnn1] = size(maile1);
      nmax = max(max(maile1));
      nbel2 = 0;
      nbnn2 = 1;
      temp = zeros(1,nmax);
      for el1=1:nbel1
        for no1=1:nbnn1
          temp(maile1(el1,no1)) = 1;
        end
      end
      maile2 = unique([maile2 ; find(temp)']);
    end
    mail2{1} = struct('MAIL',maile2,'TYPE',mot);

  case 'TRI3', 
%   On degenere du TRI6 en TRI3 (a meme nb d'elements)
    nbzone1 = TestMeshType(mail1,'TRI6');
    if (nbzone1)
      clear mail2;
      for zo1 = 1:nbzone1
        topo1 = mail1{zo1}.MAIL;
        topo2 = topo1(:,1:3);
        maile2 = struct('MAIL',topo2,'TYPE',mot);
        mail2{zo1} = maile2;
      end
    else
      mot
      error('ne fonctionne actuellement que sur TRI6')
    end

  case 'TET5b', 
%   On transforme du TET4 en TET5b (ajout du centroide en dernier)
    nbzone1 = TestMeshType(mail1,'TET4');
    if (nbzone1 == 0)
      mot
      error('ne fonctionne actuellement que sur TET4')
    end
    if ((nin1~=1) || (nou1~=1))
      nin1
      nou1
      error('Bad number of input or output arguments for TET5b option')
    end
    xcoor1 = varargin{1};
    [nbnot1,idim] = size(xcoor1);
    clear mail2; xcoor2 = zeros(0,idim);
    for zo1 = 1:nbzone1
      topo1 = mail1{zo1}.MAIL;
      nbel1 = size(topo1,1);
      nbnot1_new = nbnot1+size(xcoor2,1);
      topo2 = [topo1 (nbnot1_new+1:nbnot1_new+nbel1)'];
      xcoor_new = zeros(nbel1,idim);
      for el1 = 1:nbel1
        lno1 = topo1(el1,:);
        xcoor_new(el1,:) = 0.25 * sum(xcoor1(lno1,:),1);
      end
      xcoor2 = [xcoor2 ; xcoor_new]; clear xcoor_new
      mail2{zo1} = struct('MAIL',topo2,'TYPE','TET5b');
    end
    varargout{1} = xcoor2;

  case 'TET9b', 
%   On transforme du TET4 en TET9b
%   (ajout du centroide en 5e position puis des 4 autres sub-centroides)
    nbzone1 = TestMeshType(mail1,'TET4');
    if (nbzone1 == 0)
      mot
      error('ne fonctionne actuellement que sur TET4')
    end
    if ((nin1~=1) || (nou1~=1))
      nin1
      nou1
      error('Bad number of input or output arguments for TET9b option')
    end
    xcoor1 = varargin{1};
    [nbnot1,idim] = size(xcoor1);
    clear mail2; xcoor2 = zeros(0,idim);
    for zo1 = 1:nbzone1
      topo1 = mail1{zo1}.MAIL;
      nbel1 = size(topo1,1);
      nbnot1_new = nbnot1+size(xcoor2,1);
      topo2 = [topo1 reshape([nbnot1_new+1:nbnot1_new+5*nbel1],5,nbel1)'];
      xcoor_new = zeros(5*nbel1,idim);
      for el1 = 1:nbel1
        lno1 = topo1(el1,:);
        xcoore1 = xcoor1(lno1,:);
        x5 = 0.25 * sum(xcoore1,1); % centroid
        xcoore2 = xcoore1; xcoore2(1,:) = x5; x6 = 0.25 * sum(xcoore2,1);
        xcoore2 = xcoore1; xcoore2(2,:) = x5; x7 = 0.25 * sum(xcoore2,1);
        xcoore2 = xcoore1; xcoore2(3,:) = x5; x8 = 0.25 * sum(xcoore2,1);
        xcoore2 = xcoore1; xcoore2(4,:) = x5; x9 = 0.25 * sum(xcoore2,1);
        xcoor_new((el1-1)*5+1:(el1-1)*5+5,:) = [x5;x6;x7;x8;x9];
      end
      xcoor2 = [xcoor2 ; xcoor_new]; clear xcoor_new
      mail2{zo1} = struct('MAIL',topo2,'TYPE','TET9b');
    end
    varargout{1} = xcoor2;

  case 'TET4',
    clear mail2;
%   TET5b -> TET4
    nbzone1 = TestMeshType(mail1,'TET5b');
    if (nbzone1)
      disp('  ChangeMesh2: TET5b to TET4')
      for zo1 = 1:nbzone1
        topo1 = mail1{zo1}.MAIL;
	remap = [5 2 3 4 ; 1 5 3 4 ; 1 2 5 4 ; 1 2 3 5];
	nbel1 = size(topo1,1);
	topo2 = zeros(0,4);
	for el1 = 1:nbel1
	  topo1e = topo1(el1,:);
	  topo2 = [topo2 ; topo1e(remap)];
	end
        maile2 = struct('MAIL',topo2,'TYPE','TET4');
        mail2{zo1} = maile2;
        clear maile2 topo2 topo1;
      end
    end
%   TET9b -> TET4
    nbzone1 = TestMeshType(mail1,'TET9b');
    if (nbzone1)
      disp('  ChangeMesh2: TET9b to TET4')
      for zo1 = 1:nbzone1
        topo1 = mail1{zo1}.MAIL;
	remap = [6 2 3 4 ; 5 6 3 4 ; 5 2 6 4 ; 5 2 3 6 ; ...
	         7 5 3 4 ; 1 7 3 4 ; 1 5 7 4 ; 1 5 3 7 ; ...
	         8 2 5 4 ; 1 8 5 4 ; 1 2 8 4 ; 1 2 5 8 ;
	         9 2 3 5 ; 1 9 3 5 ; 1 2 9 5 ; 1 2 3 9 ];
	nbel1 = size(topo1,1);
	topo2 = zeros(0,4);
	for el1 = 1:nbel1
	  topo1e = topo1(el1,:);
	  topo2 = [topo2 ; topo1e(remap)];
	end
        maile2 = struct('MAIL',topo2,'TYPE','TET4');
        mail2{zo1} = maile2;
        clear maile2 topo2 topo1;
      end
    end
    if isempty(mail2)
      mot
      error('ne fonctionne actuellement que sur TET5b ou TET9b')
    end

  case 'JOIN',
    clear mail2;
%   QUA4 -> RAC2
    nbzone1 = TestMeshType(mail1,'QUA4');
    if (nbzone1)
      disp('  ChangeMesh2: QUA4 to RAC2')
      for zo1 = 1:nbzone1
        topo2 = mail1{zo1}.MAIL;
        maile2 = struct('MAIL',topo2,'TYPE','RAC2');
        mail2{zo1} = maile2;
        clear maile2 topo2;
      end
    end
%   QUA8 -> RAC3
    nbzone1 = TestMeshType(mail1,'QUA8');
    if (nbzone1)
      disp('Change from QUA8 to RAC3')
      for zo1 = 1:nbzone1
        topo1 = mail1{zo1}.MAIL;
        topo2 = topo1(:,[1 2 3 5 6 7]);
        maile2 = struct('MAIL',topo2,'TYPE','RAC3');
        mail2{zo1} = maile2;
        clear maile2 topo2 topo1;
      end
    end
    if isempty(mail2)
      mot
      error('ne fonctionne actuellement que sur QUA4 ou QUA8')
    end

  case 'LINE',
%   On degenere les elements "quadratiques" en elements "lineaires"
%   sans changer le nb d'elements. Des noeuds ne sont plus utilises.
    error('Option LINE not yet implemented')

  case 'QUAD',
%   On transforme les elements "lineaires" en elements "quadratiques"
%   sans changer le nombre d'elements. Des noeuds supplementaires
%   sont crees.
    if ((nin1~=1) || (nou1~=1))
      nin1
      nou1
      error('Bad number of input or output arguments for QUAD option')
    end
    xcoor1 = varargin{1};
    [nbnot1,idim] = size(xcoor1);
    clear mail2; xcoor2 = zeros(0,idim);
    nbzo1 = length(mail1);
    for zo1 = 1:nbzo1
      type1 = mail1{zo1}.TYPE;
      topo1 = mail1{zo1}.MAIL;
      nbel1 = size(topo1,1);
      nbnot1_new = nbnot1+size(xcoor2,1);
      switch type1
        case 'SEG2',
%         SEG2 -> SEG3
%disp('ON RANGE LES ARETES APRES LES SOMMETS ???')
          topo2 = [topo1 (nbnot1_new+1:nbnot1_new+nbel1)'];
          xcoor2 = [xcoor2 
                    0.5*(xcoor1(topo1(:,1),:) + xcoor1(topo1(:,2),:))];
          mail2{zo1} = struct('MAIL',topo2,'TYPE','SEG3');
        case 'TRI3',
%         TRI3 -> TRI6
%disp('ON RANGE LES ARETES APRES LES SOMMETS ???')
          topo2 = [topo1 (nbnot1_new+1:nbnot1_new+nbel1)' ...
                         (nbnot1_new+nbel1+1:nbnot1_new+2*nbel1)' ...
                         (nbnot1_new+2*nbel1+1:nbnot1_new+3*nbel1)'];
          xcoor2 = [xcoor2 
                    0.5*(xcoor1(topo1(:,1),:) + xcoor1(topo1(:,2),:))
                    0.5*(xcoor1(topo1(:,2),:) + xcoor1(topo1(:,3),:))
                    0.5*(xcoor1(topo1(:,3),:) + xcoor1(topo1(:,1),:))];
          mail2{zo1} = struct('MAIL',topo2,'TYPE','TRI6');
        otherwise,
          type1
          error('Type of element not implemented to be transformed into QUAD')
      end
    end
    varargout{1} = xcoor2;

  otherwise
    mot
    error('type pas prevu dans ChangeMesh')
end
