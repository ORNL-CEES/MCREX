function A = FFmatrix_fread(file) 

% =======================================================================================
% LECTURE D'UNE MATRICES CALCULEE PAR FREEFEM++
% =======================================================================================
% Fonction :   Lecture d'un fichier contenant une matrice dans l'un des formats de FF++
%              Ces fichiers sont construits par la syntaxe : fout << mat << end ;   
%   
% Syntaxe :    A = FFmatrix_fread(file)
%
% Arguments :  file ... (In)  Nom du fichier a lire 
%              A ...... (Out) Matrice lue, au forma SPARSE de matlab.
% 
% Rappels de FF++ :
% - Les matrices sont construites par la syntaxe :
%
%     matrix matA = formeA(Vh,Vh,solver=<olveur>)
%
% - Le format utilis\'e depend du choix de <Solveur> :
%
%        --------------------------------------------
%        Solveur    stockage                 Matrice
%        --------------------------------------------
%        LU         profil general           generale
%        Crout      profil symetrique        symetrique indefinie
%        Cholesky   profil symetrique        sdp
%        UMFPACK    morse non symetrique     generale
%        GMRES      morse non symetrique     generale
%        CG         morse general            sdp
%
% - Au niveau du fichier a lire, 
% ---------------------------------------------------------------------------------------
% Auteur     : Richard MICHEL, CNRS-LIPhy 
% Historique : 
% V1.0 du 21/08/12 : creation
% =======================================================================================

% --------------------------------------------------------------------------
% Types de matrices admissibles
% --------------------------------------------------------------------------
  MORSE   =  0 ;
  PROFIL  =  1 ;

  EVOL_MORSE = 100 ;
  INCONNU    = 999 ;


% --------------------------------------------------------------------------
% Controle syntaxique de l'appel
% --------------------------------------------------------------------------
  if ( (all(nargout ~= [0 1])) | (nargin ~= 1) )
      error('The calling syntax must be Th = FFmatrix_fread(MATRIXFILE)')
  end


% --------------------------------------------------------------------------
% Controle d'existence du fichier de la matrice
% --------------------------------------------------------------------------
  [fRep,fName,fExt] = fileparts(file) ;  
% if ( isempty(fRep) ) % A desactiver pour laisser la possibilit\'e de chercher dans le path
%     fRep = pwd ;
% end
  if ( isempty(fExt) )
      fExt  = '.mtx' ;
  end
  FileMat = [fRep, filesep, fName, fExt] ; 
  if (exist(FileMat) ~= 2)
      str = sprintf('The file ''%s'' does not exists.\n',FileMat) ;
      error(str) ;
  end


% --------------------------------------------------------------------------
% Ouverture du fichier
% --------------------------------------------------------------------------
  fid = fopen(FileMat) ;


% --------------------------------------------------------------------------
% Lecture de la premiere ligne et determination du type de stockage
% --------------------------------------------------------------------------
  lig_1_ref_MORSE = '# Sparse Matrix (Morse)' ;

  lig = textscan(fid,'%s',1,'delimiter','\n') ;
  lig = deblank(lig{1}) ;

  if ( strcmp(lig,lig_1_ref_MORSE) == 1 )
     type_mat = MORSE ;
  else
     type_mat = INCONNU ; % Ca serait etonnant qu'on utilise des matrices PROFIL !
  end


% --------------------------------------------------------------------------
% Traitement des formats non support\'es
% --------------------------------------------------------------------------
  if ( type_mat == INCONNU )
     msg = sprintf('ERREUR - Unsupported matrix format or evolution of the Morse format\n') ; 
     msg = [ msg sprintf('   expected line 1 = %s\n',lig_1_ref_MORSE) ] ;
     msg = [ msg sprintf('   found line 1    = %s\n',lig{1}) ] ;
     error(msg) ;
  end


% --------------------------------------------------------------------------
% Controle de non evolution pour le format MORSE
% --------------------------------------------------------------------------
% On controle les lignes 2 et 3, soit la suite et fin de l'entete
  lig_2_ref_MORSE = '# first line: n m (is symmetic) nbcoef' ;
  lig_3_ref_MORSE = '# after for each nonzero coefficient:   i j a_ij where (i,j) \in  {1,...,n}x{1,...,m}' ;

  msg_evol = '' ;

  lig = textscan(fid,'%s',1,'delimiter','\n') ;
  lig = deblank(lig{1}) ;

  if ( strcmp(lig,lig_2_ref_MORSE) ~= 1 )
     type_mat = EVOL_MORSE ;
     msg_evol = [ msg_evol sprintf('   expected line 2 = %s\n',lig_2_ref_MORSE) ] ;
     msg_evol = [ msg_evol sprintf('   found line 2    = %s\n',lig{1}) ] ;
  end

  lig = textscan(fid,'%s',1,'delimiter','\n') ;
  lig = deblank(lig{1}) ;

  if ( strcmp(lig,lig_3_ref_MORSE) ~= 1 )
     type_mat = EVOL_MORSE ;
     msg_evol = [ msg_evol sprintf('   expected line 3 = %s\n',lig_3_ref_MORSE) ] ;
     msg_evol = [ msg_evol sprintf('   found line 3    = %s\n',lig{1}) ] ;
  end


  if ( type_mat == EVOL_MORSE )
     error([sprintf('ERREUR - The Morse format has evolved\n') msg_evol]) ; 
  end


% --------------------------------------------------------------------------
% Lecture d'une matrice MORSE
% --------------------------------------------------------------------------
  if ( type_mat == MORSE )
  % Lecture du dimensionnement et de la symetrie
     lig  = textscan(fid,'%f %f %f %f',1) ;  % Lire les entiers en double
     NbL  = lig{1} ;                         % pour pouvoir utiliser SPARSE
     NbC  = lig{2} ;
     sym  = lig{3} ;
     NbNZ = lig{4} ;

  % Lecture des coefficients 
    lig = textscan(fid,'%f %f %f',NbNZ) ; % Lire les entiers en double pour pouvoir
    ptr_lig = lig{1} ;                    % utiliser SPARSE
    ptr_col = lig{2} ;
    coef_NZ = lig{3} ;
  end


% --------------------------------------------------------------------------
% Passage au format sparse de matlab 
% --------------------------------------------------------------------------
  A = sparse(ptr_lig,ptr_col,coef_NZ,NbL,NbC,NbNZ) ;

  if ( sym == 1 ) 
     A = A + tril(A,-1)' ;
  end

% --------------------------------------------------------------------------
% Fermeture du fichier
% --------------------------------------------------------------------------
  fclose(fid) ;


end
