function [m_B,detJ] = matrixB_tria_t1_epd_ld(m_Coord,e_DatElemSet,e_VG)

%*********************************************************************************
%*  Obtiene la matriz B del triแngulo lineal de 3 nodos para EPD                 *
%*  evaluada en xg.                                                              *
%*                                                                               *
%*  ARGUMENTOS DE ENTRADA:                                                       *
%*  จจจจจจจจจจจจจจจจจจจจจ                                                        *
%*  m_Coord : Coordenadas de los nodos que componen el elemento.                 *
%*            m_Coord = [x1  x2  x3                                              *
%*                       y1  y2  y3]                                             *
%*  numeraci๓n de nodos: 3                                                       *
%*                       |\                                                      *
%*                       | \                                                     *
%*                       |  \                                                    *
%*                       |   \                                                   *
%*                       1----2                                                  *
%*  xg      : Vector de posici๓n del punto de Gauss en estudio.                  *
%*            xg = [E n]                                                         *
%*  npe     : N๚mero de nodos por elemento.                                      *
%*                                                                               *
%*  LISTA DE VARIABLES:                                                          *
%*  จจจจจจจจจจจจจจจจจจ                                                           *
%*  E       : coordenada natural                                                 *
%*  n       : coordenada natural                                                 *
%*  dN_E    : Derivadas de las funciones de forma respecto a "E"                 *
%*            dN_E = [dN1_E  dN2_E  dN3_E]                                       *
%*  dN_n    : Derivadas de las funciones de forma respecto a "n"                 *
%*            dN_n = [dN1_n  dN2_n  dN3_n]                                       *
%*  dN_En   : Derivadas de las funciones de forma respecto a "E" y a "n"         *
%*  m_J     : matriz jacobiana de la transformaci๓n                              *
%*            m_J = [dx_E dy_E ; dx_n dy_n]                                        *
%*  dN_xy   : Derivadas de las funciones de forma respecto a "x" y a "y"         *
%*                                                                               *
%*  ARGUMENTOS DE SALIDA:                                                        *
%*  จจจจจจจจจจจจจจจจจจจจ                                                         *
%*  B       : Matriz B para elemento triangular de 3 nodos en estado plano.      *
%*  detJ    : determinante de la matriz jacobiana.                               *
%*********************************************************************************

ndn = e_VG.ndn;
nDime = e_VG.ndime;
%dofpe = e_VG.dofpe;
ntens = e_VG.ntens;
dofpe = e_DatElemSet.dofpe;

% E = xg(1);
% n = xg(2);

% Derivadas de Ni respecto "E" y "n":
% จจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจ
dN_E = [-1,1,0];
dN_n = [-1,0,1];
dN_En = [dN_E;dN_n];

% Matriz Jacobiana:
% จจจจจจจจจจจจจจจจ
% J11 = m_Coord(1,:)*dN_E';
% J12 = m_Coord(2,:)*dN_E';
% J21 = m_Coord(1,:)*dN_n';
% J22 = m_Coord(2,:)*dN_n';
% J = [J11 J12 ; J21 J22];
m_J = dN_En*m_Coord(1:nDime,:)';
detJ = det(m_J);

% Matriz B:
% จจจจจจจจ
m_B = zeros(ntens,dofpe);
dN_xy = m_J\dN_En;
m_B(1,1:ndn:dofpe) = dN_xy(1,:);
m_B(2,2:ndn:dofpe) = dN_xy(2,:);
m_B(4,1:ndn:dofpe) = dN_xy(2,:);
m_B(5,2:ndn:dofpe) = dN_xy(1,:);