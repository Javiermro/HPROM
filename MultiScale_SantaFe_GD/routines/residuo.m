% function res = residuo(Fext,Fint,m_LinCond,dofl,doff,lambda)
function res = residuo(Fext,Fext_inic,Fint,m_LinCond,dofl,doff,lambda)
%******************************************************************************************
%*  COMPUTO DEL RESIDUO DEL SISTEMA DISCRETIZADO DE ECUACIONES                            *
%*                                                                                        *                  
%*  ARGUMENTOS DE ENTRADA:                                                                *                  
%*  Fext      : vector de fzas externa                                                    *
%*  Fint      : vector de fzas internas                                                   *
%*  lambda    : factor de escala para fza externa                                         *
%*                                                                                        *                  
%*  ARGUMENTOS DE SALIDA:                                                                 *                  
%*  res       : vector residuo                                                            *                  
%*                                                                                        *                  
%*  Dr. A.E. Huespe, Dr. P.J.Sanchez                                                      *                  
%*  CIMEC-INTEC-UNL-CONICET                                                               *                  
%******************************************************************************************

%Considerando que el trabajo de las reacciones (fuerzas externas sobre los grados restringidos)
%y las variaciones admisibles de desplazamiento son nulas se puede escribir.
% res = Fint(dofl)-lambda*Fext(dofl)+m_LinCond'*Fint(doff);
res = Fint(dofl)-lambda*Fext(dofl)-Fext_inic(dofl)+m_LinCond'*Fint(doff);

%Expresión completa del residuo (no se conoce Fint(doff))
%res = Fint(dofl)-lambda*Fext(dofl)+m_LinCond'*(Fint(doff)-Fext(doff));

%Si se piensa que las fuerzas externas sobre los grados de libertad restringidos (reacciones) se
%actualiza en cada iteración al calcular las Fint en la iter i, quiere decir que los residuos sobre
%estos grados de libertad siempre son nulos, por lo que se puede escribir.
%res = Fint(dofl)-lambda*Fext(dofl);
