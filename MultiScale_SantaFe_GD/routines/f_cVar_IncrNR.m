% function c_GdlCond = f_cVar_IncrNR(du_iter,c_GdlCond,e_DatSet,e_VG)
function c_GdlCond = f_cVar_IncrNR(du_iter,c_GdlCond,e_DatSet,e_VarEst_new,e_VarAux,e_VG)

   %Esta función calcula los incrementos de las variables condensadas, asumiendo que son variables
   %internas del elemento.
   nSet = e_VG.nSet;
   ndime = e_VG.ndime;
   ndn = e_VG.ndn;
   ntens = e_VG.ntens;
   c_dGdlCond = cell(nSet,1);
   for iSet = 1:nSet
      nElem = e_DatSet(iSet).nElem;
      conec = e_DatSet(iSet).conec;
      eltype = e_DatSet(iSet).e_DatElem.eltype;
      % Grados de libertad y coordenadas de los nodos de los elementos del set
      m_dofElemSet = f_DofElem(reshape(conec',1,[]),ndn);
      % Incrementos de desplazamientos de los nodos de cada elemento del set.
      m_duElemSet = reshape(du_iter(m_dofElemSet),[],nElem);
      switch eltype
         case 10
            m_beta    = c_GdlCond{iSet,1};
            m_resbT   = c_GdlCond{iSet,2};
            m_kbuT    = c_GdlCond{iSet,3};
            m_invkbbT = c_GdlCond{iSet,4}; 
            m_dGdlCond = zeros(ndime,nElem);
            for iElem = 1:nElem
            %parfor iElem = 1:nElem
               m_dGdlCond(:,iElem) = -m_invkbbT(:,:,iElem)*(m_resbT(:,iElem)+...
                  m_kbuT(:,:,iElem)*m_duElemSet(:,iElem));
             
               m_beta(:,iElem)  = m_beta(:,iElem)+ m_dGdlCond(:,iElem) ;
              
            end
            c_dGdlCond{iSet}  = m_dGdlCond;
            c_GdlCond{iSet,1} = m_beta ;

         case {21,22,23}
             
            %Recuperación de variable de salto (variable interna condesada)
            beta_SDA     = c_GdlCond{iSet,1};
            Dbeta_SDA    = c_GdlCond{iSet,2};
            beta_MSL     = c_GdlCond{iSet,3};
            Dbeta_MSL    = c_GdlCond{iSet,4};
            % SDA ELEMENT
            m_kbu_SDA      = c_GdlCond{iSet,5}  ;
            m_kbbIN_SDA    = c_GdlCond{iSet,6}  ;
            m_Res_beta_SDA = c_GdlCond{iSet,7}  ;
            % MSL ELEMENT
            m_kbu_MSL      = c_GdlCond{iSet,8} ;
            m_kbbIN_MSL    = c_GdlCond{iSet,9}  ;
            m_Res_beta_MSL = c_GdlCond{iSet,10} ;

            dbeta_iter_SDA = zeros(ndime,nElem);
            dbeta_iter_MSL = zeros(ntens,nElem);
            for iElem = 1:nElem

                % INCREMENTO EN SALTO - SDA
                dbeta_iter_SDA(:,iElem) = -m_kbbIN_SDA(:,:,iElem)*...
                    (m_Res_beta_SDA(:,iElem) + m_kbu_SDA(:,:,iElem)*m_duElemSet(:,iElem));   %*du_iter_gen(dof));
                
                beta_SDA(:,iElem)     =  beta_SDA(:,iElem) + dbeta_iter_SDA(:,iElem);
                Dbeta_SDA(:,iElem)    = Dbeta_SDA(:,iElem) + dbeta_iter_SDA(:,iElem);
    
               % INCREMENTO EN SALTO - MSL
               dbeta_iter_MSL(:,iElem) = -m_kbbIN_MSL(:,:,iElem)*...
                   (m_Res_beta_MSL(:,iElem) + m_kbu_MSL(:,:,iElem)*m_duElemSet(:,iElem));   %*du_iter_gen(dof));

               beta_MSL(:,iElem)     =  beta_MSL(:,iElem) + dbeta_iter_MSL(:,iElem);
               Dbeta_MSL(:,iElem)    = Dbeta_MSL(:,iElem) + dbeta_iter_MSL(:,iElem);
    
%             %parfor iElem = 1:nElem
%                m_dGdlCond(:,iElem) = -m_invkbbT(:,:,iElem)*(m_resbT(:,iElem)+...
%                   m_kbuT(:,:,iElem)*m_duElemSet(:,iElem));


            end
 %           c_dGdlCond{iSet}  = m_dGdlCond    ;
            c_GdlCond{iSet,1} = beta_SDA      ;
            c_GdlCond{iSet,2} = Dbeta_SDA     ;
            c_GdlCond{iSet,3} = beta_MSL      ;
            c_GdlCond{iSet,4} = Dbeta_MSL     ;

      end
   end

end