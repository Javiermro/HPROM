function [vectVHElem_new]  = ...
    Reg_param_MICRO_SF(xx_MICRO,e_VarEst_new,...
    e_VarEst_old,e_VarAux,delta_uMICRO,e_VG,e_DatSet,vectVHElem_old,...
    kinf,omegaMicro,longFis,omegaMicroL)

% Calculo de la distancia S_mu en la microcelda
% *********************************************
inc_dissip_CELL       = 0;
m_valHomog_old_5      = 0;
m_valHomog_new_5      = 0;
m_valHomog_old        = 0;
m_valHomog_new        = 0;

damage_micro          = 0;

nSet = e_VG.nSet;
vectVHElem_new = vectVHElem_old;
for iSet = 1:nSet
    eltype            = e_DatSet(iSet).e_DatElem.eltype;
    switch eltype
        case 32
            nElem_MICRO       = e_DatSet(iSet).nElem;
            e_DatElemSet      = e_DatSet(iSet).e_DatElem;
            wg                = e_DatSet(iSet).e_DatElem.wg;
            npg               = e_DatSet(iSet).e_DatElem.npg;
            m_DetJT           = e_DatSet(iSet).m_DetJT;
            hvar_old          = e_VarEst_old(iSet).hvar;
            hvar_new          = e_VarEst_new(iSet).hvar;
%            area_elem         = e_DatSet(iSet).m_VolElem;
            p_IntDissip       =  e_DatElemSet.pointersVHE.i_IntDissip  ;
            m_VarHistElemNew  =  e_VarEst_new(iSet).VarHistElem ;
 %           le_Elem = e_DatSet(iSet).e_DatElem.le_Elem;
            for iElem_MICRO = 1:nElem_MICRO
                hvar_OLD        = reshape(hvar_old(:,iElem_MICRO),[] , npg);
                hvar_NEW        = reshape(hvar_new(:,iElem_MICRO),[] , npg);
                m_pesoPG        = m_DetJT(:,iElem_MICRO).*wg;
                m_valHomog_old_5  = m_valHomog_old_5+ hvar_OLD(5,:)*m_pesoPG;
                m_valHomog_new_5  = m_valHomog_new_5+ hvar_NEW(5,:)*m_pesoPG;
                m_valHomog_old  = m_valHomog_old+ hvar_OLD(7,:)*m_pesoPG;
                m_valHomog_new  = m_valHomog_new+ hvar_NEW(7,:)*m_pesoPG;
                damage_micro = damage_micro +  hvar_NEW(1,:)*m_pesoPG;
                % material que dania y ablanda con deformacion - ielem_material==1 la segunda condicion revisa que en los 4 PGS se este presentando
                % carga, de lo contrario no tiene en cuenta el elemento finito para el calculo de la longitud de fallo
                m_IntDissip     =  m_VarHistElemNew (p_IntDissip   ,iElem_MICRO)   ;
                inc_dissip_CELL = inc_dissip_CELL + m_IntDissip;
            end
    end
end

m_valHomog_old_5  = m_valHomog_old_5 /omegaMicro ;
m_valHomog_new_5  = m_valHomog_new_5 /omegaMicro;
m_valHomog_old    = m_valHomog_old /omegaMicro ;
m_valHomog_new    = m_valHomog_new /omegaMicro;
vectVHElem_new(5) = m_valHomog_new_5-m_valHomog_old_5; % CRITERIO PARA RECUPERAR EL ELEMENTO FINITO ESTANDAR EN MODELOS MULTIESCALA.
vectVHElem_new(4) = vectVHElem_old(4)+ kinf*(m_valHomog_new-m_valHomog_old); % VARIABLE DE SUAVIZADO MULTIESCALA

if longFis==0; longFis=1e-16; end
vectVHElem_new(7) = longFis ;     % S_mu
vectVHElem_new(6) = omegaMicro/longFis; %long caracteristica del modelo multiescala

vectVHElem_new(8) = damage_micro/omegaMicroL;

if (longFis==1e-16);
    FACT=0;
else
    FACT=(kinf/omegaMicro);
end
% Variable de disipacion acumulada en la microcelda
vectVHElem_new(1)  = vectVHElem_old(1) + FACT*inc_dissip_CELL;



