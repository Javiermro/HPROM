function f_TraccHomog2Dom(varargin)
    
   if nargin==1
      %Argumentos: f_TraccHomog2Dom(e_VG)
      %Se inicializa los archivos de resultados
      e_VG = varargin{1};
      fileCompleto = e_VG.fileCompleto;
      fId = fopen([fileCompleto,'.TraccDifDom'],'w');
      fprintf(fId,'#Archivo de las tracciones homogenizadas en diferentes dominios de A\n');
      fclose(fId);
   else
      %Se guarda los resultados de cada paso
      %Argumentos: f_TraccHomog2Dom(time,sigma_new,PGLOC,m_DetJT,Omega_micro,Omega_micro_loc,m_nMacro,e_VG)
      [time,sigma_new,PGLOC,m_DetJT,Omega_micro,Omega_micro_loc,m_nMacro,e_VG] = varargin{:};
      %No se hace incremental, ya que se quiere verificar equilibrio total.
      %Homogeneización en todo el dominio de la celda unitaria
      sigmaHomog = f_HomogArea(sigma_new,Omega_micro,m_DetJT,e_VG);
      %Homogeneización en el dominio localizado
      sigmaHomogLoc = f_HomogArea(sigma_new,Omega_micro_loc,bsxfun(@times,m_DetJT,PGLOC),e_VG);
      %Tracciones
      m_Tracc = m_nMacro*sigmaHomog;
      m_TraccLoc = m_nMacro*sigmaHomogLoc;   
      %Se guarda en archivos los resultados
      fileCompleto = e_VG.fileCompleto;
      fId = fopen([fileCompleto,'.TraccDifDom'],'a');
      fprintf(fId,'%f %f %f %f %f\n',time,m_Tracc,m_TraccLoc);
      fclose(fId);
   end
   
end