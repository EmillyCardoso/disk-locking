pro gallet2013_legendadomassas

  common CNT, Msun, Rsun, P0, myr, tau_d, periodo_segundos,velo_solar, Gconst, ano_segundo, Lsun

  constants

;close, 1, file

  Nmodels = 3 ;numero de modelos/massas a ser lido

  dk = 1

  imass   = [0.5,0.8,1.0]  ; array com as massas
  titlemass = ['0.5 $M_\odot$','0.8 $M_\odot$','1.0 $M_\odot$']
    color   = ['g--','r--','b--']

  ; Colocar de 0 ao último ítem da lista, ou específicar como 0,0; 1,1; 1,2 e etc.
  for jj = 0, 2 do begin
    iind  = ROUND(10.*imass[jj])   ; inteiro correspondente a 10*massa (5, 8 e 10, no caso)

    ;   Atribuindo um arquivo de dados à "file":


    if (iind lt 10) then file = 'C:\Users\Emilly\OneDrive\IDL\Dados\dados_0'+strtrim(iind,1)+'_massa.txt'
    if (iind ge 10) then file = 'C:\Users\Emilly\OneDrive\IDL\Dados\dados_'+strtrim(iind,1)+'_massa.txt'

    ;   Obtendo o número de linhas do arquivo:
    ;
    nrows = File_lines(file)
    ;
    ;   declarando uma variáel string a ser lida
    ;
    line = ""
    ;
    ;   abrindo o arquivo:
    ;
    openr, 1, file
    ;
    ;   Lendo a 1a e 2a linhas (tem um problema com a primeira linha, pq log Tc e log R0c
    ;     ora aparece junto, como no arquivo de 1 massa solar, ora aparece separado, como
    ;     nos arquivos de 0.5 e 0.8 massas solares. Isto faz com que o IDL lesse (correta-
    ;     mente a primeira linha, mas atribuisse mais colunas do que realmente tem. Para
    ;     evitar este problema, eu estou lendo as duas primeiras linhas, e pegando o se-
    ;     gundo valor de 'line' - que é 13, como deve ser. Mais facil seria simplesmente
    ;     deletar esta primeira linha dos arquivos, mas agora já foi.)
    ;     
    ;   Primeira linha:

    readf, 1, line
    ;
    ;   Segunda linha:
    ;
    readf,1,line
    ;
    ;   obtendo o número de colunas do arquivo
    ;
    cols = N_Elements(StrSplit(line, /RegEx, /Extract))
    ;
    ;   declarando o tamanho do array (e excluíndo a primeira linha)
    ;
    data = fltarr(cols,nrows-1)
    ;
    ;   Imprimindo algum diagnóstico:
    ;
    print,'O arquivo contem ',cols,' colunas e ',nrows,' linhas'

    print, nrows-1
    ;
    ;    help,data
    ;
    ;    stop

    ;
    ;   Retroagindo à primeira linha (antigo "rewind" das fitas dat...)
    ;
    Point_Lun, 1, 0
    ;
    ;;   Lendo o array
    ;;
    ;    readf,1,line   ; apenas o cabeçalho
    ;    readf,1,data   ; e agora, o resto do arquivo
    ;    free_lun,1     ; liberando a unidade lógica


    NDIM = nrows-1
    NDIM2 = 10000L


;  NDIM = 424
;  NDIM2 = 10000L

  data = dblarr(13,NDIM)
  produto = dblarr(NDIM)
  velocidade_angular = dblarr(NDIM)
  icore = dblarr(NDIM)
  ienv  = dblarr(NDIM)
  itotal = dblarr(NDIM)
  Jtotal = dblarr(NDIM)
  age   = dblarr(NDIM)
  omega_core = dblarr(NDIM)
  omega_env = dblarr(NDIM)
  delta_t = dblarr(NDIM)
  age_myr = dblarr(NDIM)
  delta_j = dblarr(NDIM)
  omega_core_primeiro_termo = dblarr(NDIM)
  omega_core_segundo_termo = dblarr(NDIM)
  omega_core_terceiro_termo = dblarr(NDIM)
  omega_env_primeiro_termo = dblarr(NDIM)
  omega_env_segundo_termo = dblarr(NDIM)
  omega_env_terceiro_termo = dblarr(NDIM)
  raio_nucleo = dblarr(NDIM)
  raio_estrela = dblarr(NDIM)
  massa_estrela = dblarr(NDIM)
  massa_nucleo = dblarr(NDIM)
  luminosidade = dblarr(NDIM)
  

  age2 = dblarr(NDIM2)
  icore2 = dblarr(NDIM2)
  ienv2 = dblarr(NDIM2)
  massa_nucleo2 = dblarr(NDIM2)
  raio_nucleo2 = dblarr(NDIM2)
  massa_estrela2 = dblarr(NDIM2)
  raio_nucleo2 = dblarr(NDIM2)
  omega_env2 = dblarr(NDIM2)
  omega_core2 = dblarr(NDIM2)
  delta_t2 = dblarr(NDIM2)
  delta_j2 = dblarr(NDIM2)
  luminosidade2 = dblarr(NDIM2)
  omega_core_primeiro_termo2 = dblarr(NDIM2)
  omega_core_segundo_termo2 = dblarr(NDIM2)
  omega_core_terceiro_termo2 = dblarr(NDIM2)
  omega_env_primeiro_termo2 = dblarr(NDIM2)
  omega_env_segundo_termo2 = dblarr(NDIM2)
  omega_env_terceiro_termo2 = dblarr(NDIM2)
  jponto2 = dblarr(NDIM2)
  pp = dblarr(NDIM2)


  readf,1,line   ; apenas o cabeçalho
  readf,1,data   ; e agora, o resto do arquivo
  free_lun,1     ; liberando a unidade lógica

  
;  lixo = ' '
;
;  file = 'C:\Users\Emilly\OneDrive\IDL\Dados\dados_10_massa.txt'
; ; close,1
;
;  openr,1,file
;  readf,1,lixo
;  readf,1,data
;  close,1

  age_year = 10.^data[1,*]
  age_myr = age_year/myr
  age_segundo = age_year * ano_segundo
  
  deltamax = age_year[NDIM-1] - age_year[0]
  

  
  for i = 0, NDIM2 - 1 do begin
    age2[i] = age_year[0] + deltamax*float(i)/float(NDIM2-1)
  endfor

  for i = 0L, NDIM-1L do begin
    massa_nucleo[i] = data[9,i]*Msun
    raio_estrela[i] = data[5,i]*Rsun
    massa_estrela[i] = data[0,i]*Msun
    raio_nucleo[i] = data[10,i]*Rsun
    icore[i] = ((data[12,i])^2.)*data[0,i]*Msun*(data[5,i]*Rsun)^2.
    ienv[i] = ((data[11,i])^2.)*data[0,i]*Msun*(data[5,i]*Rsun)^2.
    luminosidade[i] = data[3,i]
  endfor
  
 
  icore2 = interpol(icore,age_year,age2)
  ienv2 = interpol(ienv,age_year,age2)
  massa_nucleo2 = interpol(massa_nucleo,age_year,age2)
  raio_estrela2 = interpol(raio_estrela,age_year,age2)
  massa_estrela2 = interpol(massa_estrela,age_year,age2)
  raio_nucleo2 = interpol(raio_nucleo,age_year,age2)
  luminosidade2 = interpol(luminosidade,age_year,age2)
  
  age2_segundo = age2 * ano_segundo

  velocidade_angular_inicial = 2.*!pi/(P0*periodo_segundos)
  ;tau_ce_anos = 30.d+6
;  tau_ce_segundos = tau_ce_anos*365.*24.*60.*60.

  omega_env[0] = velocidade_angular_inicial
  omega_core[0] = 0.
  omega_env2[0] = velocidade_angular_inicial
  omega_core2[0] = 0.
  k1 = 22.592*data[0,0]^2. - 47.504*data[0,0] + 26.602
  
  if (jj eq 0) then tcc = 500.d+6
  if (jj eq 1) then tcc = 80.d+6
  if (jj eq 2) then tcc = 30.d+6

  tau_ce_segundos = tcc*365.*24.*60.*60.


  icount = -1
  for i = 1,NDIM2-1L do begin
    delta_t2[i] = (age2_segundo[i] - age2_segundo[i-1])
    delta_j2[i] = ((ienv2[i-1] * icore2[i-1])/(icore2[i-1] +ienv2[i-1])) * (omega_core2[i-1] - omega_env2[i-1])


    Protday = (2.*!pi/omega_env2[i-1])/(periodo_segundos)
    
    aux1 = massa_estrela2[0]/Msun
    aux2 = raio_estrela2[i]/Rsun
    aux3 = luminosidade2[i]
    
    pp[i] = Protday

    BOREAS4,aux1,aux2,aux3,Protday,0.0, $
      Rossby,fstar,Mdot,Mdot_hot,Mdot_cold,Bphoto

    Bphoto *= fstar
    Mdot *= Msun/(365.*24.*3600.)
    Vesc2  = (2.*Gconst*data[0,0]*Msun)/raio_estrela2[i]
    
;    k1 = 22.592*(massa_estrela2/Msun)^2. - 47.504*massa_estrela2/Msun + 26.602
    r_a_zero = sqrt(((0.0506^2.)*Vesc2) + (omega_env2[i-1]^2.)*(raio_estrela2[i-1]^2.))
    r_a_um = ((Bphoto^2.)*(raio_estrela2[i-1]^2.))/(r_a_zero*Mdot)
    r_a = (r_a_um^0.22)*raio_estrela2[i-1]*k1
;    r_a = (r_a_um^0.22) * aux2 * (1.8 * (Msun / aux1)^0.5)
    
;    m_atual = massa_estrela2[i-1]
;    r_atual = raio_estrela2[i-1]
;    
;    ; Cálculo com o ajuste de massa
;    fator_m = 1.8 * (Msun / m_atual)^(0.5)
;    r_a = (r_a_um^0.22) * r_atual * fator_m

    Jponto = omega_env2[i-1]*Mdot*(r_a^2.)
;    jponto2[i] = omega_env2[i-1]*Mdot*(r_a^2.)

    if (icore2[i] eq 0.) then begin
      omega_core2[i] = omega_env2[0]
      if (age2_segundo[i] LT tau_d) then begin
        omega_env2[i] = omega_env2[0]
      endif else begin
        omega_env2[i] = (omega_env2[i-1] - ((omega_env2[i-1])/ienv2[i])*(ienv2[i]-ienv2[i-1])) - Jponto*delta_t2[i]/ienv2[i-1] ; - jponto*dt/Ienv
      endelse

    endif else begin

      icount += 1
      omega_core_primeiro_termo2[i] = delta_t2[i]*delta_j2[i-1] / (icore2[i] * tau_ce_segundos)
      omega_core_segundo_termo2[i] = (((2. / 3.) * (raio_nucleo2[i]^2.)) / icore2[i]) * omega_env2[i - 1] * (massa_nucleo2[i] - massa_nucleo2[i - 1])
      omega_core_terceiro_termo2[i] = ((omega_core2[i - 1]) / icore2[i]) * (icore2[i] - icore2[i - 1])
     
      omega_core2[i] = omega_core2[i - 1] - omega_core_primeiro_termo2[i] + omega_core_segundo_termo2[i] - omega_core_terceiro_termo2[i]
 
      if (age2_segundo[i] LT tau_d) then begin
        omega_env2[i] = omega_env2[0]
      endif else begin

        omega_env_primeiro_termo2[i] = delta_t2[i]*(delta_j2[i-1] / (ienv2[i] * tau_ce_segundos))
        omega_env_segundo_termo2[i] = (((2. / 3.) * (raio_nucleo2[i - 1]^2.)) / ienv2[i]) * omega_env2[i - 1] * (massa_nucleo2[i] - massa_nucleo2[i - 1])
        omega_env_terceiro_termo2[i] = ((omega_env2[i - 1]) / ienv2[i]) * (ienv2[i] - ienv2[i - 1])
        omega_env2[i] = omega_env2[i - 1] + omega_env_primeiro_termo2[i] - omega_env_segundo_termo2[i] - omega_env_terceiro_termo2[i] - Jponto*delta_t2[i]/ienv2[i-1] ; - jponto*dt/Ienv
      endelse
    endelse

   ; Jtotal1[i] = (icore2[i]*omega_core2[i]+ienv2[i]*omega_env2[i])/Jtotal[0]
    
  ;  print, omega_env2[i] 
   

  endfor

;  print, jponto2
;  print, jponto
;  stop


  cc = plot(age2/myr, omega_core2/velo_solar, /xlog, /ylog, $
    xtitle='Idade (Manos)', ytitle='$\Omega/\Omega_\odot$', $
    yrange=[0.1,30.], name = 'Nucleo', color = 'black', $
      title=titlemass[jj], thick= 3)

  rr = plot(age2/myr, omega_env2/velo_solar,color=[255, 40, 150], /current, /overplot, name = 'Envelope',LINESTYLE='--', thick= 3)

  l = LEGEND(TARGET=[cc, rr], POSITION=[0.4, 0.8], font_size=10)
  
  endfor
  
end


pro constants
  common CNT, Msun, Rsun, P0, myr, tau_d, periodo_segundos,velo_solar, Gconst, ano_segundo, Lsun

  Lsun   = 3.826d+33
  Msun = 1.98847d+33
  Rsun = 6.957d+10
  P0 = 10.
  myr = 1.d+6
  tau_d = 5.*myr*365.*24.*60.*60.
  periodo_segundos = 86400.
  ano_segundo = 365.*24.*60.*60.
  velo_solar = 2.87e-6
  Gconst = 6.6732d-8

  return
end