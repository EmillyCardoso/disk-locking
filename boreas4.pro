
;-----------------------------------------------------------------------------
    PRO BOREAS4,Mstar,Rstar,Lstarr,Protday,FeH, $
               Rossby,fstar,Mdot,Mdot_hot,Mdot_cold,Bphoto;, accrate, P0, 
;-----------------------------------------------------------------------------
;
; boreas.pro:  Stand-alone IDL procedure to compute the mass loss rate of
;  a cool, late-type star using the models of Cranmer & Saar (2011).
;
; Inputs (all scalars):
;
; Mstar:  mass of star (in units of solar mass)
; Rstar:  radius of star (in units of solar mass)
; Lstar:  bolometric luminosity of star (in units of solar mass)
; Protday:  rotation period of star (in days)
; FeH:  metallicity (log of iron/hydrogen abundance ratio, in solar units)
;
; Outputs: (all scalars):
;
; Rossby:  dimensionless Rossby number for this star
; fstar:  dimensionless filling factor of open flux tubes at photosphere
; Mdot:  mass loss rate (in units of solar masses per year)
; Mdot_hot:  mass loss rate from just hot coronal (gas pressure) model
; Mdot_cold:  mass loss rate from just cold (wave pressure) model
; MATR:  hot wind's Alfvenic Mach number at the transition region
;
;-----------------------------------------------------------------------------

; Are the inputs all there?


    if (N_params() ne 11) then begin
      print,' '
      print,' Syntax: '
      print,' '
      print,' BOREAS,Mstar,Rstar,Lstar,Protday,FeH, $ '
      print,'    Rossby,fstar,Mdot,Mdot_hot,Mdot_cold,MATR '
      print,' '
      print,' See comments in boreas.pro for meanings of parameters. '
      print,' '
      return
    endif

; First define useful constants

    Gconst  = 6.6732d-08
    xMsun   = 1.989d+33
    xRsun   = 6.96d+10
    xLsun   = 3.826d+33
    boltzK  = 1.380622d-16
    xmHyd   = 1.67333d-24
    stefan  = 5.66961d-5
    xMdotyr = 6.30276d+25

; Some of the key parameters for the models are defined here

    alphaturb  = 0.5
    heightfac  = 0.5
    theta      = 0.333
    ffloor     = 1.0d-4
    ellperpSUN = 3.0d7
    niter      = 50

; Set up basic stellar parameters in cgs units


    xMcgs  = Mstar * xMsun
    xRcgs  = Rstar * xRsun
    Lstar = 10.^Lstarr
    xLcgs  = Lstar * xLsun
    ;print, xLcgs,xRcgs
    Teff   = (xLcgs/4./!PI/xRcgs/xRcgs/stefan)^0.25
    grav   = Gconst*xMcgs/xRcgs/xRcgs
    logg   = alog10(grav)
    ZoZsun = 10.^(FeH)

    Vesc2  = 2.*Gconst*xMcgs/xRcgs
    Vesc   = sqrt(Vesc2)

    TeffSUN    = 5770.2
    ProtdaySUN = 25.3
    gravSUN    = Gconst*xMsun/xRsun/xRsun

; Make sure parameters aren't too crazy

    icrazy = 0
    if ((Teff lt 1500.) or (Teff gt 12000.)) then icrazy = 1
    if ((logg lt -4.0) or (logg gt 7.0)) then icrazy = 2
    if ((Lstar lt 1.d-4) or (Lstar gt 1.d6)) then icrazy = 3
    if ((Mstar lt 0.001) or (Mstar gt 100.)) then icrazy = 4
    if ((FeH lt -5.0) or (FeH gt 2.0)) then icrazy = 5

    if (icrazy ne 0) then begin
      print,' '
      print,' Input parameters seem to be out of bounds! '
      print,icrazy
      ;print,xLcgs, Teff,logg,Lstar,Mstar,Feh
      stop
      print,' '
      return
    endif

; Estimate Rossby number and open-flux filling factor

    tauc    = 314.241*exp(-Teff/1952.5)*exp(-(Teff/6250.)^18.)+0.002
    taucSUN = 314.241*exp(-TeffSUN/1952.5)*exp(-(TeffSUN/6250.)^18.)+0.002

    gratio = gravSUN/grav
    if (gratio gt 1.) then tauc = tauc * (gratio^0.18)

    Rossby    = Protday / tauc
    RossbySUN = ProtdaySUN / taucSUN
    Ronorm    = Rossby/RossbySUN
    fstar     = 0.55 / (1. + (Ronorm/0.16)^2.3)^1.22

    if (fstar lt ffloor) then fstar = ffloor

; Compute photospheric mass density

    ngrid    = 18
    Teffgrid = dindgen(ngrid)/double(ngrid-1)*(6000.-2500.)+2500.
    C0fit = [-9.0329342,-8.8057005,-8.6187019,-8.5343736,-8.5792239,$
             -8.6915425,-8.8033701,-8.8791533,-8.9288180,-8.9604793,$
             -8.9954977,-9.0593624,-9.1566287,-9.2743908,-9.4120161,$
             -9.5781877,-9.7290674,-9.8972636]
    C1fit = [0.78081830,0.68284416,0.60471646,0.56629124,0.57113263,$
             0.59584083,0.61352883,0.61218030,0.59646729,0.57949132,$
             0.56963417,0.57219919,0.58595944,0.60671743,0.63103575,$
             0.65574884,0.67753323,0.69808401]

    C0now = INTERPOL(C0fit,Teffgrid,[Teff])
    C1now = INTERPOL(C1fit,Teffgrid,[Teff])

    rhophoto = 10.^(C0now(0) + C1now(0)*logg)

; Compute the photospheric equation of state (from fit to OPAL models),
; magnetic field, and scale height

    xmuavg = 1.75 + 0.5*tanh((3500.-Teff)/600.)
    gadia  = 5./3.

    Pphoto = rhophoto*boltzK*Teff/(xmuavg*xmHyd)
    Bequi  = sqrt(8.*!PI*Pphoto)
    Bphoto = 1.13*Bequi
    
    
    

    
    cs2    = gadia*boltzK*Teff/(xmuavg*xmHyd)
    Hphoto = cs2 / (gadia*grav)

    xmuavgSUN = 1.75 + 0.5*tanh((3500.-TeffSUN)/600.)
    cs2SUN    = gadia*boltzK*TeffSUN/(xmuavgSUN*xmHyd)
    HphotoSUN = cs2SUN / (gadia*gravSUN)

; Estimate surface flux of Alfven waves using a parameterized fit to the
; Musielak et al. (2002a) kink-mode flux models

    alpha_MU02 = 6.774 + 0.5057*logg
    T0_MU02    = 5624. + 600.2*logg
    F0_MU02    = exp(22.468 - 0.0871*logg)

    FluxAphoto = F0_MU02 * ((Teff/T0_MU02)^alpha_MU02) $
               * exp(-(Teff/T0_MU02)^25.)    
;    FluxAphoto /= 2.5
    if (Mstar ge 0.8) then FluxAphoto/=2.5
    if (Mstar eq 0.5) then FluxAphoto/=1.5
    if (FluxAphoto lt 1.d-10) then FluxAphoto = 1.d-10

; Set up MHD turbulence parameters at the photosphere

    Valfphoto    = Bphoto / sqrt(4.*!PI*rhophoto)
    vperpPHOT    = sqrt(FluxAphoto/rhophoto/Valfphoto)
    ellperpphoto = ellperpSUN * (Hphoto/HphotoSUN)

    Qphoto       = alphaturb*rhophoto*(vperpPHOT^3)/ellperpphoto

; Extrapolate MHD turbulence parameters up to the transition region (TR)

    fTR  = fstar^theta
    BTR  = Bphoto * (fstar/fTR)
    uinf = Vesc

    TTR       = 2.0d5
    Lammax    = 7.4d-23 + 4.2d-22*(ZoZsun^1.13)
    LammaxSUN = 7.4d-23 + 4.2d-22

; Iterate to find self-consistent solution for density and heating rate
; at the TR, assuming that the non-WKB reflection at the TR is imperfect.
; The reflection coefficient is given by an approximate version of the
; low-frequency Cranmer (2010) limit.

    ReflCoef = 0.5

    for iter=1,niter do begin
      quench  = ReflCoef*(1.+ReflCoef)/(1.+ReflCoef^2)^1.5 * sqrt(2.)
      sqbrac  = quench*Qphoto*xmHyd*xmHyd/(rhophoto^0.25)/Lammax
      rhoTR   = (sqbrac^(4./7.)) * (fstar^(2.*(1.-theta)/7.))
      QTR     = Qphoto*quench*((rhoTR/rhophoto)^0.25)*sqrt(BTR/Bphoto)
      ValfTR  = BTR / sqrt(4.*!PI*rhoTR)
      PTR     = 2.0*rhoTR*boltzk*TTR/xmHyd
      ReflCoefnew = abs((ValfTR-uinf)/(ValfTR+uinf))
      ReflCoef    = sqrt(ReflCoefnew*ReflCoef)
    endfor

; Does the heating-related energy flux at the TR exceed the flux in
; "passive propagation" of the Alfven waves?  If so, cap it!

    FluxTR0  = heightfac*QTR*xRcgs
    FluxA_TR = FluxAphoto * fstar/fTR

    if (FluxTR0 gt FluxA_TR) then FluxTR0 = FluxA_TR

; Estimate the mass loss rate for a hot coronal wind, using the
; Hansteen et al. (1995) energy balance approximation.

    velfacTR = 1.4d6 * sqrt(Lammax/LammaxSUN)
    FcondTR  = PTR * velfacTR
    fraccond = FcondTR / FluxTR0
    fccmax   = 0.9
    if (fraccond gt fccmax) then fraccond = fccmax
    FluxTR   = FluxTR0 * (1.-fraccond)

    AreaTR      = fTR * (4.*!PI*xRcgs*xRcgs)
    Mdotcgs_hot = AreaTR*FluxTR/ ((Vesc2+uinf*uinf)/2.)
    Mdot_hot    = Mdotcgs_hot / xMdotyr

    uTR_hot = Mdotcgs_hot / (rhoTR*AreaTR)
    MATR    = uTR_hot / ValfTR

; For the Holzer et al. (1983) cold wave-driven wind, first estimate the
; radius, speed, and magnetic field at the wave-modified critical point.

    enn   = 2.0       ; B \propto r^{-enn} at crit point
    beta  = 0.5*enn
    eee36 = (3./7./beta + 4./7.) / (1. + (vperpPHOT/Vesc)^2)

    rcrit = 1.75*eee36 * xRcgs
    ucrit = sqrt(Gconst*xMcgs/enn/rcrit)

    Bcrit = Bphoto * ((xRcgs/rcrit)^2) * fstar

; At the critical point, iterate on the definition of u_crit and the
; condition of wave action conservation to get density and wave amplitude.

    Areaphoto    = fstar * (4.*!PI*xRcgs*xRcgs)
    Areacrit     = 4.*!PI*rcrit*rcrit
    action_photo = rhophoto*(vperpPHOT^2)*Valfphoto*Areaphoto

    vperpcrit = 2.*ucrit
    rhocrit   = 4.*!PI*(action_photo/((vperpcrit^2)*Bcrit*Areacrit))^2

    for iter=1,niter do begin
      Valfcrit  = Bcrit / sqrt(4.*!PI*rhocrit)
      MAcrit    = ucrit / Valfcrit
      macfac    = (1. + 3.*MAcrit)/(1. + MAcrit)
      vperpcrit = 2.*ucrit / sqrt(macfac)
      rhocrit   = action_photo / $
                  ((vperpcrit^2)*Valfcrit*Areacrit*(1.+MAcrit)^2)
    endfor

    Mdotcgs_cold = rhocrit*ucrit*Areacrit
    Mdot_cold    = Mdotcgs_cold / xMdotyr

; Estimate the actual mass loss from both hot and cold processes.

    Mdot = Mdot_cold + (Mdot_hot*exp(-4.*MATR^2))

    return
    end

