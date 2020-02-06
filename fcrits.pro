pro fcrits
;-----------------------------------------
;  v1 parameters in CESM2.0 -11/23/16
;-----------------------------------------
; real(r8), parameter ::  alpha_Fr2 =2.0_r8 
; real(r8), parameter ::  Frx_crit  =2.0_r8 
; real(r8), parameter ::  Fr1_crit  =1.0_r8 

 alpha_Fr2 =2.0d
 Frx_crit  =2.0d 
 Fr1_crit  =1.0d 


 hdsp = 3000.-dindgen(3000)
 anixy = 1.0+dblarr(3000)
 nsrc = dblarr(3000) + 0.013532387d
 ubmsrc = dblarr(3000) + 10.0d
 Fr1 = dblarr(3000)
 Fr2 = dblarr(3000)
 Frx = dblarr(3000)
  
 mxdis = hdsp
 

  Fr1(*) = Fr1_crit  ; 1.00_r8
  Frx(*) = hdsp(*)*nsrc(*)/abs( ubmsrc(*) )

  ; SM-like param of high-drag regime: Peaks for Frx 
  ; around 2 and decreases for higher topo.
  ;------------------------------------------------
  Fr2(*) = Fr1(*) + Fr1(*)* alpha_Fr2 *anixy(*)

  Fr2oo= Fr2

  for i=0,3000-1 do begin
;  where(Frx > Frx_crit )
;        Fr2(*) = Fr1(*) + Fr1(*)* alpha_Fr2 *anixy(*) &
;             * (2*Frx_crit -Fr1_crit - Frx(*)) &
;             / (Frx_crit - Fr1_crit )
;  endwhere
    if (Frx(i) gt Frx_crit ) then begin
        Fr2(i) = Fr1(i) + Fr1(i)* alpha_Fr2 *anixy(i) $
              * (2*Frx_crit -Fr1_crit - Frx(i)) $
              / (Frx_crit - Fr1_crit )
    endif
  endfor

  for i=0,3000-1 do begin
;  where(Fr2 < Fr1)
;     Fr2=Fr1
;  endwhere
    if (Fr2(i) lt Fr1(i) ) then Fr2(i)=Fr1(i)
 endfor


 ddw = Fr2 * ubmsrc / nsrc
 tlb = mxdis - ddw
 oo=where(tlb lt 0)
 if min(oo) ne -1 then tlb(oo)=0.d

 dwv  = Fr1 * ( abs(ubmsrc) )/nsrc
 bwv = mxdis - dwv
 oo=where(bwv lt 0)
 if min(oo) ne -1 then dwv(oo)=hdsp(oo)



;========================================
; Following code appears to reproduce 
;   tau_free + tau_d
; in Scinocca & McFarlane (2000; Eq. 28)
; where:
;   tau_free <==> tauoro
;   tau_d    <==> taudsw
;=======================================
 taulin = nsrc * ubmsrc * (hdsp^2)
 tauoro = nsrc * ubmsrc * (dwv^2)

 beta = taulin*0

 Bmax=2.0
 oo1 = where( Frx ge 1. and Frx lt 1.5 )
 beta(oo1) = 2.*Bmax*(Frx(oo1)-1)
 oo2 = where( Frx ge 1.5 and Frx lt 3.0 )
 beta(oo2) = ( 1.+ Bmax - (1./1.5)^2 ) * ( ( ( 3.-Frx(oo2) )/1.5)^2 ) + (1./Frx(oo2))^2 -1.
 
 taudsw = taulin*0

 oo3 = where( Frx ge 1. and Frx lt 3.0 )
 taudsw(oo3) = ( 1.+beta(oo3))*taulin(oo3) - tauoro(oo3)
 


STOP

return
end
