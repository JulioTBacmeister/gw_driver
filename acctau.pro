pro acctau,mon=mon $ 
          ,use_cesm2_0_ridges=use_cesm2_0_ridges $
          ,use_devel_ridges=use_devel_ridges $
          ,use_isotropic=use_isotropic $
          ,tauc=tauc,lat=lat,lon=lon,zi=zi,zm=zm

start=long( [mon(0),mon(1),1,0] )
final=long( [mon(0),mon(1)+1,1,0] )
;final=long( [mon(0),mon(1),1,0] )

;start=[2010L,11L,1,0]
;final=[2010L,12L,1,0]


   if keyword_set(use_devel_ridges) then begin
      lncom = '/project/amp/juliob/Topo-generate-052416/Topo/cube_to_target_devel/output/fv_0.9x1.25_nc3000_Nsw042_Nrs008_Co060_Fi001_ZR_test_v3.nc ridgedata.nc'
      xtag='GW-devel'
      ridge=1
   endif
   if keyword_set(use_cesm2_0_ridges) then begin
      lncom = '/project/amp/juliob/scam/bndtopo/fv_0.9x1.25_nc3000_Nsw042_Nrs008_Co060_Fi001_ZR.nc ridgedata.nc'
      xtag='GW-cesm-2-0'
      ridge=1
   endif
   if keyword_set(use_isotropic) then begin
      lncom = '/project/amp/juliob/scam/bndtopo/fv_0.9x1.25_nc3000_Nsw042_Nrs008_Co060_Fi001_ZR.nc ridgedata.nc'
      xtag='GW-isotropic'
      ridge=0
   endif
   isotropic = 1-ridge

   dmet$ = '/project/amp/juliob/ERAinterim/ERAI_fv09x1.25/2010/'
   xtag=xtag+'-ERAI.'

     print,' Experiment tag ' 
     print,xtag
     print,' Dates '
     print,start,'->',final

; check things


date=start
stopcon=0
ird=0
while stopcon eq 0 do begin


     y$=padstring( date(0) )
     m$=padstring( date(1) , /e1 )
     d$=padstring( date(2) , /e1 )
     s$=padstring( date(3) , /e4 )


     outdir$='/project/amp/juliob/OGW_UnitTest_output/'

     out$  =  outdir$+xtag+y$+'-'+m$+'-'+d$+'-'+s$+'.dat'
      
     rdgw,f=out$,tauf=tau,ridge=ridge,lat=lat,lon=lon,zi=zi,isot=isotropic

     if keyword_set(ridge) then begin
       if ird eq 0 then begin
         s=size(tau)
         ; DIMS tau = (lon,lat,lev,irdg)
         tauc=dblarr( s(1),s(2),s(3),s(4), 31*4*1L )
       endif
       tauc(*,*,*,*,ird) = tau    
     endif else begin
       if ird eq 0 then begin
         s=size(tau)
         ; DIMS tau = (lon,lat,lev,irdg)
         tauc=dblarr( s(1),s(2),s(3), 31*4*1L )
       endif
       tauc(*,*,*,ird) = tau    
     endelse


     xx=total( abs(date-final ) )
     if xx eq 0 then stopcon=1
     date = bump_time_vec(date,dtvc=[0,0,0,6*3600L],/sssss)

     print," Got ",out$
     ird=ird+1

  endwhile

  if keyword_set(ridge) then begin
  tauc=tauc(*,*,*,*, 0:ird-1 )
  endif else begin
  tauc=tauc(*,*,*, 0:ird-1 )
  endelse

end
