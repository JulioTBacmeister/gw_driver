start=[2010L,11L,1,0]
final=[2010L,12L,1,0]

;start=[2010L,11L,15,0]
;final=[2010L,11L,16,0]

use_cesm2_0_ridges=0
use_devel_ridges=0
use_isotropic=1

   if keyword_set(use_devel_ridges) then begin
      lncom = 'ln -sf /project/amp/juliob/Topo-generate-052416/Topo/cube_to_target_devel/output/fv_0.9x1.25_nc3000_Nsw042_Nrs008_Co060_Fi001_ZR_test_v3.nc ridgedata.nc'
      nmlst = 'ridge.nml'
      xtag='GW-devel'
   endif
   if keyword_set(use_cesm2_0_ridges) then begin
      lncom = 'ln -sf /project/amp/juliob/scam/bndtopo/fv_0.9x1.25_nc3000_Nsw042_Nrs008_Co060_Fi001_ZR.nc ridgedata.nc'
      nmlst = 'ridge.nml'
      xtag='GW-cesm-2-0'
   endif
   if keyword_set(use_isotropic) then begin
      lncom = 'ln -sf /project/amp/juliob/scam/bndtopo/fv_0.9x1.25_nc3000_Nsw042_Nrs008_Co060_Fi001_ZR.nc ridgedata.nc'
      nmlst = 'isotropic.nml'
      xtag='GW-isotropic-s1'
   endif
   dmet$ = '/project/amp/juliob/ERAinterim/ERAI_fv09x1.25/2010/'
   xtag=xtag+'-ERAI.'

     print,'topo file: '
     print,' ',lncom
     print,'namelist: ' 
     print,' ',nmlst
     print,'metdata: ' 
     print,' ',dmet$
     print,'Experiment tag: ' 
     print,' ',xtag
     print,'Dates: '
     print,' ',start,'->',final

; check things
STOP

spawn,'mkdir -p output'
spawn,lncom
spawn,'ln -sf '+nmlst+' control.nml'


date=start
stopcon=0
dout$ = '/project/amp/juliob/OGW_UnitTest_output/'


spawn,'ls -l control.nml'
spawn,'more control.nml'
spawn,'ls -l ridgedata.nc'

; double check
STOP

while stopcon eq 0 do begin


     y$=padstring( date(0) )
     m$=padstring( date(1) , /e1 )
     d$=padstring( date(2) , /e1 )
     s$=padstring( date(3) , /e4 )

     fmet$ =  'ERAI_fv09_L30.cam2.i.'+y$+'-'+m$+'-'+d$+'-'+s$+'.nc'
     out$  =  dout$+xtag+y$+'-'+m$+'-'+d$+'-'+s$+'.dat'

     lncom = 'ln -sf '+dmet$+fmet$+' metdata.nc'
     print,lncom
     spawn,lncom

     xcom = './test.x'
     spawn,xcom

     mvcom = 'mv fort.511 '+out$
     print,mvcom
     spawn,mvcom


     xx=total( abs(date-final ) )
     if xx eq 0 then stopcon=1
     date = bump_time_vec(date,dtvc=[0,0,0,6*3600L],/sssss)

  endwhile


end
