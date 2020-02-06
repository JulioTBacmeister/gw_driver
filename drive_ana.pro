pro drive_ana,xtag=xtag,tauxy=tauxy,lon=lon,lat=lat

start=[2010L,10L,1,0]
final=[2010L,11L,1,0]


date=start
stopcon=0
dout$ = '/project/amp/juliob/OGW_UnitTest_output/beta07/'
;xtag='AniOutput'

; double check

itime=0
while stopcon eq 0 do begin


     y$=padstring( date(0) )
     m$=padstring( date(1) , /e1 )
     d$=padstring( date(2) , /e1 )
     s$=padstring( date(3) , /e4 )

     out$  =  dout$+xtag+'.'+y$+'-'+m$+'-'+d$+'-'+s$


     rdgw3d,fn=out$,tauf=tauf,utgw=utgw,vtgw=vtgw,lon=lon,lat=lat

     print,'  Got ',out$

     if itime eq 0 then begin
        s=size(tauf)
        tauxy = dblarr( s(1),s(2), 31L*4 )
     endif

      tauxy(*,*,itime) = tauf(*,*,46,0)

     date = bump_time_vec(date,dtvc=[0,0,0,6*3600L],/sssss)
     xx=total( abs(date-final ) )
     if xx eq 0 then stopcon=1

     itime=itime+1

  endwhile

STOP
return
end
