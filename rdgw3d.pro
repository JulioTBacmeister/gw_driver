pro rdgw3d,fn=fn,ridge=ridge,isotropic=isotropic $
    ,lon=lons,lat=lats,tauf=tauf,utgw=utgw,vtgw=vtgw $
    ,endstop=endstop

if not keyword_set(fn) then fn='fort.511'

close,1
openr,1,/f77_u,fn

        nlon=0L
        nlat = 0L 
        pver = 0L
        nrdg=0L
        ntim=0L
        itime=-1L

  while itime lt ntim do begin
        readu,1,nlon,nlat,pver,nrdg,itime,ntim
        statelat = dblarr( nlon * nlat )
        lons    = dblarr( nlon )
        lats    = dblarr( nlat )
        readu,1,lons,lats,statelat

        tau_o = dblarr( nlon * nlat, pver+1 ) 
        tauf  = dblarr( nlon * nlat, pver+1, nrdg ) 
        utgw  = dblarr( nlon * nlat, pver )
        vtgw  = dblarr( nlon * nlat, pver )

        for irdg = 0,nrdg-1 do begin
               ; After gw_drag_prof
          readu,1,tau_o  
             tauf(*,*,irdg) = tau_o
        endfor     

             readu,1,utgw
             readu,1,vtgw
     
         
             tauf = reform( tauf,nlon,nlat,pver+1,nrdg )
             utgw = reform( utgw,nlon,nlat,pver )
             vtgw = reform( vtgw,nlon,nlat,pver )

    endwhile

    if keyword_set(stops) then STOP

return
end
