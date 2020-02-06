pro rdgw,fn=fn,ridge=ridge,isotropic=isotropic $
    ,tauf=tau_stg3 ,tau0=tau_stg0, tau2=tau_stg2 $
    ,tau1=tau_stg1, tauf=tauf $
    ,tauoro=tauoro,taudsw=taudsw,taulin=taulin,tausoro=tau2oro $
    ,frx=frxo,fr1=fr1o,fr2=fr2o $
    ,wbr=wbr,tlb=tlb,bwv=bwv,hdspwv=hdspwv,hdspdw=hdspdw $     
    ,zm=zm,zi=zi,uu=u,v=v,t=t,ubm=ubm,nm=nm $
    ,lat=lats,lon=lons $
    ,mxdis=mxdis,sgh=sgh,skip_reform=skip_reform,stops=stops,old=old $
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

if keyword_set(ridge) then begin
  while itime lt ntim do begin
        readu,1,nlon,nlat,pver,nrdg,itime,ntim
        statelat = dblarr( nlon * nlat )
        lons    = dblarr( nlon )
        lats    = dblarr( nlat )
        readu,1,lons,lats,statelat

        taurx0 = dblarr( nlon * nlat, pver+1 )
        taury0 = dblarr( nlon * nlat, pver+1 )
        taurm0 = dblarr( nlon * nlat, pver+1 )

        tau_stg0 = dblarr( nlon * nlat, pver+1, nrdg )
        tau_stg1 = dblarr( nlon * nlat, pver+1, nrdg ) 
        tau_stg2 = dblarr( nlon * nlat, pver+1, nrdg ) 
        tau_stg3 = dblarr( nlon * nlat, pver+1, nrdg ) 
        tau_stg4 = dblarr( nlon * nlat, pver+1, nrdg ) 
        tau_o    = dblarr( nlon * nlat, pver+1 ) 

        tauoro = dblarr( nlon * nlat )
        tau2oro = dblarr( nlon * nlat )
        taudsw = dblarr( nlon * nlat )
        hdspwv = dblarr( nlon * nlat )
        hdspdw = dblarr( nlon * nlat )

        ubmsrc   = dblarr( nlon * nlat )
        tlb    = dblarr( nlon * nlat )
        bwv    = dblarr( nlon * nlat )
        wbr    = dblarr( nlon * nlat )
        wbrx    = dblarr( nlon * nlat )
        trn    = dblarr( nlon * nlat )
        usrc   = dblarr( nlon * nlat )
        vsrc   = dblarr( nlon * nlat )
        nsrc   = dblarr( nlon * nlat )
        rsrc   = dblarr( nlon * nlat )
        mxdis   = dblarr( nlon * nlat, nrdg )
        mxdiso  = dblarr( nlon * nlat )
        Fr1o   = dblarr( nlon * nlat )
        Fr2o   = dblarr( nlon * nlat )
        Frxo   = dblarr( nlon * nlat )


        angllo = dblarr( nlon * nlat )
        anixyo = dblarr( nlon * nlat )
        hwdtho = dblarr( nlon * nlat )
        clngto = dblarr( nlon * nlat )
        kwvrdgo = dblarr( nlon * nlat )
        gbxaro = dblarr( nlon * nlat )

        m2_diag    = dblarr( nlon * nlat, pver ) 

        effgw   = dblarr( nlon * nlat )
        effgw_x = dblarr( nlon * nlat )
 
        src_level = lonarr( nlon * nlat )
        tend_level = lonarr( nlon * nlat )
        bwv_level = lonarr( nlon * nlat )
        tlb_level = lonarr( nlon * nlat )
        
        v = dblarr( nlon * nlat, pver )
        u = dblarr( nlon * nlat, pver )
        t = dblarr( nlon * nlat, pver )
        nm   = dblarr( nlon * nlat, pver, nrdg )
        nm_o   = dblarr( nlon * nlat, pver )
        zm   = dblarr( nlon * nlat, pver )
        ubm  = dblarr( nlon * nlat, pver, nrdg )
        ubi  = dblarr( nlon * nlat, pver+1, nrdg )
        ubm_o  = dblarr( nlon * nlat, pver )
        ubi_o  = dblarr( nlon * nlat, pver+1 )
        pifc = dblarr( nlon * nlat, pver+1 )
        pmid = dblarr( nlon * nlat, pver )
        zi   = dblarr( nlon * nlat, pver+1 )
        utgw = dblarr( nlon * nlat, pver )
        vtgw = dblarr( nlon * nlat, pver )
        gwut = dblarr( nlon * nlat, pver )


        for irdg = 0,nrdg-1 do begin
           ;   After gw_rdg_src
          if keyword_set(old) then readu,1,mxdiso,angllo,anixyo,hwdtho,clngto,kwvrdgo,gbxaro
          if not keyword_set(old) then readu,1,mxdiso,angllo,anixyo,hwdtho,clngto,kwvrdgo,gbxaro
             mxdis(*,irdg) = mxdiso
          readu,1,ubmsrc,usrc,vsrc,nsrc, tlb, bwv, Fr1o, Fr2o, Frxo,rsrc
          readu,1,src_level,tend_level, bwv_level, tlb_level
          readu,1,tau_o,ubi_o,ubm_o,nm_o,zm,zi,u,v,t,pmid,pifc
             ubi(*,*,irdg) = ubi_o
             ubm(*,*,irdg) = ubm_o
             nm(*,*,irdg)  = nm_o
             tau_stg0(*,*,irdg) = tau_o

             ; After gw_rdg_below_peak
          readu,1,tauoro,taudsw, hdspwv, hdspdw,effgw,tau_o
             tau_stg1(*,*,irdg) = tau_o

             ; After gw_rdg_break_trap
          if keyword_set(old) then readu,1,tau_o,wbr,tau2oro ;,trn,m2_diag
          if not keyword_set(old) then readu,1,tau_o,wbr
          tau_stg2(*,*,irdg) = tau_o

             ; After gw_drag_prof
          readu,1,tau_o  
             tau_stg3(*,*,irdg) = tau_o
       

             tauf = tau_stg3

        endfor     


         taulin = rsrc * kwvrdgo * ubmsrc * nsrc * (mxdis^2) 

         if not keyword_set(skip_reform) then begin 
         mxdis = reform( mxdis, nlon, nlat, nrdg )
         ubm   = reform( ubm  , nlon, nlat, pver, nrdg )
         nm    = reform( nm   , nlon, nlat, pver, nrdg )
         u     = reform( u    , nlon, nlat, pver )
         v     = reform( v    , nlon, nlat, pver )
         t     = reform( t    , nlon, nlat, pver )
         zm    = reform( zm   , nlon, nlat, pver )
         zi    = reform( zi   , nlon, nlat, pver+1 )
         pmid  = reform( pmid , nlon, nlat, pver )
         pifc  = reform( pifc , nlon, nlat, pver+1 )
         tau_stg2    = reform( tau_stg2   , nlon, nlat, pver+1, nrdg )
         tau_stg3    = reform( tau_stg3   , nlon, nlat, pver+1, nrdg )
         endif
         
    endwhile

    if keyword_set(stops) then STOP

endif




if keyword_set(isotropic) then begin
  while itime lt ntim do begin
        readu,1,nlon,nlat,pver,nrdg,itime,ntim

        taurx0 = dblarr( nlon * nlat, pver+1 )
        taury0 = dblarr( nlon * nlat, pver+1 )
        taurm0 = dblarr( nlon * nlat, pver+1 )

        tau_stg0 = dblarr( nlon * nlat, pver+1 )
        tau_stg1 = dblarr( nlon * nlat, pver+1 )
        tau_stg2 = dblarr( nlon * nlat, pver+1 )
        tau_stg3 = dblarr( nlon * nlat, pver+1 )
        tau_stg4 = dblarr( nlon * nlat, pver+1 )

        tauoro = dblarr( nlon * nlat )
        taudsw = dblarr( nlon * nlat )
        hdspwv = dblarr( nlon * nlat )
        hdspdw = dblarr( nlon * nlat )
        statelat = dblarr( nlon * nlat )

        ubmsrc   = dblarr( nlon * nlat )
        sgh      = dblarr( nlon * nlat )
        sgh_sc   = dblarr( nlon * nlat )
        effgw   = dblarr( nlon * nlat )
        effgw_x = dblarr( nlon * nlat )
        lons    = dblarr( nlon )
        lats    = dblarr( nlat )

        src_level = lonarr( nlon * nlat )
      
        
        v = dblarr( nlon * nlat, pver )
        u = dblarr( nlon * nlat, pver )
        t = dblarr( nlon * nlat, pver )
        nm   = dblarr( nlon * nlat, pver )
        zm   = dblarr( nlon * nlat, pver )
        ubm  = dblarr( nlon * nlat, pver )
        ubi  = dblarr( nlon * nlat, pver+1 )
        pifc = dblarr( nlon * nlat, pver+1 )
        pmid = dblarr( nlon * nlat, pver )
        zi   = dblarr( nlon * nlat, pver+1 )
        utgw = dblarr( nlon * nlat, pver )
        vtgw = dblarr( nlon * nlat, pver )
        gwut = dblarr( nlon * nlat, pver )


        readu,1,lons,lats,statelat
        readu,1,sgh,sgh_sc
        readu,1,src_level
        readu,1,tau_stg0,ubi,ubm,nm,zm,zi,u,v,t,pmid,pifc

        readu,1,tau_stg3  
     
          
         sgh   = reform( sgh  , nlon, nlat )
         sgh_sc= reform( sgh_sc, nlon, nlat )
         ubm   = reform( ubm  , nlon, nlat, pver )
         u     = reform( u    , nlon, nlat, pver )
         v     = reform( v    , nlon, nlat, pver )
         t     = reform( t    , nlon, nlat, pver )
         nm    = reform( nm   , nlon, nlat, pver )
         zm    = reform( zm   , nlon, nlat, pver )
         zi    = reform( zi   , nlon, nlat, pver+1 )
         pmid  = reform( pmid , nlon, nlat, pver )
         pifc  = reform( pifc , nlon, nlat, pver+1 )
         tau_stg3    = reform( tau_stg3   , nlon, nlat, pver+1 )
         tau_stg0    = reform( tau_stg0   , nlon, nlat, pver+1 )
    endwhile
endif

if keyword_set(endstop) then STOP
end
