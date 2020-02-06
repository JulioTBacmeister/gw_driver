function bump_time_vec,ivc,dtvc=dtvc,dt$=dt,trim_hhmm=trim_hhmm,otg_trim=otg_trim,hhonly=hhonly,sssss=sssss

; PURPOSE:
;     Bumps time vector of the form 
;         [ year ,month , day , hour, min ]
;     or
;         [ year ,month , day , hour ]  ; if keyword_set(hhonly) 
;     or
;         [ year ,month , day , seconds ]  ; if keyword_set(sssss) 


imn=-99
ihr=-99
isc=-99L

if not keyword_set(hhonly) and not keyword_set(sssss)  then begin
   imn = ivc(4)
endif else begin
   imn = -99
endelse
  


if not keyword_set(sssss) then ihr = ivc(3)
if keyword_set(sssss) then isc = ivc(3)
idy = ivc(2)
imo = ivc(1)
iyr = ivc(0)

if keyword_set(dt) then begin
case dt of
'30min': imn=imn+30
'1hr':   ihr=ihr+1
'2hr':   ihr=ihr+2
'3hr':   ihr=ihr+3
'4hr':   ihr=ihr+4
'5hr':   ihr=ihr+5
'6hr':   ihr=ihr+6
'1dy':   idy=idy+1
'5dy':   idy=idy+5

'-1hr':  ihr=ihr-1
'-2hr':  ihr=ihr-2
'-3hr':  ihr=ihr-3
'-4hr':  ihr=ihr-4
'-5hr':  ihr=ihr-5
'-6hr':  ihr=ihr-6
endcase
endif

if keyword_set(dtvc) then begin
   if not keyword_set(hhonly) and not keyword_set(sssss) then imn=imn+dtvc(4)
   if not keyword_set(sssss) then ihr = ihr+dtvc(3)
   idy = idy+dtvc(2)
   imo = imo+dtvc(1)
   iyr = iyr+dtvc(0)
   if keyword_set(sssss) then isc = isc + dtvc(3)
endif   



if isc ge 86400L then begin
   idy = idy+1
   isc = isc - 86400L
endif

if imn ge 60 then begin
   ihr=ihr+1
   imn=imn-60
endif

if ihr ge 24 then begin
   idy=idy+1
   ihr=ihr-24
endif

if idy gt days_in_month(imo,iyr) then begin
   idy=idy-days_in_month(imo,iyr)
   imo=imo+1
 endif

if imo gt 12 then begin
   iyr=iyr+1
   imo=1
endif

if ihr lt 0 and not keyword_set(sssss) then begin
   idy=idy-1
   ihr=ihr+24
endif

if idy le 0 then begin
   imo=imo-1
   if imo le 0 then begin
      imo=imo+12 & iyr=iyr-1
   endif
   idy=idy+ days_in_month(imo,iyr)
endif

iyr=fix(iyr)
imo=fix(imo)
idy=fix(idy)
ihr=fix(ihr)
imn=fix(imn)

if not keyword_set(hhonly) and not keyword_set(sssss)  then begin
   ovc = [  iyr , imo , idy, ihr, imn ]
endif
if keyword_set(hhonly) then begin 
   ovc = [  iyr , imo , idy, ihr ]
endif
  
if keyword_set(sssss) then begin 
   ovc = [  iyr , imo , idy, isc ]
endif




return,ovc
end
