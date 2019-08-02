external EX01 "/DFS-B/DATA/pritchard/mdfowler/GPI_daily/bootstrap_ERAI/pcmin_2013_fixed.so"

begin
  filename_in = getenv ("GPIFILEIN")
  filename_out = getenv ("GPIFILEOUT")
  fin=addfile (filename_in,"r")

; 1) Make potential intensity using Kerry's fortran code:
;        SUBROUTINE PCMIN(SST,PSL,P,T,R,NA,N,PMIN,VMAX,IFL)
 
; Get the SST:

 lat = fin->lat
 lon = fin->lon
; time = fin->day

 ny = dimsizes(lat)
 nx = dimsizes(lon)
 nt = 365
 nz = dimsizes(fin->levT)
 
; dummy allocation:
 pmin = 1.
 vmaxaux = 1.
 ifl = 0 

 VMAX1 = fin->SST ; put in time x lat x lon
 VMAX = VMAX1(day|:,lat|:,lon|:)
 VMAX@long_name = "Potential intensity"
 VMAX@units = "m/s"
 VMAX@_FillValue = -32767
 VMAX(:,:,:) = VMAX@_FillValue

 TS = fin->SST
 TT = fin->T
 Q  = fin->Q
 P  = fin->levT 
 PS = fin->SLP

; Re-order so that dimensions are in order used in NCL functions (time x lev x lat x lon)
 TSre = TS(day|:,lat|:,lon|:)
 TTre = TT(day|:,levT|:,lat|:,lon|:)
 Qre = Q(day|:,levT|:,lat|:,lon|:)
 PSre = PS(day|:,lat|:,lon|:)

; Loop over lat/lon/time and compute GPI at each point
 do kx=0,nx-1
  do ky=0,ny-1
    if (.not.ismissing(TSre(1,ky,kx)) ) then ; is ocean?
; actually do stuff if ocean grid cell.
     do kt = 0,nt-1
       SST = doubletofloat(TSre(kt,ky,kx)-273.15)  ;[deg C]
       PSL = doubletofloat(PSre(kt,ky,kx)/1.e2)    ;[mb]
       Pmb = doubletofloat(P(::-1))   ;[mb]
       ; Pmb = doubletofloat(P/1.e2)     ;[mb]
       T = doubletofloat(TTre(kt,::-1,ky,kx)-273.15)  ;[deg C]
       R = doubletofloat(Qre(kt,::-1,ky,kx)*1000.)    ;[g/kg]
       EX01::PCMIN(SST,PSL,Pmb,T,R,nz,nz,pmin,vmaxaux,ifl)
       if ( ifl .eq. 1 ) then
         VMAX(kt,ky,kx) = vmaxaux ;[m/s]
       end if       
     end do
   ; end if
  end do
 ; print (kx)
 end do
  ;fout=addfile (filename_out,"c")
  ;fout->VMAX = VMAX
; Finished with GPI


;**** 2) Compute absolute vorticity (850 hPa) and relative humidity (700 hPa)
U = fin->U
V = fin->V

U250 = U(0,:,:,:)
V250 = V(0,:,:,:)
U850 = U(1,:,:,:)
V850 = V(1,:,:,:) 

U250re = U250(day|:,lat|:,lon|:)
V250re = V250(day|:,lat|:,lon|:)
U850re = U850(day|:,lat|:,lon|:)
V850re = V850(day|:,lat|:,lon|:)

; Reorder so that dimensions are time x lat x lon when calculating these things
Ppa = P*100 ; Pressure should be in Pascal 
relH = relhum(TTre,Qre,conform(TTre,Ppa,1))
RelHum700 = relH(:,25,:,:)   ;Relative humidity at 700 hPa

; relVort = uv2vr_cfd(U850re, V850re, lat, lon, 2)
relVort   = uv2vrG(U850re(:,::-1,:),V850re(:,::-1,:)) ; uv2vrG requires latitude to be in ascending order

AbsVort1 = fin->SST ; time x lat x lon
AbsVort = AbsVort1(day|:,lat|:,lon|:)
 AbsVort@long_name = "Absolute Vorticity"
 AbsVort@units     = "per sec"
 AbsVort@_FillValue = -32767
 AbsVort(:,:,:) = AbsVort@_FillValue

omega = 7.292e-5
pi = atan(1)*4.
d2r = pi/180.
latRad = lat*d2r
planetVort = doubletofloat(2.*omega*sin(latRad))

do kx=0,nx-1
  do kt=0,nt-1
    AbsVort(kt,:,kx) = relVort(kt,:,kx) + planetVort ;Vorticity at 850 hPa
  end do
end do
print("AbsVort computed")

; **** 3) Compute vertical wind shear 
VSHEAR1 = fin->SST ; time x lat x lon
VSHEAR = VSHEAR1(day|:,lat|:,lon|:)
 VSHEAR@long_name = "Magnitude of vertical wind shear between 850 and 200 hPa"
 VSHEAR@units     = "m/s"
 VSHEAR@_FillValue = -32767
 VSHEAR(:,:,:) = VSHEAR@_FillValue

VSHEAR = sqrt((U250re-U850re)^2. + (V250re-V850re)^2.)

GPI1 = fin->SST ; time x lat x lon
GPI = GPI1(day|:,lat|:,lon|:)
 GPI@long_name = "Genesis Potential Index"
 GPI@units     = "[none]"
 GPI@_FillValue = -32767
 GPI(:,:,:) = GPI@_FillValue

; *** 4) Compute GPI as the product of individual terms and save output to file
term1 = abs((10.^5) * AbsVort)^(3./2.)
term2 = (RelHum700/50.)^3.
term3 = (VMAX/70.)^3.
term4 = (1. + (0.1*VSHEAR))^-2.

; printVarSummary(term1)
; printVarSummary(term2)
; printVarSummary(term3)
; printVarSummary(term4) 

GPI = term1*term2*term3*term4

fout=addfile (filename_out,"c")
filedimdef(fout,"time",-1,True)
fout->VMAX = VMAX
fout->ABSVORT = AbsVort
fout->VSHEAR = VSHEAR
fout->RELHUM = RelHum700
fout->GPI = GPI

end
