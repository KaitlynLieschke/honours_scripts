; -----------------------------------------------------------
; NAME:
;       PLOT_SFIT_AKS
;
; PURPOSE:
;       compare SFIT2 FTS data to model output, specifically Wollongong vs
;         GEOS-Chem stations output
;       Involves removing unmeasured days/times from model output,
;         applying FTS averaging kernels to model output,
;         plotting both datasets as a time series.
;
; CALLING SEQUENCE:
;       PLOT_SFIT_AKS, diag=diag, trac=trac
;
; INPUTS:
;       DIAG: Diagnostic to be compared
;       TRAC: Tracer to be compared
;
; OUTPUT:
;       csv file containing of GEOS_Chem output vs measured data
;
; SUBROUTINES:
;       ppbv_to_molecs_per_cm2, relevel
;
; NOTES:
;       averaging kernels are stored from top of the atmosphere to the
;       surface (decreasing altitude), model output is stored from
;       surface to the top of the atmosphere (increasing
;       altitude). Measurement pressure levels are reversed for releveling
;       
;       To run with different species change ftsfile name on L48,
;       dataset_name on L81, L87, L93, L99 and L105 and variable name
;       on L84 and L127
;
; EXAMPLE:
;       plot_sfit_aks, diag='IJ-AVG-$',trac=20
;
; CREATED:
;       30.06.15 by Kaitlyn Lieschke
; -----------------------------------------------------------

pro plot_sfit_aks, diag=diag, trac=trac

; Specify fts file and path
ftsfile = '/short/m19/kjl574/sfit/hcho_sfit2.dat'
; Specify apriori file and path
aprifile = '/short/m19/kjl574/sfit/hcho_apriori.dat'
; Specify averaging kernel file and path
akfile = '/short/m19/kjl574/sfit/hcho_ak.dat'

; READ IN FTS DATA
; Open fts file
openr,lun,ftsfile,/get_lun
; Set up variables
hdr=''
tmpdate=0L
data_date=0L
decyr=0.0
tmphcho=0.0
data_hcho=0.0
tmperr=0.0
data_err=0.0
;Read fts file
readf,lun,hdr
while ~eof(lun) do begin
  readf,lun,tmpdate,decyr,tmphcho,tmperr,format='(i8,3x,f9.4,3x,f9.2,5x,f9.2)'
  data_date=[data_date,tmpdate]
  data_hcho=[data_hcho,tmphcho]
  data_err=[data_err,tmperr]
endwhile
;Close file and free memory
close,lun
free_lun,lun
; Remove zero values from start of arrays
data_date=data_date[1:*]
data_hcho=data_hcho[1:*]
data_err=data_err[1:*]

; Read in apriori file
; Open fts file
openr,lun,aprifile,/get_lun
; Set up variables
tmpalt=0.0
data_altb=0.0
tmpaltmid=0.0
altmid=0.0
vmrapri=0.0
vmr=0.0
sigvmr=0.0
tmpapri=0.0
data_ap=0.0
col=0.0
;Read fts file
readf,lun,hdr
while ~eof(lun) do begin
  readf,lun,tmpalt,tmpaltmid,vmrapri,vmr,sigvmr,tmpapri,col,$
        format='(2x,f5.2,2x,f5.2,3x,f10.4,3x,f10.4,3x,f10.4,3x,f10.4,3x,f10.4)'
  data_altb=[data_altb,tmpalt]
  altmid=[altmid,tmpaltmid]
  data_ap=[data_ap,tmpapri]
endwhile
;Close file and free memory
close,lun
free_lun,lun
; Remove zero values from start of arrays
data_altb=data_altb[1:*]
altmid=altmid[1:*]
data_ap=data_ap[1:*]

; Read in averaging kernel file
; Open fts file
openr,lun,akfile,/get_lun
; Set up variables
alt=0.0
tmpak=0.0
data_ak=0.0
tmpaksd=0.0
ak_sd=0.0
;Read fts file (no header on this file)
while ~eof(lun) do begin
  readf,lun,alt,tmpak,tmpaksd,$
        format='(3x,f5.2,2x,f10.3,2x,f10.3)'
  data_ak=[data_ak,tmpak]
  ak_sd=[ak_sd,tmpaksd]
endwhile
;Close file and free memory
close,lun
free_lun,lun
; Remove zero values at start of arrays
data_ak=data_ak[1:*]

; Re-name fts data for later use
fts_final = data_hcho

; Apriori data is accumulated from top down so needs to be separated:
data_apri = fltarr(n_elements(data_ap))
data_apri[0] = data_ap[0]
for l = 1,43 do begin
   data_apri[l] = data_ap[l]-data_ap[l-1]
endfor

; Altitude boundaries are stored in two vectors, mids and bottoms
; Reorganise data_altb to be a single vector
data_alt = fltarr(45)
data_alt[1:44] = data_altb
data_alt[0] = altmid[0]+(altmid[0]-data_altb[0])
; Create variable to store pressure boundaries
data_pres = fltarr(45)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;FIX PRESS
; Convert altitude boundaries at each grid box to pressure boundaries
; Establish values for use in calculations
psurf = 1013.25
H = (8.314*240)/(0.0289751*9.8) ;---scale height
data_pres=fltarr(n_elements(data_alt))
; Calculate pressure edges
data_pres = psurf*EXP(-(data_alt*1000)/H)

; Convert date to tau0 date/time to match GEOS-Chem
data_tau = nymd2tau(long(data_date), 000000L)
; Save yymmdd values for later use
date_final = data_date

; Create a variable to hold original and final model data
gc_final = fltarr(n_elements(data_tau))
gc_orig = fltarr(n_elements(data_tau))
gc_rebin = fltarr(n_elements(data_tau))

; Create variable to mark movement of data points through months
init = 0

; READ IN MODEL DATA
;===============================================================
; Create variable to count fts measurements / model data read in for
; entire year (not reset at start of each month)
count = 0

; Loop through model files (one file per month)
for month = 0,11 do begin
   mm = month+1
   ; If all data_tau values have been matched then break loop
   if count eq n_elements(data_date) then break
   ; Set filename and path for model data, first convert month to string
   smm = string(mm, format = '(i02)')
   gcfile = '/short/m19/kjl574/kjl-honours-run/2008_run/stations.20080102'

   ; Determine first time point of model output
   mytau = nymd2tau(080102, 010000L)

   ; Determine number of time points within the month and create array of
   ; this size for use in for loops
   if month lt 11 then begin
      endtau = nymd2tau(20080202, 010000L)
      nhrs = endtau-mytau
      ndays = (endtau-mytau)/24
   endif
   if month eq 11 then begin
      endtau = mytau+744
      nhrs = 744
      ndays = 31
   endif
   hrarr = lonarr(nhrs)
   dyarr = lonarr(ndays)

   ; Read model file but only save data for dates matching fts measurements
   ; Create variable to mark passing through the fts files (reset at
   ; start of each month)
   t = 0
   ; Create variables to store data
   wdata = fltarr(47)
   wpress = fltarr(47)

   ; If first fts value is earlier than first model value, go to next fts value
   while (count eq 0 and data_tau[0] lt mytau) do begin
      data_tau = data_tau[1:*]
      date_final = date_final[1:*]
      fts_final = fts_final[1:*]
      data_err = data_err[1:*]
      gc_final = gc_final[1:*]
      gc_orig = gc_orig[1:*]
      gc_rebin = gc_rebin[1:*]
   endwhile

   ; Loop through each hour of the model data
   for h=0,n_elements(hrarr)-1 do begin

   ; Only save model data that matches a measurement date/time
      if mytau eq data_tau[t] then begin
         ; Specify tracerinfo.dat file to read
         ctm_tracerinfo,/force,file='/short/m19/kjl574/kjl-honours-run/2008_run/tracerinfo.dat'
         ; Extract data for specific diagnistic, tracer and time
         ctm_get_data, di, diag, filename=gcfile, tracer=trac, tau0=mytau, status=2
         ; Save Wollongong data into variable (create variable first)
         ; NB di[1,2] represent Darwin and Lauder data
         ctm_read_data, wdata_tmp, di[0], result = r
         ; Extract pressure edges of grid boxes at specific time
         ctm_get_data, pi, 'PEDGE-$', filename=gcfile, tracer=1, tau0=mytau, status=2
         ctm_read_data, wpress_tmp, pi[0], result=r
         ; Reform dimensions of temporary variables from (1,1,47) to (47)
         wdata_tmp = reform(wdata_tmp)
         wpress_tmp = reform(wpress_tmp)
         ; Save data into permanent variables (create variables first)
         wdata = [[wdata],[wdata_tmp]]
         wpress = [[wpress],[wpress_tmp]]
         ; Roll over to next fts tau0 value
         t = t + 1
         count = count + 1
         if n_elements(data_tau) eq t then break
      endif
      mytau = mytau + 1

      ;If all months have been done then break
      if n_elements(data_tau) lt t then begin
         if month eq 11 and data_tau[t] gt endtau then break
      endif

   endfor
   ; Free memory blocked by excessive use of GAMAP package
   ctm_cleanup

   ; If data point exist then continue processing
   if n_elements(wdata) gt 47 then begin
      ; Remove 0 values from beginning of arrays
      wdata = wdata[*,1:*]
      wpress = wpress[*,1:*]
      ; Determine number of datasets stored
      nsets = n_elements(wdata)/47

      ; Remove used tau values from data_tau array in preparation for next month
      if nsets lt n_elements(data_tau) then data_tau = data_tau[nsets:*]

      ; Grid box pressures are recorded only at the bottom of each grid box
      ; Add top pressure edge 0.010
      gc_pres = fltarr(48,nsets)
      gc_pres[0:46,*] = wpress
      gc_pres[47,*] = 0.010

      ; Convert units of wdata from ppbv to molecs/cm2
      gc_molecs = fltarr(47,nsets)
      for z=0,nsets-1 do begin
         temp_pres = reform(gc_pres[*,z])
         temp_data = reform(wdata[*,z])
         wdatat = ppbv_to_molecs_per_cm2(temp_pres,temp_data)
         gc_molecs[*,z] = wdatat
      endfor

      ; Save original model data just in case
      gc_orig[init:(init+nsets-1)] = total(gc_molecs,1)

      ; Allocate fts pressure levels for model month to new variable
      fts_pres = data_pres
      fts_ak = data_ak
      fts_apri = data_apri
      ; Reverse fts pressure levels to match model output
      fts_pres_rev = reverse(fts_pres,1)
      fts_ak_rev = reverse(fts_ak,1)
      fts_apri_rev = reverse(fts_apri,1)

      ; Create variable to store regridded data
      gc_regrid = fltarr(44,nsets)
      ; need ro regrid model output, using relevel.pro from Jesse Greenslade
      ;print,'Pre-regrid: ',gc_molecs[*,0];==================================================
      for z=0,nsets-1 do begin
         gc_tmp = relevel(gc_pres[*,z],fts_pres_rev,gc_molecs[*,z],/keepvmr)
         gc_regrid[*,z] = gc_tmp
      endfor

      ; Save regridded modeled data just in case
      gc_rebin[init:(init+nsets-1)] = total(gc_regrid,1)

      ;print,'Post-regrid: ',gc_regrid[*,0];=================================================
      ; Create variable to store data convolved with averaging kernels
      gc_post = fltarr(44,nsets)
      ; Make averaging kernel vector into a square matrix
      for z=0,nsets-1 do begin
         mat_ak = diag_matrix(fts_ak_rev)
         ; Apply averaging kernels to model data (AK*GC + APR - AK*APR)
         ; OR (AK*GC + (I-AK)APR)        OR         (APR + AK(GC-APR)
         gc_post[*,z] = mat_ak # gc_regrid[*,z] + fts_apri_rev - mat_ak # fts_apri_rev
      endfor

      ;print,'Post-ak: ',gc_post[*,0];======================================================
      ; Total model to give total column value
      gc_final[init:(init+nsets-1)] = total(gc_post,1)
      
; Plot model values before and after averaging kernels were applied
plot, reverse(gc_post[*,0]),fts_pres[*,0],thick=2,color=1,yrange=[1000,0],$
title='GC HCHO pre and post AKs', ytitle='pressure',xtitle='HCHO abundance',$
xrange=[-2d15,6d15]
oplot,fts_apri[*,0],fts_pres[*,0],thick=2,color=!myct.blue
oplot,reverse(gc_regrid[*,0]),fts_pres[*,0],thick=2,color=!myct.red
stop

; Plot averaging kernels
;plot, fts_ak[*,0],fts_pres[*,0],thick=2,color=1,yrange=[1000,0],$
;title='HCHO Averaging Kernels Profile',ytitle='Pressure'
;stop
; Plot averaging kernels (log pressure)
;plot, fts_ak[*,0],alog(fts_pres[*,0]),thick=2,color=1,$;yrange=[1000,0],$
;title='HCHO Averaging Kernels Profile',ytitle='Pressure'
;stop


      ; Reset init and month values for following month
      init = init + nsets
   endif
endfor

; Remove excess values in arrays
date_final=date_final[0:count-1]
gc_final=gc_final[0:count-1]
fts_final=fts_final[0:count-1]
data_err=data_err[0:count-1]
gc_orig = gc_orig[0:count-1]
gc_rebin = gc_rebin[0:count-1]

; Write date, model+fts data and random+systematic errors into text file
; Define speciesv
if trac eq 4 then spec='co'
if trac eq 20 then spec='hcho'
if trac eq 21 then spec='c2h6'
if trac eq 67 then spec='ch3oh'
if trac eq 68 then spec='hcn'

; Create header for csv file
hdr = ['DATE','GC','FTS','ERR','GCORIG','GCREBIN']
; Write csv file containing all necessary data
write_csv,spec+'_'+year+'_gc_fts_data.csv',date_final[0:n_elements(gc_final)-1],gc_final,fts_final,data_err,gc_orig,gc_rebin,header = hdr

stop
end
