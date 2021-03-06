; -----------------------------------------------------------
; NAME:
;       plot_pre_post_ak
;
; PURPOSE:
;       plot difference between raw GEOS-Chem output and output
;       convolved with averaging kernels
;
; CALLING SEQUENCE:
;       plot_pre_post_ak, diag=diag, trac=trac
;
; INPUTS:
;       DIAG: Diagnostic to be compared
;       TRAC: Tracer to be compare
;
; OUTPUT:
;       plot of GEOS_Chem output vs output convolved with averaging kernels
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
;       plot_pre_post_ak, diag='IJ-AVG-$',trac=4
;
; CREATED:
;       30.06.15 by Kaitlyn Lieschke
; -----------------------------------------------------------

pro plot_pre_post_ak, diag=diag, trac=trac

; Select fts file and path
if trac eq 4 then begin
 ftsfile = '/short/m19/kjl574/sfit/fts_co_woll_2008.hdf'
 spec='co'
 capspec='CO'
endif
if trac eq 20 then begin
 ftsfile = '/short/m19/kjl574/sfit/hcho_sfit2.dat'
 spec='hcho'
 capspec='HCHO'
endif
if trac eq 68 then begin
 ftsfile = '/short/m19/kjl574/sfit/fts_hcn_woll_2008'
 spec='hcn'
 capspec='HCN'
endif

; READ IN FTS DATA
; Open file and assign it a file ID
fileid = hdf_sd_start(ftsfile, /read)
; Find the number of file attributes and datasets
hdf_sd_fileinfo, fileid, numsds, numatts
; Identify the dataset to save from the file, select it and assign it a dataset ID
dataset_name = 'DATETIME'
dataset_index = hdf_sd_nametoindex(fileid, dataset_name)
datasetid = hdf_sd_select(fileid, dataset_index)
; Import the setected dataset
hdf_sd_getdata, datasetid, data_date

; Repeat for next dataset (Pressure centres)
dataset_name = 'PRESSURE_INDEPENDENT'
dataset_index = hdf_sd_nametoindex(fileid, dataset_name)
datasetid = hdf_sd_select(fileid, dataset_index)
hdf_sd_getdata, datasetid, data_presc

; Repeat for next dataset (Altitude centres)
dataset_name = 'ALTITUDE'
dataset_index = hdf_sd_nametoindex(fileid, dataset_name)
datasetid = hdf_sd_select(fileid, dataset_index)
hdf_sd_getdata, datasetid, data_altc

; Repeat for next dataset (Altitude boundaries)
dataset_name = 'ALTITUDE.BOUNDARIES'
dataset_index = hdf_sd_nametoindex(fileid, dataset_name)
datasetid = hdf_sd_select(fileid, dataset_index)
hdf_sd_getdata, datasetid, data_altb

; Repeat for next dataset (CO ppbv)
dataset_name = capspec+'.COLUMN_ABSORPTION.SOLAR'
dataset_index = hdf_sd_nametoindex(fileid, dataset_name)
datasetid = hdf_sd_select(fileid, dataset_index)
hdf_sd_getdata, datasetid, data_co

; Repeat for next dataset (CO avg kernels)
dataset_name = capspec+'.COLUMN_ABSORPTION.SOLAR_AVK'
dataset_index = hdf_sd_nametoindex(fileid, dataset_name)
datasetid = hdf_sd_select(fileid, dataset_index)
hdf_sd_getdata, datasetid, data_avk

; Repeat for next dataset (Systematic error)
dataset_name = capspec+'.COLUMN_ABSORPTION.SOLAR_UNCERTAINTY.SYSTEMATIC.STANDARD'
dataset_index = hdf_sd_nametoindex(fileid, dataset_name)
datasetid = hdf_sd_select(fileid, dataset_index)
hdf_sd_getdata, datasetid, data_syserr

; Repeat for nest dataset (Random error)
dataset_name = capspec+'.COLUMN_ABSORPTION.SOLAR_UNCERTAINTY.RANDOM.STANDARD'
dataset_index = hdf_sd_nametoindex(fileid, dataset_name)
datasetid = hdf_sd_select(fileid, dataset_index)
hdf_sd_getdata, datasetid, data_ranerr

; Repeat for nest dataset (Random error)
dataset_name = capspec+'.COLUMN.PARTIAL_ABSORPTION.SOLAR_APRIORI'
dataset_index = hdf_sd_nametoindex(fileid, dataset_name)
datasetid = hdf_sd_select(fileid, dataset_index)
hdf_sd_getdata, datasetid, data_apri

; Altitude boundaries are stored in two vectors, tops and bottoms
; Reorganise data_altb to be a single vector
data_alt = fltarr(49)
data_alt[0:47] = data_altb[*,0]
data_alt[48] = data_altb[47,1]
; Create variable to store pressure boundaries
data_pres = fltarr(49,n_elements(data_date))

; Convert altitude boundaries at each grid box to pressure boundaries
for z=0,n_elements(data_presc)/48-1 do begin
log_pres = interpol(alog(data_presc[*,z]),data_altc,data_alt,/spline)
data_pres[*,z] = EXP(log_pres)
endfor

; Round to the nearest day and average all values within the same day
data_date = ceil(data_date)
data_dt = tapply(data_date,data_date,'mean',/nan)
fts_final = tapply(data_co,data_date,'mean',/nan)
sys_final = tapply(data_syserr,data_date,'mean',/nan)
ran_final = tapply(data_ranerr,data_date,'mean',/nan)

data_ak = fltarr(48,n_elements(data_dt))
data_ap = fltarr(48,n_elements(data_dt))
data_pr = fltarr(49,n_elements(data_dt))
for ll = 0,47 do begin
data_ak[ll,*] = tapply(reform(data_avk[ll,*]),data_date,'mean',/nan)
data_ap[ll,*] = tapply(reform(data_apri[ll,*]),data_date,'mean',/nan)
endfor
for ll = 0,48 do begin
data_pr[ll,*] = tapply(reform(data_pres[ll,*]),data_date,'mean',/nan)
endfor

; Convert mjd2000 date/time value to tau0 date/time to match GEOS-Chem
; Note time is UTC time
data_tau = data_dt*24+131472
; Convert tau values to yymmdd hhmmss for later use
date_final = tau2yymmdd(data_tau,/nformat)

; Create a variable to hold final model and measurement data and errors
gc_final = fltarr(n_elements(data_tau))
; Create variable to mark movement of data points through months
init = 0

; READ IN MODEL DATA
;===============================================================
; Create variable to count fts measurements / model data read in
count = 0

; Set filename and path for model data
gcfile = '/short/m19/kjl574/kjl-honours-run/2008_run/stations.20080102'

; Determine first time point of model output
mytau = nymd2tau(20080102, 010000L)

; Determine number of time points within the month and create arrays
endtau = nymd2tau(20080202, 010000L)
nhrs = endtau-mytau
ndays = (endtau-mytau)/24

hrarr = lonarr(nhrs)
dyarr = lonarr(ndays)

; Read model file but only save data for dates matching fts measurements
; Create variable to mark passing through the fts files
t = 0
; Create variables to store data
wdata = fltarr(47)
wpress = fltarr(47)

; If first fts value is earlier than first model value, go to next fts value
if (count eq 0 and data_tau[0] lt mytau) then begin
data_tau = data_tau[1:*]
date_final = date_final[1:*]
fts_final = fts_final[1:*]
sys_final = sys_final[1:*]
ran_final = ran_final[1:*]
gc_final = gc_final[1:*]
data_ak = data_ak[*,1:*]
data_ap = data_ap[*,1:*]
data_pr = data_pr[*,1:*]
endif

; Loop through each hour of the model data
for h=0,n_elements(hrarr)-1 do begin
; If all tau values for the month ahve been matched then break
if t eq n_elements(data_tau) then break
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
endif
mytau = mytau + 1
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
`endfor

; Allocate fts pressure levels for model month to new variable
fts_pres = data_pres[*,0:nsets-1]
fts_ak = data_avk[*,0:nsets-1]
fts_apri = data_apri[*,0:nsets-1]
; Reverse fts pressure levels to match model output
fts_pres_rev = reverse(fts_pres,1)
fts_ak_rev = reverse(fts_ak,1)
fts_apri_rev = reverse(fts_apri,1)

; Create variable to store regridded data
gc_regrid = fltarr(48,nsets)
; need ro regrid model output, using relevel.pro from Jesse Greenslade
;print,'Pre-regrid: ',gc_molecs[*,0];==================================================
for z=0,nsets-1 do begin
gc_tmp = relevel(gc_pres[*,z],fts_pres_rev[*,z],gc_molecs[*,z],/keepvmr)
gc_regrid[*,z] = gc_tmp
endfor
print,'Post-regrid: ',gc_regrid[*,0];=================================================
; Create variable to store data convolved with averaging kernels
gc_post = fltarr(48,nsets)
; Make averaging kernel vector into a square matrix
for z=0,nsets-1 do begin
mat_ak = diag_matrix(fts_ak_rev[*,z])
; Apply averaging kernels to model data (AK*GC + APR - AK*APR)
; OR (AK*GC + (I-AK)APR)        OR         (APR + AK(GC-APR)
gc_post[*,z] = mat_ak # gc_regrid[*,z] + fts_apri_rev[*,z] - mat_ak # fts_apri_rev[*,z]
endfor

print,'Post-ak: ',gc_post[*,0];======================================================

; Plot model values before and after averaging kernels were applied
plot,reverse(gc_post[*,0]),fts_pres[*,0],thick=2,color=1,yrange=[1000,0],$
title='GC '+spec+' pre and post AKs',ytitle='pressure',xtitle=spec+' abundance',$
xrange=[-0.5d13,2.5d14]
oplot,fts_apri[*,0],fts_pres[*,0],thick=2,color=!myct.blue
oplot,reverse(gc_regrid[*,0]),fts_pres[*,0],thick=2,color=!myct.red
stop

; Plot averaging kernels (log pressure)
;plot, fts_ak_rev[*,0],alog(fts_pres_rev[*,0]),thick=2,color=1,title=capspec+' Averaging Kernels',$
;ytitle='Log Pressure'
; Plot averaging kernels
;plot, reverse(fts_ak_rev[*,0]),fts_pres[*,0],thick=2,color=1,yrange=[1000,0],$
;title=capspec+' Averaging Kernels',ytitle='Pressure',xrange=[0.8,1.1]
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; Plot:
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; fts_ak_rev[*,0]
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; gc_regrid[*,0]
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; gc_post[*,0]

; Total measurements and model to give total column value
gc_final[init:(init+nsets-1)] = total(gc_post,1)

; Reset init and month values for following month
init = init + nsets
endif

; Write date, model+fts data and random+systematic errors into text file


stop
end
