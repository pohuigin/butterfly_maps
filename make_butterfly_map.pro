;----------------------------------------------------------------------------->
;+
; PROJECT:  	Physical study of global magnetic field using 
;				SolarMonitor Active Region Tracker (SMART)
;
; PROCEDURE:    make_butterfly_map
;
; PURPOSE:    	gdddrd
;
; USEAGE:     	rdrgdrg
;
; INPUT:        DATE=DATE			- drgdrg
;				FILELIST=FILELIST	- drgdrgdrg
;
; KEYWORDS:   	
;				ERRORS		- drgdrgd
;
; OUTPUT:    
;   	    	drgdrgd
;   	    	
; EXAMPLE:    	drgdrgdr
;             
; AUTHOR:     	6-May-2011 P.A.Higgins - Written
;				29-Feb-2012 P.A.H - Added median filtering option
;
; CONTACT:		pohuigin@gmail.com
; VERSION   	0.1
;-
;----------------------------------------------------------------------------->

;----------------------------------------------------------------------------->

;Re-map images to a Lat-Lon grid
;setting NOREMAP just returns the averaged magnetogram and nothing else.
function make_latlon_img, infiles, mag=setmag, posmag=setpos, negmag=setneg, $
	euv=seteuv, tmid=tmid, logscale=logscale, missval=missval, $
	latbin=latbin, lonbin=lonbin, lonbound=lonbound, latbound=latbound, $
	setrebin=setrebin, rebinfact=rebinfact, dofilter=dofilter, $
	outfluxmap=outfluxmap, outmask=projmask, outhg=hgcoord, outavgimg=avgimg, $
	noremap=noremap

rsunmm=695.5 ;radius of sun in Mm
radpasec=2.*!pi/(360.*3600.)
missval=-9999.

ff=infiles
nff=n_elements(ff)

for j=0,nff-1 do begin

	if keyword_set(seteuv) then map=make_butterfly_euv(ff[j], logscale=logscale, index=outind)

	if keyword_set(setmag) then map=make_butterfly_mag(ff[j], index=outind, dofilter=dofilter)

	if j eq 0 then map_arr=map else map_arr=[map_arr,map]

endfor

inindex=outind

map_rot=drot_map(map_arr,time=anytim(tmid,/vms), miss=missval)

if n_elements(map_rot) lt 2 then imgavg=map_rot.data $
	else imgavg=average(map_rot.data,3,miss=missval)

;Make positive only or negative only version of the magnetogram
if keyword_set(setmag) then begin
	if keyword_set(setpos) then imgavg = imgavg > 0.
	if keyword_set(setneg) then imgavg = imgavg < 0.
endif

;AVERAGE MAP STACK IN Z
;SET MISS VAL, AVERAGE WITHOUT THE MISS VAL


;MAKE GRID ON IMAGE IN HGLAT LON, PROJECT TO LAT-LON GRID, SEE IF IT IS CORRECT!!!
;!!! Need to make this depend on Lat-Lon resolution !!! use input latbin and lonbin
;projdim=[361,181]

nlat=ceil(abs(latbound[1]-latbound[0])/latbin)+1.
nlon=ceil(abs(lonbound[1]-lonbound[0])/lonbin)+1.
projdim=[nlon,nlat]
projlim=[lonbound,latbound]
lon = (findgen(projdim[0])*lonbin+projlim[0]) # replicate(1,projdim[1])
lat = replicate(1,projdim[0]) # (findgen(projdim[1])*latbin+projlim[2])
projmask=fltarr(projdim[0],projdim[1])+1.

;Create World Coordinate System structure
wcs = fitshead2wcs( inindex )
COORD = WCS_GET_COORD(WCS)

;Assume image rotation is 0 degrees
wcs.ROLL_ANGLE=0.

;Get coordinate arrays
xx=reform(coord[0,*,*])
yy=reform(coord[1,*,*])
rr=(xx^(2.)+yy^(2.))^(.5)

;Get radius of Sun in arcsecs, particular to the data source
dsunmm=wcs.position.DSUN_OBS/1d6 ;put dist. to sun in Mm
rsunmm=WCS_RSUN(units='Mm') ;radius of sun in Mm
rsunasec=atan(rsunmm/dsunmm)/radpasec

;stop

if keyword_set(setmag) then begin
	;Create cosine correction map
	coscor=rr
	;coscor=coscor/max(coscor)
	coscor=1./cos(asin(coscor/rsunasec))
	coscor[where(rr gt rsunasec)]=1
	;plot,1./!dtor*asin(rr[512:*,512]/rsunasec),coscor[512:*,512] & vline,60 & hline,2
	;Limit correction to maximum area a pixel can cover on edge of disk
	thetalim=asin(1.-(wcs.cdelt)[0]/(rsunasec)) ;should it be half a pixel or a full pixel? asin(1.-(wcs.cdelt)[0]/(2.*rsunasec))
	coscor=coscor < 1./cos(thetalim)
	
	;Cosine correct the image
	imgavg=imgavg*coscor
	fluxmap=imgavg*(map.dx*map.dy*(rsunmm/map.rsun)^(2.))*coscor^(2.)
endif

avgimg=imgavg
if keyword_set(noremap) then return,avgimg

;!!!TO DO!!! (DONE! 20100526 - P.A.H.)
;output mask of where pixels cover
;output number of images used to make slice
;output transformed image
;output number of pixels included in each latitude slice (collapse masks in longitude??)
; do coordinate conversion here:

if keyword_set(setrebin) then begin
	;Rebin the image to 1/4 to increase SNR
	;Optimized for images of 1k x 1k!! Must be evenly divisable by BINFACTOR
	if n_elements(rebinfact) ne 1 then binfactor=4. else binfactor=rebinfact
	wcs.cdelt=wcs.cdelt*binfactor
	wcs.crpix=wcs.crpix/binfactor
	wcs.naxis=wcs.naxis/binfactor
	;imgavgrb=rebin(imgavg,(wcs.naxis)[0],(wcs.naxis)[1])
	imgavgrb=reduce(imgavg,binfactor,binfactor,/average,miss=missval)
	if keyword_set(setmag) then fluxmaprb=reduce(fluxmap,binfactor,binfactor,/average,miss=missval)
endif else begin
	imgavgrb=imgavg
	fluxmaprb=fluxmap
endelse
undefine,COORD

;Remap image to Lat-Lon grid 
;wcs_convert_to_coord, wcs, coord, ’HG’, lon, lat
HGLN=lon & HGLT=lat
WCS_CONVERT_TO_COORD, WCS, COORD, 'HG', HGLN, HGLT
pixel = wcs_get_pixel( wcs, coord )
proj = reform( interpolate( imgavgrb, pixel[0,*,*], pixel[1,*,*] ))
if keyword_set(setmag) then fluxproj = reform( interpolate( fluxmaprb, pixel[0,*,*], pixel[1,*,*] ))

hgcoord=[[[HGLN]],[[HGLt]]]

;Replace missing pixels
wmiss=where(proj eq missval)
if wmiss[0] ne -1 then begin
	proj[wmiss]=missval
	projmask[wmiss]=0.
	if keyword_set(setmag) then fluxproj[wmiss]=missval
endif


;proj=wcspixcoord( dum, outx, outy, xx=xx, yy=yy, outdex=indproj, index=inindex, data=imgavg, $
;	/project, /hgs, resolution=[1.,1.], /cos_correction)

;OUTPUT MASK OF IMAGE AVERAGING

;OUTPUT N IMAGES FOR SLICE

outproj=proj
if keyword_set(setmag) then outfluxmap=fluxproj

return,outproj

end

;----------------------------------------------------------------------------->

;bin=[lonbin,latbin]
pro make_butterfly_map, files, euv=seteuv, mag=setmag, posmag=setpos, negmag=setneg, $
	path=path, res1tore=res1tore, fnamemod=fnamemod, lonbound=setlonbound, latbound=setlatbound, $
	bin=inbin, dofilter=dofilter

if n_elements(path) ne 1 then path='~/science/data/butterfly2/'

test=1.
save,test,file=path+'test.save'

missval=-9999.

;stop

if not keyword_set(res1tore) then begin
	mreadfits,files,ind_arr,/nodat
	save,ind_arr,file=path+'all_file_index.sav'
endif else restore,path+'all_file_index.sav',/ver

;stop

;ftimes=anytim(ind_arr.date_obs) ;anytim(file2time(files))
;ntim=n_elements(ftimes)
;tmm=minmax(ftimes)
;utimes=ftimes[uniq(ftimes)]
;nuf=n_elements(utimes)
;tmm=minmax(utimes)

;Set parameters for each data type
if keyword_set(setmag) then begin
	;Define parameters
	nfperd=0 ;number of files to use for each time bin
	tbin=3600.*24. ;time bin in seconds
	if n_elements(inbin) eq 2 then begin
		lonbin=inbin[0]
		latbin=inbin[1]
	endif else begin
		lonbin=1. ;longitude bin in degrees
		latbin=1. ;latitude bin in degrees
	endelse
	if n_elements(setlatbound) eq 2 then latbound=setlatbound $
		else latbound=[-90,90] ;minimum and maximum latitudes
	if n_elements(setlonbound) eq 2 then lonbound=setlonbound $
		else lonbound=[-60.,60.]
	;Filter data set
	wgood=where(ind_arr.missvals eq 0 and ind_arr.interval ge 300.)
	wbad=where(ind_arr.missvals eq 0)
	setrebin=0
	rebinfact=1
endif
if keyword_set(seteuv) then begin
	;Define parameters
	nfperd=1. ;number of files to use for each time bin
	tbin=3600.*24. ;time bin in seconds
	latbin=1. ;latitude bin in degrees
	latbound=[-90,90] ;minimum and maximum latitudes
	lonbin=1. ;longitude bin in degrees
	lonbound=[-5.,5.]
	logscale=0;1
	setrebin=0
	rebinfact=1

	;Filter data set of images with missing blocks
	wgood=where(ind_arr.naxis eq 2 and ind_arr.naxis1 ge 512 and ind_arr.naxis2 ge 512)
	wbad=where(ind_arr.naxis ne 2 or ind_arr.naxis1 lt 512 or ind_arr.naxis2 lt 512)
endif

;Save list of "bad" files
if wbad[0] ne -1 then begin
	ffbad=files[wbad]
	ind_arr_bad=ind_arr[wbad]
	ftimesbad=anytim(ind_arr_bad.date_obs)
endif

;Get rid of "bad" files
ffgood=files[wgood]
ind_arr=ind_arr[wgood]

;Find file times
ftimes=anytim(ind_arr.date_obs)
ntim=n_elements(ftimes)
tmm=minmax(ftimes)

;stop

;Time goes from 1-jan of year of first file to year of last file + 1 (incase all files from same year)
tbound=[anytim('1-jan-'+strmid(anytim(tmm[0],/CCSDS),0,4)), anytim('1-jan-'+strtrim(strmid(anytim(tmm[1],/CCSDS),0,4)+1,2))]

nt=ceil((tbound[1]-tbound[0])/tbin)+1.
nlat=ceil((latbound[1]-latbound[0])/latbin)+1.

;Initialize arrays
xarr60=findgen(nt)*tbin+tbound[0]
yarr60=findgen(nlat)*latbin+latbound[0]
NIMGARR60=fltarr(nt)
MASKMAP60=fltarr(nt,nlat)

;Initialize data-type-specific arrays
if keyword_set(setmag) then begin
	;Initialize arrays
	MEANUNSB60=fltarr(nt,nlat)+missval
	MEANNETB60=fltarr(nt,nlat)+missval
	MEANFLUX60=fltarr(nt,nlat)+missval
	exparr60=fltarr(nt)
endif
if keyword_set(seteuv) then begin
	;Initialize arrays
	MEANINT60=fltarr(nt,nlat)+missval
	exparr60=fltarr(nt)
endif

for i=0.,nt-1. do begin
	
	thistbin=[tbound[0]+i*tbin,tbound[0]+(i+1.)*tbin]
	
	;Check for files within the current bin
	wthisf=where(ftimes ge thistbin[0] and ftimes lt thistbin[1])
	if wthisf[0] eq -1 then begin
		if wbad[0] eq -1 then continue
		if keyword_set(seteuv) then continue
		
		;If no files found at all then use "bad" files
		wthisf=where(ftimesbad ge thistbin[0] and ftimesbad lt thistbin[1])	
		if wthisf[0] eq -1 then continue
		
		thisff=ffbad[wthisf]
		thisind=ind_arr_bad[wthisf]
	endif else begin
		thisff=ffgood[wthisf]
		thisind=ind_arr[wthisf]
	endelse
	
	;take the differential rotation time to be half way between the beginning and end of the time bin 
	tmid=mean(thistbin)
	
	;pull out projected image, mask of pixel coverage, number of images used to make slice, and number of pixels included for each latitude slice
	imglatlon=make_latlon_img(thisff, euv=seteuv, mag=setmag, posmag=setpos, negmag=setneg, $
		tmid=tmid, logscale=logscale, missval=missval, $
		latbin=latbin, lonbin=lonbin, latbound=latbound, lonbound=lonbound, $
		outmask=thismask, outfluxmap=thisfluxmap, $
		setrebin=setrebin, rebinfact=rebinfact, dofilter=dofilter)
	
	;stop
	
	;ALREADY CROPPED IN PROJECTION ROUTINE
	;Crop to +-60
	;imglatlon = imglatlon[120:240,*]
	
	NIMGARR60[i]=n_elements(thisff)
	MASKMAP60[i,*]=total(thismask,1)
	
	;Initialize data-type-specific arrays
	if keyword_set(setmag) then begin
		;Fill arrays
		MEANUNSB60[i,*]=average(abs(imglatlon),1,miss=missval)
		MEANNETB60[i,*]=average(imglatlon,1,miss=missval)
		MEANFLUX60[i,*]=average(thisfluxmap,1,miss=missval)
		exparr60[i]=mean(thisind.interval)
		
	endif
	if keyword_set(seteuv) then begin
		;Fill arrays
		MEANINT60[i,*]=average(imglatlon,1,miss=missval)
		exparr60[i]=mean(thisind.EXPTIME)
	endif

if keyword_set(seteuv) then begin
	if logscale then begin
		if i mod 10 eq 0 then plot_image,MEANINT60[0:i,*]>(0)<5,/nosq
	endif else begin
		if i mod 10 eq 0 then plot_image,MEANINT60[0:i,*]>(0)<100,/nosq
	endelse
endif
if keyword_set(setmag) then $
	if i mod 100 eq 0 then plot_image,MEANNETB60[0:i,*]>(-10)<10,/nosq

;AVERAGE PROJECTED IMAGE IN LONGITUDE
;INSERT INTO BUTTERFLY MAPS

;TEMP!!!!!!!!
;plot_image,MEANNETB60 > (-10.) < (10.),/nosq,/noerase


endfor

plot_image,MEANNETB60 > (-10.) < (10.),/nosq,/noerase

;stop

if keyword_set(setmag) then begin
	info={tbin:tbin,tminmax:tmm,tbound:tbound,latbin:latbin,lonbin:lonbin,latbound:latbound,lonbound:lonbound,missval:missval}
	pol='net'
	if keyword_set(setpos) then pol='pos'
	if keyword_set(setneg) then pol='neg'
	if n_elements(fnamemod) eq 1 then fmod='_'+strtrim(fnamemod,2) else fmod=''

	magfname=path+'butterfly_mag_'+time2file(ind_arr[0].date_obs,/date)+'_'+time2file((reverse(ind_arr.date_obs))[0],/date)+'_lon'+strjoin(strtrim(fix(lonbound),2),'')+'_pol'+pol+fmod+'.sav'
	;magfname=path+'butterfly_mag_'+time2file(systim(/utc))+'_lon'+strjoin(strtrim(fix(lonbound),2),'')+'_pol'+pol+fmod+'.sav'
	print,magfname
	save,xarr60,yarr60,NIMGARR60,MASKMAP60,exparr60,MEANUNSB60,MEANNETB60,file=magfname
	;save,xarr60,yarr60,NIMGARR60,MASKMAP60,exparr60,MEANUNSB60,MEANNETB60,file=path+'butterfly_mag_'+time2file(systim(/utc))+'.sav',/comp
endif
if keyword_set(seteuv) then begin
	info={tbin:tbin,tminmax:tmm,tbound:tbound,latbin:latbin,lonbin:lonbin,latbound:latbound,lonbound:lonbound,missval:missval}
	save,xarr60,yarr60,NIMGARR60,MASKMAP60,exparr60,MEANINT60,file=path+'butterfly_euv_'+time2file(systim(/utc))+'.sav',/comp

endif



;stop

end
