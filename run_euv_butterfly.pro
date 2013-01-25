pro run_euv_butterfly, res1tore=res1tore, res2tore=res2tore

datapath='/Volumes/My Book/Larisza/EIT/lz/'
savpath='~/science/papers/active_regions_3_diffusion/sav/'

ffff=''

yyyy=['1996','1997','1998']
mm=string(indgen(12)+1,form='(I02)')

if not keyword_set(res1tore) then begin
	for i=0,n_elements(yyyy)-1 do begin
		for j=0,n_elements(mm)-1 do begin
			ff=file_search(datapath+yyyy[i]+'/'+mm[j]+'/*')
			if ff[0] eq '' then continue
			ffu=ff[uniq(strmid(ff,43,8))]
			ffff=[ffff,ffu]
		endfor
		print,yyyy[i]
	endfor
	ffff=ffff[1:*]
	
	save,ffff,file='euv_butterfly_fits_filename.sav'
endif else restore,'euv_butterfly_fits_filename.sav',/ver

stop

make_butterfly_map, ffff, /euv, path=savpath, res1tore=res2tore

avgint=total(MEANINT60[*,80:99],2)/20.
avgint[where(NIMGARR60 eq 0)]=1.
avgimg=transpose(congrid(transpose(avgint),181,1097))

stop

end