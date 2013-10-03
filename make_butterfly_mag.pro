;----------------------------------------------------------------------------->

;Process Magnetogram images before remapping to Lat-Lon grid
function make_butterfly_mag, infile, index=outind, dofilter=dofilter

mreadfits,infile,ind,dat

if keyword_set(dofilter) then begin
	dat=filter_image(dat,/MEDIAN)
endif

index2map,ind,dat,map

;map=rot_map(map,-map.roll_angle)
img=rot(dat,-map.roll_angle)
map.roll_angle=0.
map.data=dat

outind=ind
outmap=map

return, outmap

end

;----------------------------------------------------------------------------->