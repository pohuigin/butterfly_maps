;----------------------------------------------------------------------------->

;Process EUV images before remapping to Lat-Lon grid
function make_butterfly_euv, infile, logscale=logscale, index=outind, err=err

file=infile

;stop

read_eit,file,index,data
eit_prep,index,data=data,ind,img

;eit_prep,file,head,img
;ind=hdr2struct(head)

;img=img/ind.EXPTIME

if keyword_Set(logscale) then imgscl=alog10(abs(img)+.0001) $
	else imgscl=img

index2map,ind,imgscl,map

;map=rot_map(map,-map.roll_angle)
imgscl=rot(imgscl,-map.roll_angle)
map.roll_angle=0.

map.data=imgscl

outind=ind
outmap=map

return,outmap

end

;----------------------------------------------------------------------------->