pro polarxy, type=type, loc=loc, start=start, finish=finish, r0=r0 $
             , ct=ct,plotrange0=plotrange0, log=log, xrange=xrange, nopert=nopert, basic=basic $
             ,xtickinterval=xtickinterval,mp=mp,yrange=yrange, scale=scale, nonaxi=nonaxi, cart=cart
if not keyword_set(ct) then ct=5
if not keyword_set(finish) then finish=start
common consts, pi, nrad, time
!p.font = 0
pi=!dpi

;;;;;;;;;;;;;;;
;WHERE IS DATA;
;;;;;;;;;;;;;;;
location=strcompress(loc,/remove_all)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;GET DIMENSIONS AND TIME UNIT;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
dims=(read_ascii(filepath('dims.dat',root_dir='.',subdir=[location]))).(0)
nout=fix(dims(5))
nrad=fix(dims(6))
nsec=fix(dims(7))
nlines = file_lines(filepath('planet0.dat',root_dir='.',subdir=[location]))
info=dblarr(11,nlines)
openr,3,filepath('planet0.dat',root_dir='.',subdir=[location])
readf,3,info
close,3
if not keyword_set(r0) then begin 
   a0=info(1,0)
endif else a0 = r0 
dt=info(7,1)
p0=2.*pi*(a0)^(3./2.)


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;CREATE ARRAY FOR AZIMUTH AND RADIAL;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
azi=dblarr(nsec)
azi1=dblarr(nsec)
radtmp=dblarr(nrad+1)
rad=dblarr(nrad)
openr,1,filepath('used_rad.dat',root_dir='.',subdir=[location])
readf,1,radtmp
close,1
for i=0, nsec-1 do azi(i)=2.*pi*i/nsec
rad(0:nrad-1)=(radtmp(0:nrad-1)+radtmp(1:nrad))/2.

if keyword_set(cart) then begin
   dazi = azi(1) - azi(0)
   dlogr = alog(radtmp(nrad)/radtmp(0))/nrad

   dnx1 = 1024
   xaxis = -rad(nrad-1) + 2*rad(nrad-1)*dindgen(dnx1)/(dnx1-1.0)
   yaxis = xaxis
   dataxy = dblarr(dnx1, dnx1)
endif


if keyword_set(mp) then begin
   f0 = (mp/3d0)^(1d0/3d0)
   xtitle = textoidl('(r-r_p)/r_h')

   azi1 = azi/!dpi - 1d0
   ytitle = textoidl('(\phi - \phi_p)/2\pi')
endif else begin
   xtitle = textoidl('r/r_0')
   
   azi1   = azi/(2d0*!dpi)
   ytitle = textoidl('\phi/2\pi')
endelse

;;;;;;;;;;;;;;;;;;;;;;;;
;DO POLAR CONTOUR PLOTS;
;;;;;;;;;;;;;;;;;;;;;;;;
data = dblarr(nsec,nrad)
data0=dblarr(nsec,nrad)
dataplot = dblarr(nrad,nsec)
loadct,ct, bottom=0

if not keyword_set(nopert) then begin
   if not keyword_set(basic) then begin 
      basic = 0
   endif else basic = 0
   
   openr,2,filepath(strcompress('gas'+type+string(basic)+'.dat',/remove_all),root_dir='.',subdir=[location])
   readu,2,data0
   close,2
   
   dataT = transpose(data0)
   
   dataplot0 = dataT
endif

for k=start, finish do begin
   ks=string(k,format='(I03)')
   
   openr,2,filepath(strcompress('gas'+type+string(k)+'.dat',/remove_all),root_dir='.',subdir=[location]) 
   readu,2,data
   close,2
   if keyword_set(scale) then data *= scale
   dataT = transpose(data)


   plx=info(1,k)
   ply=info(2,k)
   plrad=sqrt(plx*plx+ply*ply)
   phi=pltphi(plx,ply)
   
   if keyword_set(mp) then begin
      rhill = f0*plrad
      rplot = (rad - plrad)/rhill
   endif else begin
      rplot=rad/r0
   endelse
   

   temp=min(abs(azi-phi),grid)   
   if(phi gt !dpi) then begin
      temp2=min(abs(azi-(phi-!dpi)),grid2)
      dataplot(0:nrad-1,0:nsec-1-grid2) = dataT(0:nrad-1,grid2:nsec-1)
      dataplot(0:nrad-1,nsec-grid2:nsec-1) = dataT(0:nrad-1, 0:grid2-1)
   endif
   
   if(phi lt !dpi) then begin
      temp2=min(abs(azi-(phi+!dpi)),grid2)
      dataplot(0:nrad-1,nsec-grid2:nsec-1) = dataT(0:nrad-1,0:grid2-1)
      dataplot(0:nrad-1,0:nsec-1-grid2) = dataT(0:nrad-1, grid2:nsec-1)
   endif
   
   if(phi eq !dpi) then dataplot = dataT
   
   if not keyword_set(nopert) then begin
      dataplot /= dataplot0
      if not keyword_set(log) then  dataplot -= 1.0
   endif
   
   if keyword_set(log) then dataplot=alog10(dataplot)
 
   
   if not keyword_set(plotrange0) then begin 
      temp = min(abs(rplot - xrange(0)),r1)
      temp = min(abs(rplot - xrange(1)),r2)
      plotrange=[min(dataplot(r1:r2,*)),max(dataplot(r1:r2,*))]
   endif else plotrange=plotrange0
   
   levels=plotrange(0)+(plotrange(1)-plotrange(0))*(dindgen(32)/31.)
   levels2= (plotrange(1)/3d0)*(dindgen(2)/1.)
   time=string(info(7,k)/p0,format='(F6.1)')
   title = time+textoidl('P_0, ')+textoidl('r_p='+string(plrad,format='(f4.1)'))
   
   if not keyword_set(cart) then begin
      set_plot, 'ps'
      device, filename=filepath(strcompress('polarxy_'+type+ks+'.ps',/remove_all),root_dir='.',subdir=[location])$
              ,/color, bits_per_pixel=8,xsize=12, ysize=14
      contour,dataplot,rplot,azi1,/fill,levels=levels,title=title, $
              xtitle=xtitle, ytitle=ytitle,charsize=1.5,xmargin=[7.0,6.0], ymargin=[3.5,2], $
              xtickinterval=xtickinterval,xrange=xrange, yrange=yrange, xminor=5, xstyle=1
;contour,dataplot,rplot, azi/pi - 1.0 ,levels=levels2, /overplot
;xyouts, 0.55, 0.02, 'FARGO 2D', charsize=3.5,charthick=8, color=255
      colorbar, position=[0.865, 0.15, 0.90, 0.90],/vertical,/right,range=plotrange,format='(f5.2)'
;oplot,[0.0,0.0],[0.0,0.0]/pi,psym=7,symsize=1.5,color=!D.Table_size*1.0
;oplot,[0.,0.],[0.0,0.0]/pi,psym=7,symsize=1.5,color=!D.Table_size*0.5
;xyouts,plrad,phi/pi,'X',charsize=1.5
      device,/close
   endif else begin
      for jj=0, dnx1-1 do begin
         y = yaxis(jj)
         for ii=0, dnx1-1 do begin
            x = xaxis(ii)
            
            r_t = sqrt(x^2 + y^2)
            azi_t = pltphi(x, y)
            
            if( (r_t ge xrange(0)) and (r_t le rad(nrad-1)) )then begin
               
               temp = min(abs(rad - r_t),   x0)
               temp = min(abs(azi - azi_t),  y0)
               
               dr = rad(x0)*dlogr
               ip = x0 + (r_t - rad(x0))/(dr/2.0) + 0.0
               jp = y0 + (azi_t - azi(y0))/(dazi/2.0) + 0.0
               
               dataxy(ii,jj) = bilinear(dataplot, ip, jp)
            endif else begin
               dataxy(ii,jj) = 1d10
            endelse
         endfor
      endfor
      
      set_plot, 'ps'
      device, filename=filepath(strcompress('polarxy2_'+type+ks+'.ps',/remove_all),root_dir='.',subdir=[location])$
              ,/color, bits_per_pixel=8,xsize=14, ysize=12
      contour, dataxy, xaxis, yaxis, /isotropic,/fill,levels=levels,title=title $
               ,ymargin=[2,2],xmargin=[4,4], charsize=1., xtickinterval=xtickinterval, ytickinterval=xtickinterval
      colorbar, position=[0.85, 0.09, 0.9, 0.91],/vertical,/right,range=plotrange,format='(f5.2)'
;      xyouts, -2.4, -2.4, 'FARGO 2D', charsize=1.5,charthick=8, color=0
      device,/close
                                ;stop
   endelse
   
   print, 'done '+ks
endfor
end
