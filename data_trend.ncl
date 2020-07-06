;*************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/time_axis_labels.ncl"
;************************************************
begin
;************************************************
; read in netCDF file
;************************************************
 dir="/Users/irenecionni/"
 a = addfile(dir+"/obs_c3s/Tier3/CDS-GTO-ECV/OBS_ESACCI-OZONE_sat_L3_Amon_toz_199601-201612.nc" ,"r")
;************************************************
; read in zonal winds
;************************************************
  tropoz=a->toz                                
  dim=dimsizes(tropoz)
  opt = False
    avail1=new((/dim(1),dim(2)/),typeof(tropoz))
     do ilat=0,dim(1)-1
      do ilon=0,dim(2)-1
          avail1(ilat,ilon)=num(.not.ismissing(tropoz(:,ilat,ilon)))/tofloat(dim(0))
      end do
    end do

  tropoz1=rmMonAnnCycTLL(tropoz)   ;remove annual cycle
  trend1  = trend_manken(tropoz1, opt, 0)  ; ===> px(2,nlat,mlon)
  trend1(1,:,:)=(/trend1(1,:,:)*12/) 
  trend1(0,:,:)=where(avail1.gt.0.7,trend1(0,:,:),trend1@_FillValue)
  trend1(1,:,:)=where(avail1.gt.0.7,trend1(1,:,:),trend1@_FillValue)                                   ; if meta data is desired
  copy_VarCoords(tropoz1(0,:,:),trend1(0,:,:))
  trend1!0= "prob_trend"               ; ===> size 2

  wks = gsn_open_wks("png","trend_1_"+"CDS-GTO-ECV")         ; send graphics to PNG file
  gsn_define_colormap(wks,"NEO_div_vegetation_c")
  plot=new(2,"graphic")
  res = True                                    ; plot mods desired
  res@gsnDraw                  = False            ; don't draw yet
  res@gsnFrame                 = False            ; don't advance frame yet

  res@cnFillOn         = True               ; color Fill 
  res@cnFillMode       = "RasterFill"       ; Raster Mode
  res@cnLinesOn        =  False             ; Turn off contour lines
  res@cnLineLabelsOn       =False
  res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
  res@cnMinLevelValF       = -1               ; set min contour level
  res@cnMaxLevelValF       = 1              ; set max contour level
  res@cnLevelSpacingF      = 0.2               ; set contour spacing
  ;res@cnLabelBarEndStyle="ExcludeOuterBoxes"
  res@tiMainString         = "Trend "    ; plot title
  fmt = "%3.2f"
  res@gsnLeftString      = "[DU/year]"
  ;res@gsnRightString     = "mean: "+sprintf(fmt,avg_clim1)+" std:"+sprintf(fmt,avg_std1)           ; plot center string
  res@lbOrientation = "vertical"
  ;res@lbLabelBarOn = False
  plot(0) = gsn_csm_contour_map(wks,trend1(1,:,:), res)  ; create plot
  res@cnMinLevelValF       = 0               ; set min contour level
  res@cnMaxLevelValF       = 1               ; set max contour level
  res@cnLevelSpacingF      = 0.05               ; set contour spacing
  res@tiMainString         = "Trend Significance"  
  res@gsnLeftString      = "[1]"
  plot(1) = gsn_csm_contour_map(wks,trend1(0,:,:), res)
  resP=True
  gsn_panel(wks,plot,(/2,1/),resP)
end
