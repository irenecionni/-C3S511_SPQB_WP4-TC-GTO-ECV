;*************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/time_axis_labels.ncl"
;************************************************
;functions
; #############################################################################
undef("coswgt_area_average")
function coswgt_area_average(field:numeric)
;
; Arguments
;    field: numeric field.
;
; Return value
;    The area average using cosine lat weights.
;
; Description
;    Computes the area average using cosine lat weights and lon weights=1.
;
local funcname, scriptname, lat, wgt_lat, lon, lon_size, wgt_lon, ave
begin

  funcname = "coswgt_area_average"
  scriptname = "diag_scripts/shared/statistics.ncl"

  lat = field&lat
  wgt_lat = tofloat(NormCosWgtGlobe(lat))

  lon = field&lon
  lon_size = dimsizes(lon)
  wgt_lon = new((/lon_size(0)/), float)
  wgt_lon = 1.0

  ave = wgt_areaave_Wrap(field, wgt_lat, wgt_lon, 0)
  return(ave)

end
; #############################################################################
undef("coswgt_area_std")
function coswgt_area_std(field:numeric)
;
; Arguments
;    field: numeric field.
;
; Return value
;    The area average using cosine lat weights.
;
; Description
;    Computes the area average using cosine lat weights and lon weights=1.
;
begin

  funcname = "coswgt_area_std"

   lat = field&lat
   lon = field&lon
   re   = 6.37122e06
   rad  = 4.0 * atan(1.0) / 180.0
   con  = re * rad
   clat = cos(lat * rad)           ; cosine of latitude
   dlon = (lon(2) - lon(1))        ; assume dlon is constant
   dlat = (lat(2) - lat(1))        ; assume dlat is constant

   dx   = con * dlon * clat        ; dx at each latitude
   dy   = con * dlat               ; dy is constant
   dydx = dy * dx                  ; dydx(nlat)
   nlat=dimsizes(lat)
   mlon=dimsizes(lon)
   wgt  = new((/nlat, mlon/), typeof(dydx))
   wgt=conform (wgt, dydx, 0)



   ave = wgt_areaave2(field, wgt, 0)
   d2 = (field - ave) ^ 2 ;
     ; Calculate weighted standard deviation
   ss=sum(d2)
   arg = wgt_areasum2(d2, wgt, 0)

   v1 = sum(where(ismissing(field), wgt@_FillValue, wgt))
   v2 = sum(where(ismissing(field), wgt@_FillValue, wgt) ^ 2)
   dd=(v1 / (v1 ^ 2 - v2))

   out = sqrt( arg*dd)

  return(out)


end
; ########################################################
begin
;************************************************
; read in netCDF file
;************************************************
a = addfile("/Users/irenecionni/obs_c3s/Tier3/CDS-GTO-ECV/OBS_ESACCI-OZONE_sat_L3_Amon_toz_199601-201612.nc","r")
;************************************************
; read in zonal winds
;************************************************
  toz=a->toz                                
  dim=dimsizes(toz)

  avail1=new((/dim(1),dim(2)/),typeof(toz))
     do ilat=0,dim(1)-1
      do ilon=0,dim(2)-1
          avail1(ilat,ilon)=num(.not.ismissing(toz(:,ilat,ilon)))/tofloat(dim(0))    
      end do
    end do

    avail2=new((/dim(0),dim(2)/),typeof(toz))
     do itim=0,dim(0)-1
      do ilon=0,dim(2)-1
          avail2(itim,ilon)=num(.not.ismissing(toz(itim,:,ilon)))/tofloat(dim(1))
      end do
    end do

    avail3=new((/dim(0),dim(1)/),typeof(toz))
     do itim=0,dim(0)-1
      do ilat=0,dim(1)-1
          avail3(itim,ilat)=num(.not.ismissing(toz(itim,ilat,:)))/tofloat(dim(2))
      end do
    end do
  avail1=where(avail1.le.0,0,avail1)
  avail2=where(avail2.le.0,0,avail2)
  avail3=where(avail3.le.0,0,avail3)
  avail1!0="lat"
  avail1!1="lon"
  avail1&lat=toz&lat
  avail1&lon=toz&lon

  avail2!0="time"
  avail2!1="lon"
  avail2&time=toz&time
  avail2&lon=toz&lon

  avail3!0="time"
  avail3!1="lat"
  avail3&time=toz&time
  avail3&lat=toz&lat

; weighted average and stdeviation
  avg_avail1=coswgt_area_average(avail1)
  avg_std1=coswgt_area_std(avail1)
  avg2_avail=dim_avg_n_Wrap(avail2,0)
  avg3_avail=dim_avg_n_Wrap(avail3,0)
  avail1=where(avail1.le.0,avail1@_FillValue,avail1)
  avail2=where(avail2.le.0,avail2@_FillValue,avail2)
  avail3=where(avail3.le.0,avail3@_FillValue,avail3)
;************************************************
; create plot
;************************************************
  wks = gsn_open_wks("png","avail_1_"+"CDS-GTO-ECV")        
  gsn_define_colormap(wks,"MPL_YlGn")

  res = True                                    ; plot mods desired
  res@cnFillOn         = True               ; color Fill 
  res@cnFillMode       = "RasterFill"       ; Raster Mode
  res@cnLinesOn        =  False             ; Turn off contour lines
  res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
  res@cnMinLevelValF       = 0.               ; set min contour level
  res@cnMaxLevelValF       = 1.               ; set max contour level
  res@cnLevelSpacingF      = 0.1               ; set contour spacing
  res@cnLabelBarEndStyle="ExcludeOuterBoxes"
  res@tiMainString         = "Data Availability as a fraction [1]"    ; plot title
  fmt = "%3.2f"
  res@gsnLeftString      = "mean: "+sprintf(fmt,avg_avail1)
  res@gsnRightString     = " std:"+sprintf(fmt,avg_std1)           ; plot center string

  plot = gsn_csm_contour_map(wks,avail1, res)  ; create plot
;************************************************
  wks1 = gsn_open_wks("png","avail_2_"+"CDS-GTO-ECV")    
  plot2=new(2,"graphic")
  gsn_define_colormap(wks1,"MPL_YlGn")
  res1 = True                                    ; plot mods desired
  res1@gsnDraw                  = False            ; don't draw yet
  res1@gsnFrame                 = False            ; don't advance frame yet

  res1@cnFillOn         = True               ; color Fill
  res1@cnFillMode       = "RasterFill"       ; Raster Mode
  res1@cnLinesOn        =  False             ; Turn off contour lines
  res1@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
  res1@cnMinLevelValF       = 0.               ; set min contour level
  res1@cnMaxLevelValF       = 1.               ; set max contour level
  res1@cnLevelSpacingF      = 0.1               ; set contour spacing
  res1@cnLabelBarEndStyle="ExcludeOuterBoxes"
  res1@lbOrientation = "vertical"
  res1@tiMainString         = "Data Availability as a fraction [1]"    ; plot title
  fmt = "%3.2f"
;  res@gsnLeftString      = "mean: "+sprintf(fmt,avg_avail)
;  res@gsnRightString     = " std:"+sprintf(fmt,avg_std)           ; plot center string
  res1@vpHeightF= 0.4                    ; change aspect ratio of plot
  res1@vpWidthF = 0.7                  
  resTick                = True
  resTick@ttmFormat      = "%N/%Y" 
  resTick@ttmAxis        = "YL"
  resTick@ttmMajorStride = 12
  res1@tmYLLabelFontHeightF=0.012
   res1@tmXBLabelFontHeightF=0.012
; Set resources necessary to customize Y axis labels
  time_axis_labels( avail2&time,  res1, resTick )
  flip_avail2=lonFlip(avail2)
  plot2(0) = gsn_csm_contour(wks1,flip_avail2, res1)  ; create plot
  flip_avg2_avail=lonFlip(avg2_avail)
  res11=True
  res11@gsnDraw                  = False            ; don't draw yet
  res11@gsnFrame                 = False            ; don't advance frame yet
  ;res@tiMainString  = "An xy plot Example"   ; title
  res11@tiYAxisString = "Mean Availability [1]"           ; y axis title
  res11@vpHeightF= 0.4                    ; change aspect ratio of plot
  res11@vpWidthF = 0.7
  ;res@tiXAxisString = "Time"                 ; x axis title
  res11@trXMinF      = -180. 
  res11@trXMaxF      = 180.
  res11@trYMinF      = 0.
  res11@trYMaxF      = 1.
  res11@xyLineColor      = "NavyBlue"
  res11@xyLineThicknessF = 3.5
  res11@tmYLLabelFontHeightF=0.012
   res11@tmXBLabelFontHeightF=0.012
;  res11@gsnAttachPlotsXAxis=True   
  plot2(1) = gsn_csm_xy(wks1,flip_avg2_avail&lon,flip_avg2_avail, res11)  ; create plot
;  newplot =  gsn_attach_plots(plot2(0),plot2(1),res1,res11)
 ;draw(plot2(0))

   gsn_panel(wks1,plot2,(/2,1/),False)
 frame(wks1)
;************************************************
  wks2 = gsn_open_wks("png","avail_3_"+"CDS-GTO-ECV")                ; send graphics to PNG file
  plot3=new(2,"graphic")
  gsn_define_colormap(wks2,"MPL_YlGn")
  res2 = True                                    ; plot mods desired
  res2@gsnDraw                  = False            ; don't draw yet
  res2@gsnFrame                 = False            ; don't advance frame yet

  res2@cnFillOn         = True               ; color Fill
  res2@cnFillMode       = "RasterFill"       ; Raster Mode
  res2@cnLinesOn        =  False             ; Turn off contour lines
  res2@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
  res2@cnMinLevelValF       = 0.               ; set min contour level
  res2@cnMaxLevelValF       = 1.               ; set max contour level
  res2@cnLevelSpacingF      = 0.1               ; set contour spacing
  ;res2@lbOrientation = "vertical"
  res2@tiMainString         = "Data Availability as a fraction [1]"    ; plot title
  ;res2@lbLabelBarOn = False
  res2@cnLabelBarEndStyle="ExcludeOuterBoxes"
  fmt = "%3.2f"
;  res@gsnLeftString      = "mean: "+sprintf(fmt,avg_avail)
;  res@gsnRightString     = " std:"+sprintf(fmt,avg_std)           ; plot center string
  res2@vpHeightF= 0.4                    ; change aspect ratio of plot
  res2@vpWidthF = 0.7
  res2@tmYLMode        = "Explicit"	
  res2@tmYLValues      = (/-60,-30,0,30,60/)
  res2@tmYLLabels      = (/"60S","30S","Eq","30N","60N"/)
  res2@tmYLLabelFontHeightF=0.015
  res2@tmXBLabelFontHeightF=0.015
  res2Tick                = True
  res2Tick@ttmFormat      = "%N/%Y"
  res2Tick@ttmAxis        = "XB"
  res2Tick@ttmMajorStride = 24
    res2@tmYLLabelFontHeightF=0.012
   res2@tmXBLabelFontHeightF=0.012
; Set resources necessary to customize Y axis labels
  time_axis_labels( avail3&time,  res2, res2Tick )
  plot3(0) = gsn_csm_lat_time(wks2,avail3(lat|:,time|:), res2)  ; create plot
  res21=True
  res21@gsnDraw                  = False            ; don't draw yet
  res21@gsnFrame                 = False            ; don't advance frame yet
  ;res@tiMainString  = "An xy plot Example"   ; title
  res21@tiXAxisString = "Mean Availability [1]"           ; y axis title
  res21@vpHeightF= 0.4                    ; change aspect ratio of plot
  res21@vpWidthF = 0.3
  ;res@tiXAxisString = "Time"                 ; x axis title
  res21@trYMinF      = -90.
  res21@trYMaxF      = 90.
  res21@trXMinF      = 0
  res21@trXMaxF      = 1
  res21@tmYLMode        = "Explicit"
  res21@tmYLValues      = (/-60,-30,0,30,60/)
  res21@tmYLLabels      = (/"60S","30S","Eq","30N","60N"/)
  res21@tmYLLabelFontHeightF=0.015
  res21@tmXBMode        = "Explicit"
  res21@tmXBValues      = (/0.,0.2,0.4,0.6,0.8,1.0/)
  res21@tmXBLabels      = (/"0.","0.2","0.4","0.6","0.8","1."/)
  res21@tmXBLabelFontHeightF=0.015

;  res1@gsnAttachPlotsXAxis=True
  res21@xyLineColor      = "NavyBlue"
  res21@xyLineThicknessF  =3.5
  plot3(1) = gsn_csm_xy(wks2,avg3_avail,avg3_avail&lat, res21)  ; create plot
;  newplot =  gsn_attach_plots(plot2(0),plot2(1),res1,res11)
 ;draw(plot2(0))
  resP=True
  ;resP@gsnPanelLabelBar    = True                ; add common colorbar
   gsn_panel(wks2,plot3,(/1,2/),resP)
 frame(wks2)
;************************************************
 ;delete(res11) 
 ;ts_regions=new((/2,dimsizes(tropoz&time)/),typeof(avail_ts_NAF))
;
;   ts_regions(1,:)=avail_ts_NAF
;   ts_regions(0,:)=avail_ts_Ar
;   ts_regions!0="region"
;   ts_regions!1="time" 
;   ts_regions&time=tropoz&time
; wks4 = gsn_open_wks("png","avail_ts_"+"CDS-GTO-ECV")                ; send graphics to PNG file
; ; plot4=new(2,"graphic")
;
; res11=True
;  ;res11@gsnDraw                  = False            ; don't draw yet
;  ;res11@gsnFrame                 = False            ; don't advance frame yet
;  res11@tiMainString  = "Availability time series"   ; title
;  res11@tiYAxisString = "Availability [1]"           ; y axis title
;  res11@vpHeightF= 0.4                    ; change aspect ratio of plot
;  res11@vpWidthF = 0.6
;  res11@trYMinF      = 0.
;  res11@trYMaxF      = 1.
;  res11@trXMinF      = min(tropoz&time)
;  res11@xyLineColors     = (/"NavyBlue","red"/)
;  res11@xyLineThicknessF = 3
;  res11@pmLegendDisplayMode    = "Always"            ; turn on legend
;  res11@pmLegendSide           = "Top"               ; Change location of
;  res11@pmLegendParallelPosF   = 1.2                 ; move units right
;  res11@pmLegendOrthogonalPosF =- 0.7                ; move units down
;  res11@pmLegendWidthF         = 0.08                ; Change width and
;  res11@pmLegendHeightF        = 0.18                ; height of legend.
;  res11@lgPerimOn              = False               ; turn off box around
;  res11@lgLabelFontHeightF     = .01                 ; label font height
;  res11@xyExplicitLegendLabels = (/"Arabian Peninsula", "Sahara Desert"/)         ; create explicit labels
;  resTick                = True
;  resTick@ttmFormat      = "%N/%Y"
;  resTick@ttmAxis        = "XB"
;  resTick@ttmMajorStride = 12
;  time_axis_labels( tropoz&time,  res11, resTick )
;  plot4 = gsn_csm_xy(wks4,tropoz&time,ts_regions(region|:,time|:), res11)  ; create plot


end
