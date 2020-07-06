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
;***********************************************
dir="/Users/irenecionni/"
 a = addfile(dir+"/obs_c3s/Tier3/CDS-GTO-ECV/OBS_ESACCI-OZONE_sat_L3_Amon_toz_199601-201612.nc" ,"r")
;************************************************
; read in zonal winds
;************************************************
  tropoz=a->toz                                
  dim=dimsizes(tropoz)
  clim1=clmMonTLL(tropoz)
  anom_clim1= calcMonAnomTLL(tropoz,clim1)
  anom_ann=month_to_annual(anom_clim1, 1)
  time=tropoz&time
  date= cd_calendar(time, 0)
  yy=toint(date(:,0)) 
  dim=dimsizes(anom_ann)
  ts1_avg=coswgt_area_average(tropoz)
 

  ts_distribution=new((/7,dimsizes(tropoz&time)/),typeof(tropoz))
  do itim=0,dimsizes(tropoz&time)-1
     stat=stat_dispersion(tropoz(itim,:,:), False )
 
     ts_distribution(1:2,itim)=stat(23:24)
     ts_distribution(3,itim)=stat(6)
     ts_distribution(4,itim)=stat(10)
     ts_distribution(5:6,itim)=stat(25:26)
     delete(stat)
  end do

  ts_distribution(0,:)=ts1_avg
  ts_distribution!0="stat"
  ts_distribution!1="time"
  ts_distribution&stat=(/"Mean","Lower 1.0%","Lower 5.0%","LowQuartile","HighQuartile","Upper 5.0%","Upper 1.0%"/)
  wks = gsn_open_wks("png","var_1_"+"OBS_CDS-CDS-GTO-ECV")                ; send graphics to PNG file
  gsn_define_colormap(wks,"NEO_div_vegetation_c")
  plot=new(dim(0),"graphic")
  
  do itim=0,dim(0)-1
 
     avg_clim1=coswgt_area_average(anom_ann(itim,:,:))
     avg_std1=coswgt_area_std(anom_ann(itim,:,:))

     res = True                                    ; plot mods desired
     res@gsnDraw                  = False            ; don't draw yet
     res@gsnFrame                 = False            ; don't advance frame yet

     res@cnFillOn         = True               ; color Fill 
     res@cnFillMode       = "RasterFill"       ; Raster Mode
     res@cnLinesOn        =  False             ; Turn off contour lines
     res@cnLineLabelsOn       =False
     res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
     res@cnMinLevelValF       = -40               ; set min contour level
     res@cnMaxLevelValF       = 40               ; set max contour level
     res@cnLevelSpacingF      = 4.               ; set contour spacing
     ;res@cnLabelBarEndStyle="ExcludeOuterBoxes"
     res@tiMainString         = "Annual map "+yy(12*itim)    ; plot title
     fmt = "%3.2f"
     res@gsnLeftString      = "[DU]"
     res@gsnRightString     = "mean: "+sprintf(fmt,avg_clim1)+" std:"+sprintf(fmt,avg_std1)           ; plot center string

     res@lbLabelBarOn = False
     plot(itim) = gsn_csm_contour_map(wks,anom_ann(itim,:,:), res)  ; create plot
  end do
  resP=True
  resP@gsnPanelLabelBar    = True                ; add common colorbar
  gsn_panel(wks,plot,(/7,3/),resP)



;***********************time serie************************************************


  wks1 = gsn_open_wks("png","ts_1_"+"OBS_CDS-CDS-GTO-ECV")                ; send graphics to PNG file

; Set resources necessary to customize Y axis labels
  
  res11=True
  ;res11@gsnDraw                  = False            ; don't draw yet
  ;res11@gsnFrame                 = False            ; don't advance frame yet
  ;res@tiMainString  = "An xy plot Example"   ; title
  res11@tiYAxisString = "Total Ozone Column [DU]"           ; y axis title
  res11@vpHeightF= 0.4                    ; change aspect ratio of plot
  res11@vpWidthF = 0.6
  ;res@tiXAxisString = "Time"                 ; x axis title
  res11@trXMinF      = min(tropoz&time) 
  ;res11@trXMaxF      = 180.
  res11@trYMinF      = 250.
  res11@trYMaxF      = 350.
  res11@xyLineColor      = "NavyBlue"
  res11@xyLineThicknessF = 3.5
  resTick                = True
  resTick@ttmFormat      = "%N/%Y"
  resTick@ttmAxis        = "XB"
  resTick@ttmMajorStride = 24
  res11@tmYLLabelFontHeightF=0.01
  res11@tmXBLabelFontHeightF=0.01
  time_axis_labels( ts1_avg&time,  res11, resTick )
  plot2 = gsn_csm_xy(wks1,time,ts1_avg, res11)  ; create plot
  

;************************************************
;************************************************
;************************************************
  wks3 = gsn_open_wks("png","ts_stat_"+"OBS_CDS-CDS-GTO-ECV")                ; send graphics to PNG file

; Set resources necessary to customize Y axis labels
  delete(res11)
  res11=True
  res11@tiYAxisString = "Total Ozone Column [DU]"           ; y axis title
  res11@vpHeightF= 0.4                    ; change aspect ratio of plot
  res11@vpWidthF = 0.6
  res11@trXMinF      = min(tropoz&time)
  res11@trYMinF      = 100.
  res11@trYMaxF      = 500.
  ;res11@trYMinF      = 100.
  res11@trXMaxF      = max(ts_distribution&time)
  res11@xyLineColors     = (/"Black","LightBlue","NavyBlue","Purple","Purple","NavyBlue","LightBlue"/)
  res11@xyMonoLineThickness =False
  res11@xyLineThicknesses = (/3,2,2,2,2,2,2./)
  res11@pmLegendDisplayMode    = "Always"            ; turn on legend
  res11@pmLegendSide           = "Top"               ; Change location of
  res11@pmLegendParallelPosF   = 1.2                  ; move units right
  res11@pmLegendOrthogonalPosF = -0.6                ; move units down
  res11@pmLegendWidthF         = 0.1                ; Change width and
  res11@pmLegendHeightF        = 0.18                ; height of legend.
  res11@lgPerimOn              = False               ; turn off box around
  res11@lgLabelFontHeightF     = .01                 ; label font height
  res11@xyExplicitLegendLabels = ts_distribution&stat         ; create explicit labels
  resTick                = True
  res11@tmYLLabelFontHeightF=0.01
  res11@tmXBLabelFontHeightF=0.01
  resTick@ttmFormat      = "%N/%Y"
  resTick@ttmAxis        = "XB"
  resTick@ttmMajorStride =24 
  time_axis_labels( ts1_avg&time,  res11, resTick )
  plot3 = gsn_csm_xy(wks3,time,ts_distribution(stat|:,time|:), res11)  ; create plot
end
