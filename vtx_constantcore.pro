; Function VTX_ConstantCore
;
; Calculates semi-track length as a function of time and etching structure vT(x).
;
; INPUT:
;   time: Time since start of grain mount etching (seconds)
;   xInt: distance along track at which etchant channel intersects (µm)
;   ts: Time at which track started etching (seconds)
;   vB: Bulk etching rates (µm/s)
;   vTmax: Maximum track etching rate, at center of track (µm/s)
;   delXTo: Length of constant-etching region of track (µm)
;   delXTo_B: Distance from constant-vT core to track tip (µm)
;
; KEYWORD PARAMETERS:
;   END_VTVB: Set this to a variable to return the vT/vB ratio at the end of the etched semi-track
;
; OUTPUT:
;   Returns half-length of track (etching from xInt toward tip in positive x direction) (µm)
;
; MODIFICATION HISTORY:
;   Written by: Rich Ketcham, June 2020

Function VTX_ConstantCore, time, xInt, ts, vB, vTmax, delXTo, delXTo_B, END_VTVB=end_vTvB
; Set up variables
  xB1 = 0
  xT1 = xB1 + delXTo_B
  xT2 = xT1 + delXTo
  xB2 = xT2 + delXTo_B
  A = abs((vTMax-vB)/delXTo_B)
  
  if (xInt GE xT2) then begin   ; Zone 3
    if (time LT ts) then halfLen = 0 $
    else if (time LT ts+(1/A)*ALog((vTmax+A*(xT2-xInt))/vB)) then halfLen = (vTmax/A+xT2-xInt)*(1-EXP(-A*(time-ts))) $
    else halfLen = (time-ts-(1/A)*ALog((vTmax+A*(xT2-xInt))/vB))*vB+xB2-xInt
  endif else if (xInt GE xT1) then begin    ; Zone 2
    if (time LT ts) then halfLen = 0 $
    else if (time LT ts+(xT2-xInt)/vTmax) then halfLen = vTmax*(time-ts) $
    else if (time LT ts+(xT2-xInt)/vTmax-(1/A)*Alog(vB/vTmax)) then halfLen = xT2-xInt+(vTmax/A)*(1-EXP(-A*(time-ts-(xT2-xInt)/vTmax))) $
    else halfLen = xB2-xInt+vB*(time-ts-(xT2-xInt)/vTmax+(1/A)*Alog(vB/vTmax))
  endif else if (xInt GE xB1) then begin    ; Zone 1
    if (time LT ts ) then halfLen = 0 $
    else if (time LT ts-(1/A)*Alog(1+A*(xInt-xT1)/vTmax)) then halfLen = (vTmax/A+xInt-xT1)*(EXP(A*(time-ts))-1) $
    else if (time LT ts-(1/A)*Alog(1+A*(xInt-xT1)/vTmax)+delXTo/vTmax) then halfLen = xT1-xInt+vTmax*(time-ts+(1/A)*Alog(1+A*(xInt-xT1)/vTmax)) $
    else if (time LT ts-(1/A)*Alog(1+A*(xInt-xT1)/vTmax)+delXTo/vTmax-(1/A)*Alog(vB/vTmax)) then halfLen = xT2-xInt-(vTmax/A)*(EXP(-A*(time-ts+(1/A)*Alog(1+A*(xInt-xT1)/vTmax)-delXTo/vTmax))-1) $
    else halfLen = xB2-xInt+vB*(time-ts+(1/A)*Alog(1+A*(xInt-xT1)/vTmax)-delXTo/vTmax+(1/A)*Alog(vB/vTmax))
  endif

  if (Arg_Present(end_vTvB) OR (N_Elements(end_vTvB) GT 0)) then begin
    endPoint = xInt + halfLen
    if (endPoint LE xT1) then endVT = vB + (vTmax-vB)*(endPoint-xB1)/(xT1-xB1) $
    else if (endPoint LE xT2) then endVT = vTmax $
    else if (endPoint LE xB2) then endVT = vTmax - (vTmax-vB)*(endPoint-xT2)/(xB2-xT2) $
    else endVT = vB
    end_vTvB = endVT/vB
  endif

  return, halfLen
End

; Function VTXL_ConstantCore
;
; Calculates time required for a semi-track to reach a given length, given etching structure vT(x).
;
; INPUT:
;   length: Length to calculate (µm)
;   xInt: distance along track at which etchant channel intersects (µm)
;   ts: Time at which track started etching (seconds)
;   vB: Bulk etching rates (µm/s)
;   vTmax: Maximum track etching rate, at center of track (µm/s)
;   delXTo: Length of constant-etching region of track (µm)
;   delXTo_B: Distance from constant-vT core to track tip (µm)
;
; KEYWORD PARAMETERS:
;   END_VTVB: Set this to a variable to return the vT/vB ratio at the end of the etched semi-track
;   
; OUTPUT:
;   Returns etching time required to reach input semi-track length (seconds)
;
; MODIFICATION HISTORY:
;   Written by: Rich Ketcham, June 2020

Function VTXL_ConstantCore, length, xInt, ts, vB, vTmax, delXTo, delXTo_B, END_VTVB=end_vTvB
  xB1 = 0
  xT1 = xB1 + delXTo_B
  xT2 = xT1 + delXTo
  xB2 = xT2 + delXTo_B
  A = abs((vTMax-vB)/delXTo_B)
  
  if (xInt GE xT2) then begin   ; Zone 3
    if (length LT xB2-xInt) then time = (-1/A)*Alog(1-length/(vTmax/A+xT2-xInt)) $
    else time = (-1/A)*(Alog(vTmax/A+xT2-xB2)-Alog(vTmax/A+xT2-xInt))+(length-(xB2-xInt))/vB
  endif else if (xInt GE xT1) then begin    ; Zone 2
    if (length LT xT2-xInt) then time = length/vTmax $
    else begin
      time = (xT2-xInt > 0)/vTmax
      if (length LT xB2-xInt) then time += (-1/A)*Alog(1+(xT2-(xInt+length))/(vTmax/A)) $
      else time += (-1/A)*Alog(vB/vTmax)+(length-(xB2-xInt))/vB
    endelse
  endif else if (xInt GE xB1) then begin    ; Zone 1
    if (length LT xT1-xInt) then time = (1/A)*Alog(1+length/(vTmax/A+xInt-xT1)) $
    else begin
      time = (-1/A)*Alog(1+(xInt-xT1)/(vTmax/A))
      if (length LT xT2-xInt) then time += (length-(xT1-xInt))/vTmax $
      else begin
        time += delXTo/vTmax
        if (length LT xB2-xInt) then time += (-1/A)*Alog(1+(xT2-(xInt+length))/(vTmax/A)) $
        else time += (-1/A)*Alog(vB/vTmax)+(length-(xB2-xInt))/vB
      endelse
    endelse
  endif
  
  if Arg_Present(end_vTvB) then begin
    endPoint = xInt + length
    if (endPoint LE xT1) then endVT = vB + (vTmax-vB)*(endPoint-xB1)/(xT1-xB1) $
    else if (endPoint LE xT2) then endVT = vTmax $
    else if (endPoint LE xB2) then endVT = vTmax - (vTmax-vB)*(endPoint-xT2)/(xB2-xT2) $
    else endVT = vB
    end_vTvB = endVT/vB
  endif

  return, time + ts
End


; Pro PlotTipSequence
;
; Plots a sequence of track tips for a given etching model and set of etching times.
; Track widths are calculated as a prescribed spacing, and the tip location for 
; each etch time is calculated exactly.
; 
; At this point, it assumes just two etching rates, vT(x) along the track and vB for
; everything else.  A more sophisticated model incorporating track angle and using 
; the recent etching scheme by Jonckheere and others (2019, 2020) could make something
; more realistic-looking.  This would likely consist of getting a value for vR(phi),
; the radial growth rate, and having the track itersect the slow-etching intersection
; of basal and prismatic planes at the appropriate angle.
;
; INPUT:
;   vB: Bulk etching rate (µm/s)
;   vTmax: Maximum along-track etching rate (µm/s)
;   delXTo: Length of core zone with constant etch rate (µm)
;   delXTo_B: Length of tip zones (µm)
;
; KEYWORD_PARAMETERS:
;   TIMES: Array of times to calculate tips for (seconds, default=5,10,15,20,25)
;   X_STEP: Length grid spacing at which widths are calculated (µm)
;   
; OUTPUT
;   Creates an IDL PLOT of the sequence of track tips, and prints the vT/vB for 
;   the tip corresponding to each time step.
;
; MODIFICATION HISTORY
;   WRITTEN BY: Richard Ketcham, 18 June 2020


Pro PlotTipSequence, vB, vTmax, delXTo, delXTo_B, TIMES=times, X_STEP=xStep
  if (N_Params() LT 1) then vB = 0.022
  if (N_Params() LT 2) then vTmax = 0.83
  if (N_Params() LT 3) then delXTo = 10.5
  if (N_Params() LT 4) then delXTo_B = 3.0
  if NOT Keyword_Set(times) then times = (IndGen(5)+1)*5.
  if NOT Keyword_Set(xStep) then xStep = 0.1
  
  midPt = 0.5*delXTo + delXTo_B
  
  numTimes = N_Elements(times)
  lenTimes = FltArr(numTimes)
; Calculate tip locations for each etch time.
  for i=0,numTimes-1 do lenTimes[i] = VTX_ConstantCore(times[i], midPt, 0, vB, vTmax, delXTo, delXTo_B)
  
; Generate lengths at default spacing with tips locations spliced in
  lengths = [(IndGen(Ceil(lenTimes[numTimes-1]/xStep)) + 1)*xStep, lenTimes]
  lengths = lengths[Sort(lengths)]
  numLens = N_Elements(lengths)
  timeToLen = FltArr(numLens)    ; time that each length is reached
  tipVTVBs = FltArr(numLens)     ; Store tip information
  tipAlphas = FltArr(N_Elements(times))   
  
  for len=0,numLens-1 do begin
    timeToLen[len] = VTXL_ConstantCore(lengths[len], midPt, 0, vB, vTmax, delXTo, delXTo_B, END_VTVB=tipvTvB)
    tipVTVBs[len] = tipvTvB
  endfor

  widths = FltArr(numLens,numTimes)
  lensUsed = FltArr(numTimes)
  for timeInd = 0,numTimes-1 do begin
    lensUsed[timeInd] = (Where(lengths EQ lenTimes[timeInd]))[0]
    widths[0:lensUsed[timeInd],timeInd] = (times[timeInd]-timeToLen[0:lensUsed[timeInd]])*vB > 0
    tipAlphas[timeInd] = 2.*!RADEG*atan(widths[lensUsed[timeInd]-3,timeInd]/(lengths[lensUsed[timeInd]]-lengths[lensUsed[timeInd]-3]))
  endfor
  
  tipPlot = plot(lengths[0:lensUsed[0]],widths[0:lensUsed[0],0], XTITLE="Length (µm)", YTITLE="Thickness (µm)", ASPECT_RATIO=1, YTICKLEN=0.01)
  tipPlot = plot(lengths[0:lensUsed[0]],-widths[0:lensUsed[0],0], /OVERPLOT)
  for i=1,N_Elements(times)-1 do begin
    tipPlot = plot(lengths[0:lensUsed[i]], widths[0:lensUsed[i],i],/OVERPLOT)
    tipPlot = plot(lengths[0:lensUsed[i]], -widths[0:lensUsed[i],i],/OVERPLOT)
  endfor
  print, 'Tip vT/vB for each time: ',tipVTVBs[lensUsed]
  print, 'Alpha angles: ',tipAlphas
  
End


; Pro CalcRevelationRates
;
; Calculates confined track revelation due to etching a polished grain srface as a function of depth and 
; time due to etching model and track type.
; 
; METHOD
; 
; As semi-tracks penetrate to greater depths during etching, they intersect more volume beneath the 
; polished surface, allowing them to intersect confined tracks.  Once a track reaches a given 
; depth, it begins to widen, gradually increasing its likelihood of intersecting a confined track.
;
; Internal surface, external surface, and Cf semi-tracks will all have different effects on revelation rates.  Roughly  
; half of semi-tracks viewed on an internal surface will originate from above the polished surface, while all tracks 
; viewed on an external surface will originate below it.  Californium tracks will all have similar orientations and  
; penetration depths, whereas Uranium will be randomly orientated, with dip biasing based on 
; relative frequency (think of areas between latitude lines on a globe) and relative likelihood
; of intersecting the polished grain surface.
; 
; The program output is 
;
; INPUT:
;   vB: Bulk etching rate (µm/s)
;   vTmax: Maximum along-track etching rate (µm/s)
;   delXTo: Length of core zone with constant etch rate (µm)
;   delXTo_B: Length of tip zones (µm)
;   
; OUPUT:
; If "SHOW" is set, the routine creates plots of semi-track penetration and relative confined track revelation.
; If "REVCDF" is set to a variable, it returns a CDF of revelation probability as a function of time and depth, as 
; a FltArr with dimensions [depth, time].  First row is for depthInc, and first column is for 1*timeInc; "zero"
; depth and time are not included.  By generating a random number, and finding the 1D index of the first
; number greater than it in the CDF, the time and depth can be extracted based on the 2D array indices.
;
; KEYWORD PARAMETERS
;   EXTERNAL: Set to model tracks observed on external rather than internal surface
;   CF: Fraction of semi-tracks that are Cf tracks (0-1; default = 0)
;   SHOW: Set to create plots of outputs: semi-track penetration and relative confined track revelation.
;   DEPTH_INC: Depth increments at which to calculate revelation (µm)
;   TIME_INC: Time increments at which to calculate revelation (seconds)
;   DURATION: Duration of etching episode to reveal confined tracks (seconds)
;   REVCDF: Set to a named variable to receive the CDF of revelation rate with depth and time (FltArr[16/depthInc,20/timeInc])
;   NUM_TRACKS: Number of semi-tracks to generate (default 100,000)
;
; TO DO
; - Permit distribution of latent track lengths (Probably just by passing in arrays of vB, vTmax, etc.)
; 
; MODIFICATION HISTORY
;   WRITTEN BY: Richard Ketcham, June 2020


Pro CalcRevelationRates, vB, vTmax, delXTo, delXTo_B, EXTERNAL=external, CF=cf, SHOW=show, REVCDF=revCDF, $
  DEPTH_INC=depthInc, TIME_INC=timeInc, DURATION=duration, NUM_TRACKS=numTracks

  if (N_Params() LT 1) then vB = 0.022
  if (N_Params() LT 2) then vTmax = 0.83
  if (N_Params() LT 3) then delXTo = 10.5
  if (N_Params() LT 4) then delXTo_B = 3.0
  if NOT Keyword_Set(depthInc) then depthInc = 1  ; micrometers
  if NOT Keyword_Set(timeInc) then timeInc = 1  ; seconds
  if NOT Keyword_Set(duration) then duration = 20.
  if NOT Keyword_Set(numTracks) then numTracks = 100000L

  ts = 0.
  trueLen = delXTo + 2.*delXTo_B

  defTickLen = 0.02  ; Tickmark setting for all graphs; percentage of graph width or height

  if Keyword_Set(cf) then begin
    numCfTracks = Round(numTracks * cf)
    numUTracks = numTracks - numCfTracks
    cfLen = 5.9   ; microns; based on measurements of DUR tracks at UTCT (75 degree dip)
    cfLenSD = 1.4 ; microns
    cfDip = 75./!RADEG
    cfVTmax = 1.35 ; Default Cf etching model based on unannealed induced, constant-core model
    cfDelXTo = 3.69
    cfDelXTo_B = 6.69
    cfTrueLen = cfDelXTo + 2.*cfDelXTo_B
  endif else begin
    numUTracks = numTracks
    numCfTracks = 0
  endelse
  
  numDepths = Round(16/depthInc)
  numTimes = Round(duration/timeInc)
  
; Dip bias: area weighting (cos) * intersection (2 sin); integral of 2sin(d)cos(d) = sin^2(d)
  dips = ASin(sqrt(RandomU(seed,numTracks)))
  if Keyword_Set(cf) then dips[numUTracks:numTracks-1] = cfDip
  
; Intersection points.  Random along tracks for internal, on short side for external
  xInts = trueLen * RandomU(seed,numTracks)   
  if Keyword_Set(external) then xInts *= 0.5
  if Keyword_Set(cf) then xInts[numUTracks:numTracks-1] = cfTrueLen - cfLen + RandomN(seed,numCfTracks)*cfLenSD

  ndAtTime = FltArr(numDepths,numTimes)  
  revAtTime = FltArr(numDepths, numTimes)
  
  td = FltArr(numTracks,numDepths) ; time required to get to depths from depthInc to depthInc*numDepths microns for each track
  for d=0,numDepths-1 do begin
    for i=0L,numUTracks-1 do td[i,d] = VTXL_ConstantCore((d+1.)*depthInc/sin(dips[i]), xInts[i], ts, vB, vTmax, delXTo, delXTo_B)
    for i=numUTracks,numTracks-1 do td[i,d] = VTXL_ConstantCore((d+1.)*depthInc/sin(dips[i]), xInts[i], ts, vB, cfVTmax, cfDelXTo, cfDelXTo_B)
    for t=0,numTimes-1 do begin
      ndAtTime[d,t] = Total(td[*,d] LE (t+1)*timeInc)/numTracks
      revAtTime[d,t] = Total(((t+1)*timeInc-td[*,d]) > 0)
    endfor
    
  endfor
  
  nrevAtTime = revAtTime/Max(revAtTime)
  revCDF = nrevAtTime
  for i=1,N_Elements(revCDF)-1 do revCDF[i] += revCDF[i-1]
  revCDF /= revCDF[N_Elements(revCDF)-1]
  
  if Keyword_Set(show) then begin
    print, nrevAtTime
    ts = Round(1./timeInc)
    ndPlot = plot(ndAtTime[*,0],(IndGen(numDepths)+1)*depthInc, XTITLE="Relative Semi-Track Penetration", YTITLE="Depth (µm)", YRANGE=[20,0], XRANGE=[0,1], $
        XTICKLEN=defTickLen, YTICKLEN=defTickLen) ;, FONT_NAME="Calibri")
    for i=ts,numTimes-1,ts do ndPlot = plot(ndAtTime[*,i],(IndGen(numDepths)+1)*depthInc,/OVERPLOT)

    revPlot = plot(nRevAtTime[*,0],(IndGen(numDepths)+1)*depthInc, XTITLE="Relative Confined Track Revelation", YTITLE="Depth (µm)", YRANGE=[20,0], XRANGE=[0,1], $
        XTICKLEN=defTickLen, YTICKLEN=defTickLen) ;, FONT_NAME="Calibri")
    for i=ts,numTimes-1,ts do revPlot = plot(nRevAtTime[*,i],(IndGen(numDepths)+1)*depthInc,/OVERPLOT)
    
    plot, revCDF
  endif
  
end

; Pro GenerateTrackSet
;
; Generates a set of confined fission track lengths based on an etching structure, track distribution,
; etching schedule, and selection criteria
;
; METHOD:
; First call CalcRevelationRates to estimate 
; the relative probability of intersecting and revealing a confined track at depth.  Then models a 
; set of near-horizontal tracks, determining which ones remain confined (i.e. don't intersect the 
; surface during etching), and then which ones the observer will measure, based either on etching
; rates at track tips (standard procedure) or length (early-step-etch procedure).
;
; INPUT:
;   vB: Bulk etching rate (µm/s) 
;   vTmax: Maximum along-track etching rate (µm/s)
;   delXTo: Length of core zone with constant etch rate (µm)
;   delXTo_B: Length of tip zones (µm)
;
; KEYWORD PARAMETERS:
;   SHOW: Set to produce a series of plots of various outputs
;   EXTERNAL: Set to calculate revelation model using tracks measured on external rather than internal surface, 
;   CF: Set to calculate revelation model using the entered proportion of Cf tracks (0-1, default 0)
;   DEPTH_INC: Depth increments at which to calculate revelation (µm, default 1)
;   TIME_INC: Time increments at which to calculate revelation (seconds, default 1)
;   NUM_TRACKS: How many random tracks to generate (default 10,000)
;   MAX_DIP: Maximum dip for a confined track measurement (degrees, default 10)
;   EX_TIP_ZONE: Length of tip region to exclude from intersection (µm, default 1)
;   MAX_TIP_VTVB: Maximum etch rate ratio vT/vB for tip to be considered well etched (default 12)
;   UNDERETCH_BIAS: Power law term for simulating under-etched track selection bias (default 3.)
;   SURFACE_BUFFER: Minimum depth for track tip not to intersect surface (µm, default 0.5)
;   ETCH_TIME: Etching time if a single-step etch (seconds, default 20)
;   STEPS: Sequence of etching times, for step etch model (FltArr)
;   STEP_RESULTS: Set to variable to return mean confined track length for each etching step.
;   OPT_FIRST_STEP: Optimize tracks selected in first step to match mean track length (µm)
;   SD_LEN: Standard deviation of generated latent track lengths (µm)
;   CONTINUOUS_CULL: Set to cull semi-tracks only at the time steps they become exposed, matching Jonckheere et al. (2017) protocol
;
; TO DO
; - Add option for length biasing in track selection when there are a range of lengths. Matters <0.02 µm at SD=0.5, though
;
; MODIFICATION HISTORY
;   WRITTEN BY: Richard Ketcham, June 2020

Pro GenerateTrackSet, vB, vTmax, delXTo, delXTo_B, SHOW=show, EXTERNAL=external, NUM_TRACKS=nTracks, $
  DEPTH_INC=depthInc, TIME_INC=timeInc, MAX_TIP_VTVB=maxTipVTVB, ETCH_TIME=etchTime, SURFACE_BUFFER=surfaceBuffer, $
  MAX_DIP=maxDip, EX_TIP_ZONE = exTipZone, STEPS=steps, STEP_RESULTS=stepResults, CF=cf, OPT_FIRST_STEP=optFirstStep, $
  SD_LEN=sdLen, UNDERETCH_BIAS=underetchBias, CONTINUOUS_CULL=continuousCull

  if (N_Params() LT 1) then vB = 0.022
  if (N_Params() LT 2) then vTmax = 0.83
  if (N_Params() LT 3) then delXTo = 10.5
  if (N_Params() LT 4) then delXTo_B = 3.0
  
  if NOT Keyword_Set(nTracks) then nTracks = 10000  
  if NOT Keyword_Set(depthInc) then depthInc = 1  ; Depth incerment for revelation model (µm)
  if NOT Keyword_Set(timeInc) then timeInc = 1  ; Time increment for revelation model (seconds)
  if NOT Keyword_Set(maxDip) then maxDip = 10   ; Maximum dip for a confined track measurement (degrees)
  if NOT Keyword_Set(exTipZone) then exTipZone = 1.    ; Length of tip region to exclude from intersection (µm)
  if NOT Keyword_Set(maxTipVTVB) then maxTipVTVB = 12.   ; Maximum admissible vT/vB for a well-etched track tip
  if NOT Keyword_Set(underetchBias) then underetchBias = 3. ; Power law term for simulating under-etched track selection bias 
  if NOT Keyword_Set(surfaceBuffer) then surfaceBuffer = 0.5  ; Minimum depth for track tip not to intersect surface (µm).
                                                              ; Can also think of making it 2*vB*duration
  duration = Keyword_Set(etchTime) ? etchTime : Keyword_Set(steps) ? steps[0] : 20.
  
  numTimes = Round(duration/timeInc)
  numDepths = Round(16./depthInc)
  defTickLen = 0.02  ; Tickmark setting for all graphs; percentage of graph width or height
  
  CalcRevelationRates,  vB, vTmax, delXTo, delXTo_B, REVCDF=revCDF, EXTERNAL=external, CF=cf, DEPTH_INC=depthInc, TIME_INC=timeInc, DURATION=duration, SHOW=show
  cdfsz = Size(revCDF)
  cont = FltArr(cdfsz[2],cdfsz[1])  

; Randomize starting times and depths, based on revelation rate model
  tsGen = RandomU(seed,nTracks)
  ts = FltArr(nTracks)
  zInt = ts
  for i=0, nTracks-1 do begin
    ind = (Where(tsGen[i] LE revCDF))[0]
    inds = Array_Indices(revCDF,ind)
    ts[i] = (inds[1]+1)*timeInc
    zInt[i] = (inds[0]+1)*depthInc
    cont[inds[1],inds[0]] += 1
  endfor
  
  if Keyword_Set(show) then begin
;    tsdplot = Plot(ts, zInt, XTITLE="Start Time (s)", YTITLE="Depth (um)", YRANGE=[20,0], SYMBOL="+", LINESTYLE=6)
    ct = ColorTable(72, /reverse)
    titleLine = "All Tracks " + (Keyword_Set(external) ? "(External Surface)" : "") + " (N=" + StrTrim(String(nTracks),2) + ")"
    contPlot = Contour(cont,(IndGen(numTimes)+1)*timeInc,(IndGen(numDepths)+1)*depthInc,/FILL,RGB_TABLE=ct, YRANGE=[20,0], XTITLE="Start Time (s)", $
      YTITLE="Depth (µm)", TITLE=titleLine, POSITION=[0.08,0.17,0.96,0.93], XTICKLEN=defTickLen, YTICKLEN=defTickLen) ;  DIMENSION=[350,480]
    cb = ColorBar(TITLE='Number of Tracks',POSITION=[0.24,0.06,0.80,0.09])
  endif
  
; Randomize track dips (in radians)
  dips = ASin((RandomU(seed,nTracks))*sin(maxDip/!radeg))
;  plot, Histogram(dips*!radeg)
  
; Generate latent track lengths, and variation if requested
  if Keyword_Set(sdLen) then begin
    disp = RandomN(seed,nTracks)*sdLen
    coreFrac = delXTo/(delXTo + 2.*delXTo_B)  ; Distribute variation evenly between core and tip zones
    delXTos = disp*coreFrac + delXTo
    delXTo_Bs = disp*(1-coreFrac)/2. + delXTo_B
    latentLengths = delXTos + 2.*delXTo_Bs
  endif else latentLengths = FltArr(nTracks) + delXTo + 2.*delXTo_B
  
; Randomize intersection points
  xInts = RandomU(seed,nTracks)*(latentLengths - 2.*exTipZone) + exTipZone
  
; Now cull latent tracks that intersect surface
; Dips are defined to go down from intersection point, so xInt is length of up-facing segment
  if NOT Keyword_Set(continuousCull) then begin
    confinedTracks = Where(zInt GT xInts*sin(dips)+surfaceBuffer, nConfTracks)
    ts = ts[confinedTracks]
    zInt = zint[confinedTracks]
    dips = dips[confinedTracks]
    xInts = xInts[confinedTracks]
    latentLengths = latentLengths[confinedTracks]
    if Keyword_Set(sdLen) then begin
      delXTos = delXTos[confinedTracks]
      delXTo_Bs = delXTo_Bs[confinedTracks]
    endif
    
  ; Generate a new contour diagram with only confined tracks
    confCont = cont*0.
    for i=0, nConfTracks-1 do confCont[(ts[i]/timeInc)-1,Round((zint[i]/depthInc))-1] += 1
    if Keyword_Set(show) then begin
  ;    tsdplot = Plot(ts, zInt, XTITLE="Start Time (s)", YTITLE="Depth (um)", YRANGE=[20,0], SYMBOL="+", LINESTYLE=6)
      titleLine = "Confined Tracks " + (Keyword_Set(external) ? "(External Surface)" : "") + " (N=" + StrTrim(String(nConfTracks),2) + ")"
      confContPlot = Contour(confCont,(IndGen(numTimes)+1)*timeInc,(IndGen(numDepths)+1)*depthInc,/FILL,RGB_TABLE=ct, YRANGE=[20,0], XTITLE="Start Time (s)", $
        YTITLE="Depth (µm)", TITLE=titleLine, POSITION=[0.08,0.17,0.96,0.93], XTICKLEN=defTickLen, YTICKLEN=defTickLen)
      cb = ColorBar(TITLE='Number of Tracks',POSITION=[0.24,0.06,0.80,0.09])
    endif
  endif else nConfTracks = nTracks
  
; Calculate confined track lengths and tip etch ratios
  tip1s = FltArr(nConfTracks)
  tip2s = tip1s
  etchedLengths = tip1s
  hLen1s = tip1s
  for i=0, nConfTracks-1 do begin
    if Keyword_Set(sdLen) then begin
      hLen2 = VTX_ConstantCore(duration, xInts[i], ts[i], vB, vTmax, delXTos[i], delXTo_Bs[i], END_VTVB=tip2)
      hLen1 = VTX_ConstantCore(duration, latentLengths[i]-xInts[i], ts[i], vB, vTmax, delXTos[i], delXTo_Bs[i], END_VTVB=tip1)
    endif else begin
      hLen2 = VTX_ConstantCore(duration, xInts[i], ts[i], vB, vTmax, delXTo, delXTo_B, END_VTVB=tip2)
      hLen1 = VTX_ConstantCore(duration, latentLengths[i]-xInts[i], ts[i], vB, vTmax, delXTo, delXTo_B, END_VTVB=tip1)
    endelse 
    tip2s[i] = tip2
    tip1s[i] = tip1
    hLen1s[i] = hLen1
    etchedLengths[i] = hLen1 + hLen2
  endfor
  
  if Keyword_Set(continuousCull) then begin ; Cull only tracks that become invalid during this step
    confinedTracks = Where(zInt GT hlen1s*sin(dips)+2.0*duration*vB, nConfTracks)
    ts = ts[confinedTracks]
    zInt = zint[confinedTracks]
    dips = dips[confinedTracks]
    xInts = xInts[confinedTracks]
    latentLengths = latentLengths[confinedTracks]
    etchedLengths = etchedLengths[confinedTracks]
    tip1s = tip1s[confinedTracks]
    tip2s = tip2s[confinedTracks]
    hLen1s = hLen1s[confinedTracks]
    if Keyword_Set(sdLen) then begin
      delXTos = delXTos[confinedTracks]
      delXTo_Bs = delXTo_Bs[confinedTracks]
    endif
  endif
  
  if Keyword_Set(show) then begin
    elHist = Histogram(etchedLengths,BINSIZE=1, MIN=0)
    elHist[0] -= Total(etchedLengths EQ 0)
    elHistPlot = BarPlot(elHist,/HISTOGRAM, TITLE="Etched Lengths", $
      XTITLE="Length (µm)", YTITLE="Number", XTICKLEN=defTickLen, YTICKLEN=defTickLen)
    tipBinSize = 1
    maxTip = tip1s > tip2s
    maxTipHist = Histogram(maxTip,BINSIZE=tipBinSize, MIN=0)
    histSize = N_Elements(maxTipHist)
    useLog = maxTipHist[histSize-1]/maxTipHist[histSize-2] GT 4
    tip1HistPlot = BarPlot(IndGen(histSize)*tipBinSize,maxTipHist,/HISTOGRAM, TITLE="Tip Etching Ratios", $
      XTITLE="Maximum Tip vT/vB", YTITLE="Number", XTICKLEN=defTickLen, YTICKLEN=defTickLen, YLOG=useLog)
  endif
  
; Cull the under-etched ones
  if (duration LT 20) then $
    selectedTracks = Where(((((etchedLengths - 4.5)/5.5)>0)^underetchBias GT RandomU(seed,nConfTracks)) AND (ts LE (duration-3.)), nSelTracks) $
  else selectedTracks = Where((tip1s LE maxTipVTVB) AND (tip2s LE maxTipVTVB) AND (etchedLengths GT 4.0),nSelTracks)
  selLengths = etchedLengths[selectedTracks] 
  latentLengths = latentLengths[selectedTracks]
  dips = dips[selectedTracks]
  xInts = xInts[selectedTracks]
  ts = ts[selectedTracks]
  zInt = zint[selectedTracks]  
  tip1s = tip1s[selectedTracks]
  tip2s = tip2s[selectedTracks]
  hLen1s = hLen1s[selectedTracks]
  if Keyword_Set(sdLen) then begin
    delXTos = delXTos[selectedTracks]
    delXTo_Bs = delXTo_Bs[selectedTracks]
  endif

  
; If we need to fit a certain mean track length, do that culling, too
  if Keyword_Set(optFirstStep) then begin
    lenOrder = Sort(selLengths)
    numToGet = 200 < nSelTracks/5
    closestInd = (Where(selLengths[lenOrder] GT optFirstStep))[0]
    newSel = [closestInd]
    while N_Elements(newSel) LT numToGet-1 do begin
      if (Mean(selLengths[lenOrder[newSel]]) GT optFirstStep) then newSel = [newSel,Floor(RandomU(seed)*closestInd)] $
      else newSel = [newSel,closestInd + Floor(RandomU(seed)*(nSelTracks-closestInd))]
    endwhile
;    print, Mean(selLengths[newSel]), N_Elements(newSel)
    print, 'Bias = ',Mean(newSel)/(0.5*nSelTracks)   
    selLengths = selLengths[lenOrder[newSel]]
    latentLengths = latentLengths[lenOrder[newSel]]
    dips = dips[lenOrder[newSel]]
    xInts = xInts[lenOrder[newSel]]
    ts = ts[lenOrder[newSel]]
    zInt = zInt[lenOrder[newSel]]
    tip1s = tip1s[lenOrder[newSel]]
    tip2s = tip2s[lenOrder[newSel]]
    hLen1s = hLen1s[lenOrder[newSel]]
    if Keyword_Set(sdLen) then begin
      delXTos = delXTos[lenOrder[newSel]]
      delXTo_Bs = delXTo_Bs[lenOrder[newSel]]
    endif
    nSelTracks = N_Elements(newSel)
  endif
  
  print, "Etch step ",1, mean(selLengths), stddev(selLengths), nSelTracks, mean(latentLengths), mean(ts), mean(zInt)
  if Arg_Present(stepResults) then begin
    stepResults = FltArr(N_Elements(steps))
    stepResults[0] = mean(selLengths)
  endif
  
; Generate a final contour diagram with only selected tracks
  selCont = cont*0.
  for i=0, nSelTracks-1 do selCont[(ts[i]/timeInc)-1,(zint[i]/depthInc)-1] += 1
  if Keyword_Set(show) then begin
    titleLine = "Selected Tracks " + (Keyword_Set(external) ? "(External Surface)" : "") + " (N=" + StrTrim(String(nSelTracks),2) + ")"
    confContPlot = Contour(selCont,(IndGen(numTimes)+1)*timeInc,(IndGen(numDepths)+1)*depthInc,/FILL,RGB_TABLE=ct, YRANGE=[20,0], XTITLE="Start Time (s)", $
      YTITLE="Depth (µm)", TITLE=titleLine, POSITION=[0.08,0.17,0.96,0.93], XTICKLEN=defTickLen, YTICKLEN=defTickLen)
    cb = ColorBar(TITLE='Number of Tracks',POSITION=[0.24,0.06,0.80,0.09])
    elHistPlot = BarPlot(Histogram(selLengths,BINSIZE=1, MIN=0),/HISTOGRAM, TITLE="Etched and Selected Lengths", $
      XTITLE="Length (µm)", YTITLE="Number", XTICKLEN=defTickLen, YTICKLEN=defTickLen)
    elHistPlot = BarPlot(elHist,/HISTOGRAM,/OVERPLOT,TRANSPARENCY=50)
    
    lenInc = depthInc
    minLen = Min(selLengths, MAX=maxLen)
    minLenBin = Floor(minLen/lenInc)*lenInc
    numLenBins = Ceil((maxLen-minLenBin)/lenInc)+1
    lenDepthCont = FltArr(numLenBins,cdfsz[1])
    for i=0, nSelTracks-1 do lenDepthCont[(selLengths[i]-minLenBin)/lenInc,(zint[i]/depthInc)-1] += 1 
    lenDepthContPlot = Contour(lenDepthCont,(IndGen(numLenBins))*lenInc+minLenBin,(IndGen(numDepths)+1)*depthInc,/FILL,RGB_TABLE=ct, YRANGE=[20,0], $
      XTITLE="Length (µm)", YTITLE="Depth (µm)", TITLE=titleLine, POSITION=[0.08,0.17,0.96,0.93], XTICKLEN=defTickLen, YTICKLEN=defTickLen)
    cb = ColorBar(TITLE='Number of Tracks',POSITION=[0.24,0.06,0.80,0.09])
    
    zHist = Histogram(zInt,BINSIZE=depthInc,MIN=0)
    depthHistPlot = BarPlot(IndGen(N_Elements(zHist))*depthInc,zHist,/HISTOGRAM, TITLE="Intersection Depths", $
      XTITLE="Depth (µm)", YTITLE="Number", XTICKLEN=defTickLen, YTICKLEN=defTickLen)

; Comparable data plot
;   elHistPlot = BarPlot([Histogram(se2,BINSIZE=1, MIN=0),0],/HISTOGRAM, TITLE="Measured",XTITLE="Length (µm)", YTITLE="Number",XTICKLEN=0.02,YTICKLEN=0.02)

;    lenDepthPlot = Plot(selLengths,zInt,LINESTYLE=6,SYMBOL='x',XTITLE="Length (µm)",YTITLE="Depth (µm)",TITLE="Selected Length vs. Depth", $
;      YRANGE=[20,0])
  endif

  if Keyword_Set(steps) then begin
    if (N_Elements(steps) EQ 1) then return
    for step = 1,N_Elements(steps)-1 do begin
      for i=0, nSelTracks-1 do begin
        if Keyword_Set(sdLen) then begin
          hLen2 = VTX_ConstantCore(steps[step], xInts[i], ts[i], vB, vTmax, delXTos[i], delXTo_Bs[i], END_VTVB=tip2)
          hLen1 = VTX_ConstantCore(steps[step], latentLengths[i]-xInts[i], ts[i], vB, vTmax, delXTos[i], delXTo_Bs[i], END_VTVB=tip1)
        endif else begin
          hLen2 = VTX_ConstantCore(steps[step], xInts[i], ts[i], vB, vTmax, delXTo, delXTo_B, END_VTVB=tip2)
          hLen1 = VTX_ConstantCore(steps[step], latentLengths[i]-xInts[i], ts[i], vB, vTmax, delXTo, delXTo_B, END_VTVB=tip1)
        endelse
        tip2s[i] = tip2
        tip1s[i] = tip1
        hLen1s[i] = hLen1
        selLengths[i] = hLen1 + hLen2
      endfor
      if Keyword_Set(continuousCull) then begin ; Cull tracks that have become invalid this step
        confinedTracks = Where(zInt GT hLen1s*sin(dips)+2.0*steps[step]*vB, nConfTracks_2)
        ts = ts[confinedTracks]
        dips = dips[confinedTracks]
        zInt = zInt[confinedTracks]
        xInts = xInts[confinedTracks]
        selLengths = selLengths[confinedTracks]
        latentLengths = latentLengths[confinedTracks]
        tip1s = tip1s[confinedTracks]
        tip2s = tip2s[confinedTracks]
        hLen1s = hLen1s[confinedTracks]
        nSelTracks = N_Elements(selLengths)
        if Keyword_Set(sdLen) then begin
          delXTos = delXTos[confinedTracks]
          delXTo_Bs = delXTo_Bs[confinedTracks]
        endif
      endif

      print, "Etch step ",step+1, mean(selLengths), stddev(selLengths), N_Elements(selLengths), mean(latentLengths), mean(ts), mean(zInt)
      if Arg_Present(stepResults) then stepResults[step] = mean(selLengths)
      
      if Keyword_Set(show) AND (steps[step] EQ 20) then begin
        elHistPlot = BarPlot(Histogram(selLengths,BINSIZE=1, MIN=0),/HISTOGRAM, TITLE="Selected Lengths After 20 s", $
          XTITLE="Length (µm)", YTITLE="Number", XTICKLEN=defTickLen, YTICKLEN=defTickLen)
      endif
    endfor
    
  endif
End


; Function VTXFit1a
;
; Evaluation function for fitting of constant-core model for experiments SE1 and SE2. 
;
; METHOD:
; Calls GenerateTrackSet for each data set, and compare results to measurements (hardwired).
;
; INPUT:
;   c: FltArr(3): vTmax (µm/s), delXTo (µm), mean latent length (µm)
;
; COMMON BLOCKS
;   fitRecord: An optional set of modifiers to vary fit from default parameters, including:
;     fileUnit: File unit number for capturing each function call
;     tdRes: Time and depth resolution for track revelation CDF; 0.5 means every 0.5 seconds in etch time and 0.5 microns in depth (default 0.5)
;     nTracks: How many random internal tracks to generate (default 100000L)
;     exTip: Length of tip region to exclude from intersection (µm, default 1)
;     maxDip: Maximum dip of internal tracks (degrees, default 25)
;     maxTipVTVB: Maximum etch rate ratio vT/vB for tip to be considered well etched (default 12)
;     underetchBias: Power law term for simulating under-etched track selection bias (default 3.)
;     sdLen: Standard deviation of generated latent track lengths (µm, default 0)
;     chiSq: Set to 1 to return reduced chi squared; otherwise return sum of squared differences (0-1, default 0)
;
; MODIFICATION HISTORY
;   WRITTEN BY: Richard Ketcham, June 2020

Function VTXFit1a, c

  Common fitRecord, fileUnit, tdRes, nTracks, exTip, maxDip, maxTipVTVB, underetchBias, sdLen, chisq

  se2mtl = [9.89, 16.19, 16.99, 17.18]
  se2semtl = [0.175, 0.071, 0.073, 0.076]
  se1mtl = [15.77, 16.34, 16.92]
  se1semtl = [0.079, 0.084, 0.091]
  dof = 6
; Data including 15-second etching step
;  se2mtl = [9.89, 15.31, 16.19, 16.99, 17.18]
;  se2semtl = [0.175, 0.081, 0.071, 0.073, 0.076]
;  dof = 7

  delVTo_B = (c[2]-c[1])/2.

; Set defaults
  tdr = N_Elements(tdRes) ? tdRes : 0.5 
  ext = N_Elements(exTip) ? exTip : 2.
  nt = N_Elements(nTracks) ? nTracks : 100000L
  md = N_Elements(maxDip) ? maxDip : 25.
  mtvtvb = N_Elements(maxTipVTVB) ? maxTipVTVB : 12.
  ueb = N_Elements(underetchBias) ? underetchBias : 3.
  sdl = N_Elements(sdLen) ? sdLen : 0.   

  GenerateTrackSet, 0.022, c[0], c[1], delVTo_B, STEPS = [10,20,25,30], STEP_RESULTS=se2mod, DEPTH_INC=tdr, TIME_INC=tdr, EX_TIP=ext, NUM_TRACKS=nt, MAX_DIP=md, UNDERETCH_BIAS=ueb, SD_LEN=sdl
  GenerateTrackSet, 0.022, c[0], c[1], delVTo_B, STEPS = [20,25,30], STEP_RESULTS=se1mod, DEPTH_INC=tdr, TIME_INC=tdr, EX_TIP=ext, NUM_TRACKS=nt,MAX_DIP=md, MAX_TIP_VTVB=mtvtvb, SD_LEN=sdl

  if (chisq eq 1) then sumsq = (Total(((se1mtl - se1mod)/se1semtl)^2) + Total(((se2mtl - se2mod)/se2semtl)^2))/dof $
  else sumsq = Total((se1mtl - se1mod)^2) + Total((se2mtl - se2mod)^2)
  
  if (N_Elements(fileUnit) EQ 1) then PrintF, fileUnit, [c, sumSq] else print, [c, sumSq]
  return, sumsq
End


Pro RunFirstFit

  numIter = 20
  
  ; Set defaults
  tdr = 0.2
  ext = 2.
  nt = 100000L
  md = 25.
  mtvtvb = 12.
  ueb = 3.
  sdl = 0.5

  se1_2mtl = [15.77, 16.34, 16.92, 9.89, 16.19, 16.99, 17.18]
  se1_2semtl = [0.079, 0.084, 0.091, 0.175, 0.071, 0.073, 0.076]
  dof = 6
  c = [1.33396, 4.11996, 16.9621]
  delVTo_B = (c[2]-c[1])/2.
  se1res = FltArr(8,numIter)
  for i=0,numIter-1 do begin
    GenerateTrackSet, 0.022, c[0], c[1], delVTo_B, STEPS = [20,25,30], STEP_RESULTS=se1mod, DEPTH_INC=tdr, TIME_INC=tdr, EX_TIP=ext, NUM_TRACKS=nt,MAX_DIP=md, MAX_TIP_VTVB=mtvtvb, SD_LEN=sdl
    GenerateTrackSet, 0.022, c[0], c[1], delVTo_B, STEPS = [10,20,25,30], STEP_RESULTS=se2mod, DEPTH_INC=tdr, TIME_INC=tdr, EX_TIP=ext, NUM_TRACKS=nt, MAX_DIP=md, UNDERETCH_BIAS=ueb, SD_LEN=sdl
    se1res[0:2,i] = se1mod
    se1res[3:6,i] = se2mod
    se1res[7,i] = Total(((se1_2mtl - se1res[0:6,i])/se1_2semtl)^2)/dof
  endfor

  se3mtl = [9.11, 13.03, 14.89, 15.43, 15.69, 14.43]
  se3semtl = [0.303, 0.284, 0.107, 0.110, 0.111, 0.083]
  dof = 5.
  c = [0.824168, 8.67620, 15.4636]
  delVTo_B = (c[2]-c[1])/2.
  se3res = FltArr(7,numIter)
  for i=0,numIter-1 do begin
    GenerateTrackSet, 0.022, c[0], c[1], delVTo_B, STEPS = [10,15,20,25,30], STEP_RESULTS=se3mod, DEPTH_INC=tdr, TIME_INC=tdr, EX_TIP=ext, NUM_TRACKS=nt, MAX_DIP=md, UNDERETCH_BIAS=ueb, SD_LEN=sdl, /CF
    GenerateTrackSet, 0.022, c[0], c[1], delVTo_B, STEPS = [20], STEP_RESULTS=tk2020mod, DEPTH_INC=tdr, TIME_INC=tdr, EX_TIP=ext, NUM_TRACKS=nt,MAX_DIP=md, MAX_TIP_VTVB=mtvtvb, SD_LEN=sdl, /CF
    se3res[0:4,i] = se3mod
    se3res[5,i] = tk2020mod
    se3res[6,i] = Total(((se3mtl - se3res[0:5,i])/se3semtl)^2)/dof
  endfor
  
  se4mtl = [12.33,12.64,12.93,11.25]
  se4semtl = [0.072, 0.070, 0.069,0.123]
  dof = 3.
  c = [3.35674,0.760329,12.4503]
  delVTo_B = (c[2]-c[1])/2.
  se4res = FltArr(5,numIter)
  for i=0,numIter-1 do begin
    GenerateTrackSet, 0.022, c[0], c[1], delVTo_B, STEPS = [20,25,30], STEP_RESULTS=se4mod, DEPTH_INC=tdr, TIME_INC=tdr, EX_TIP=ext, NUM_TRACKS=nt, MAX_DIP=md, UNDERETCH_BIAS=ueb, SD_LEN=sdl
    GenerateTrackSet, 0.022, c[0], c[1], delVTo_B, STEPS = [10], STEP_RESULTS=eae1mod, DEPTH_INC=tdr, TIME_INC=tdr, EX_TIP=ext, NUM_TRACKS=nt, MAX_DIP=md, UNDERETCH_BIAS=ueb, SD_LEN=sdl
    se4res[0:2,i] = se4mod
    se4res[3,i] = eae1mod
    se4res[4,i] = Total(((se4mtl - se4res[0:3,i])/se4semtl)^2)/dof
  endfor

  se5mtl = [12.58,13.48,13.87,14.15,11.76]
  se5semtl = [0.094, 0.090, 0.089, 0.089,0.140]
  dof = 4.
  c = [1.76682,11.3476,13.3913]
  delVTo_B = (c[2]-c[1])/2.
  se5res = FltArr(6,numIter)
  for i=0,numIter-1 do begin
    GenerateTrackSet, 0.022, c[0], c[1], delVTo_B, STEPS = [15,20,25,30], STEP_RESULTS=se5mod, DEPTH_INC=tdr, TIME_INC=tdr, EX_TIP=ext, NUM_TRACKS=nt, MAX_DIP=md, UNDERETCH_BIAS=ueb, SD_LEN=sdl
    GenerateTrackSet, 0.022, c[0], c[1], delVTo_B, STEPS = [10], STEP_RESULTS=eae2mod, DEPTH_INC=tdr, TIME_INC=tdr, EX_TIP=ext, NUM_TRACKS=nt, MAX_DIP=md, UNDERETCH_BIAS=ueb, SD_LEN=sdl
    se5res[0:3,i] = se5mod
    se5res[4,i] = eae2mod
    se5res[5,i] = Total(((se5mtl - se5res[0:4,i])/se5semtl)^2)/dof
  endfor

  se6mtl = [13.77,15.13,15.34,15.65,13.38]
  se6semtl = [0.132, 0.078, 0.079, 0.078, 0.125]
  dof = 4
  c = [4.03698,0.188842,15.1048]
  delVTo_B = (c[2]-c[1])/2.
  se6res = FltArr(6,numIter)
  for i=0,numIter-1 do begin
    GenerateTrackSet, 0.022, c[0], c[1], delVTo_B, STEPS = [15,20,25,30], STEP_RESULTS=se6mod, DEPTH_INC=tdr, TIME_INC=tdr, EX_TIP=ext, NUM_TRACKS=nt, MAX_DIP=md, UNDERETCH_BIAS=ueb, SD_LEN=sdl
    GenerateTrackSet, 0.022, c[0], c[1], delVTo_B, STEPS = [10], STEP_RESULTS=eae3mod, DEPTH_INC=tdr, TIME_INC=tdr, EX_TIP=ext, NUM_TRACKS=nt,MAX_DIP=md, UNDERETCH_BIAS=ueb, SD_LEN=sdl
    se6res[0:3,i] = se6mod
    se6res[4,i] = eae3mod
    se6res[5,i] = Total(((se6mtl - se6res[0:4,i])/se6semtl)^2)/dof
  endfor

  print, se1res
  print, se3res
  print, se4res
  print, se5res
  print, se6res

End

Pro RunSecondFit

  numIter = 20

  ; Set defaults
  tdr = 0.2
  ext = 2.
  nt = 100000L
  md = 25.
  mtvtvb = 12.
  ueb = 3.
  sdl = 0.5

  delVTo = 0.0

  se1_2mtl = [15.77, 16.34, 16.92, 9.89, 16.19, 16.99, 17.18]
  se1_2semtl = [0.079, 0.084, 0.091, 0.175, 0.071, 0.073, 0.076]
  dof = 6
  c = [1.69570,17.0191]
  delVTo_B = (c[1]-delVTo)/2.
  se1res = FltArr(8,numIter)
  for i=0,numIter-1 do begin
    GenerateTrackSet, 0.022, c[0], delVTo, delVTo_B, STEPS = [20,25,30], STEP_RESULTS=se1mod, DEPTH_INC=tdr, TIME_INC=tdr, EX_TIP=ext, NUM_TRACKS=nt,MAX_DIP=md, MAX_TIP_VTVB=mtvtvb, SD_LEN=sdl
    GenerateTrackSet, 0.022, c[0], delVTo, delVTo_B, STEPS = [10,20,25,30], STEP_RESULTS=se2mod, DEPTH_INC=tdr, TIME_INC=tdr, EX_TIP=ext, NUM_TRACKS=nt, MAX_DIP=md, UNDERETCH_BIAS=ueb, SD_LEN=sdl
    se1res[3:6,i] = se2mod
    se1res[7,i] = Total(((se1_2mtl - se1res[0:6,i])/se1_2semtl)^2)/dof
  endfor

  se3mtl = [9.11, 13.03, 14.89, 15.43, 15.69, 14.43]
  se3semtl = [0.303, 0.284, 0.107, 0.110, 0.111, 0.083]
  dof = 5.
  c = [1.54455,15.6235]
  delVTo_B = (c[1]-delVTo)/2.
  se3res = FltArr(7,numIter)
  for i=0,numIter-1 do begin
    GenerateTrackSet, 0.022, c[0], delVTo, delVTo_B, STEPS = [10,15,20,25,30], STEP_RESULTS=se3mod, DEPTH_INC=tdr, TIME_INC=tdr, EX_TIP=ext, NUM_TRACKS=nt, MAX_DIP=md, UNDERETCH_BIAS=ueb, SD_LEN=sdl, /CF
    GenerateTrackSet, 0.022, c[0], delVTo, delVTo_B, STEPS = [20], STEP_RESULTS=tk2020mod, DEPTH_INC=tdr, TIME_INC=tdr, EX_TIP=ext, NUM_TRACKS=nt,MAX_DIP=md, MAX_TIP_VTVB=mtvtvb, SD_LEN=sdl, /CF
    se3res[0:4,i] = se3mod
    se3res[5,i] = tk2020mod
    se3res[6,i] = Total(((se3mtl - se3res[0:5,i])/se3semtl)^2)/dof
  endfor

  se4mtl = [12.33,12.64,12.93,11.25]
  se4semtl = [0.072, 0.070, 0.069,0.123]
  dof = 3.
  c = [3.58683, 12.4494]
  delVTo_B = (c[1]-delVTo)/2.
  se4res = FltArr(5,numIter)
  for i=0,numIter-1 do begin
    GenerateTrackSet, 0.022, c[0], delVTo, delVTo_B, STEPS = [20,25,30], STEP_RESULTS=se4mod, DEPTH_INC=tdr, TIME_INC=tdr, EX_TIP=ext, NUM_TRACKS=nt, MAX_DIP=md, UNDERETCH_BIAS=ueb, SD_LEN=sdl
    GenerateTrackSet, 0.022, c[0], delVTo, delVTo_B, STEPS = [10], STEP_RESULTS=eae1mod, DEPTH_INC=tdr, TIME_INC=tdr, EX_TIP=ext, NUM_TRACKS=nt, MAX_DIP=md, UNDERETCH_BIAS=ueb, SD_LEN=sdl
    se4res[0:2,i] = se4mod
    se4res[3,i] = eae1mod
    se4res[4,i] = Total(((se4mtl - se4res[0:3,i])/se4semtl)^2)/dof
  endfor

  se5mtl = [12.58,13.48,13.87,14.15,11.76]
  se5semtl = [0.094, 0.090, 0.089, 0.089,0.140]
  dof = 4.
  c = [3.54415, 13.5862]
  delVTo_B = (c[1]-delVTo)/2.
  se5res = FltArr(6,numIter)
  for i=0,numIter-1 do begin
    GenerateTrackSet, 0.022, c[0], delVTo, delVTo_B, STEPS = [15,20,25,30], STEP_RESULTS=se5mod, DEPTH_INC=tdr, TIME_INC=tdr, EX_TIP=ext, NUM_TRACKS=nt, MAX_DIP=md, UNDERETCH_BIAS=ueb, SD_LEN=sdl
    GenerateTrackSet, 0.022, c[0], delVTo, delVTo_B, STEPS = [10], STEP_RESULTS=eae2mod, DEPTH_INC=tdr, TIME_INC=tdr, EX_TIP=ext, NUM_TRACKS=nt, MAX_DIP=md, UNDERETCH_BIAS=ueb, SD_LEN=sdl
    se5res[0:3,i] = se5mod
    se5res[4,i] = eae2mod
    se5res[5,i] = Total(((se5mtl - se5res[0:4,i])/se5semtl)^2)/dof
  endfor

  se6mtl = [13.77,15.13,15.34,15.65,13.38]
  se6semtl = [0.132, 0.078, 0.079, 0.078, 0.125]
  dof = 4
  c = [4.08988,15.1013]
  delVTo_B = (c[1]-delVTo)/2.
  se6res = FltArr(6,numIter)
  for i=0,numIter-1 do begin
    GenerateTrackSet, 0.022, c[0], delVTo, delVTo_B, STEPS = [15,20,25,30], STEP_RESULTS=se6mod, DEPTH_INC=tdr, TIME_INC=tdr, EX_TIP=ext, NUM_TRACKS=nt, MAX_DIP=md, UNDERETCH_BIAS=ueb, SD_LEN=sdl
    GenerateTrackSet, 0.022, c[0], delVTo, delVTo_B, STEPS = [10], STEP_RESULTS=eae3mod, DEPTH_INC=tdr, TIME_INC=tdr, EX_TIP=ext, NUM_TRACKS=nt,MAX_DIP=md, UNDERETCH_BIAS=ueb, SD_LEN=sdl
    se6res[0:3,i] = se6mod
    se6res[4,i] = eae3mod
    se6res[5,i] = Total(((se6mtl - se6res[0:4,i])/se6semtl)^2)/dof
  endfor
  
  print, se1res
  print, se3res
  print, se4res
  print, se5res
  print, se6res

End


; Function VTXFit1spon
;
; Evaluation function for fitting of constant-core model for experiments SE3 and TK20
; See description of VTXFit1a for further information
Function VTXFit1spon, c

  Common fitRecord, fileUnit, tdRes, nTracks, exTip, maxDip, maxTipVTVB, underetchBias, sdLen, chisq

  se3mtl = [9.11, 13.03, 14.89, 15.43, 15.69]
  se3semtl = [0.303, 0.284, 0.107, 0.110, 0.111]
  tk2020mtl = [14.43]
  tk2020semtl = [0.083]
  dof = 5.

  delVTo_B = (c[2]-c[1])/2.
  
  ; Set defaults
  tdr = N_Elements(tdRes) ? tdRes : 0.5
  ext = N_Elements(exTip) ? exTip : 2.
  nt = N_Elements(nTracks) ? nTracks : 100000L
  md = N_Elements(maxDip) ? maxDip : 25.
  mtvtvb = N_Elements(maxTipVTVB) ? maxTipVTVB : 12.
  ueb = N_Elements(underetchBias) ? underetchBias : 3.
  sdl = N_Elements(sdLen) ? sdLen : 0.   ; or 0.5?

; Both experiments Cf-irradiated, so include CF flag
  GenerateTrackSet, 0.022, c[0], c[1], delVTo_B, STEPS = [10,15,20,25,30], STEP_RESULTS=se3mod, DEPTH_INC=tdr, TIME_INC=tdr, EX_TIP=ext, NUM_TRACKS=nt, MAX_DIP=md, UNDERETCH_BIAS=ueb, SD_LEN=sdl, /CF
  GenerateTrackSet, 0.022, c[0], c[1], delVTo_B, STEPS = [20], STEP_RESULTS=tk2020mod, DEPTH_INC=tdr, TIME_INC=tdr, EX_TIP=ext, NUM_TRACKS=nt,MAX_DIP=md, MAX_TIP_VTVB=mtvtvb, SD_LEN=sdl, /CF

  if (chisq eq 1) then sumsq = (Total(((se3mtl - se3mod)/se3semtl)^2) + Total(((tk2020mtl - tk2020mod)/tk2020semtl)^2))/dof $
  else sumsq = Total((se3mtl - se3mod)^2) + Total((tk2020mtl - tk2020mod)^2)

  if (N_Elements(fileUnit) EQ 1) then PrintF, fileUnit, [c, sumSq] else print, [c, sumSq]
  return, sumsq
End

; Function VTXFit1SE6
;
; Evaluation function for fitting of constant-core model for experiments SE6 and EAE3
; See description of VTXFit1a for further information
Function VTXFit1SE6, c

  Common fitRecord, fileUnit, tdRes, nTracks, exTip, maxDip, maxTipVTVB, underetchBias, sdLen, chisq

  se6mtl = [13.77, 15.13, 15.34, 15.65]
  se6semtl = [0.132, 0.078, 0.079, 0.078]
  eae3mtl = [13.38]
  eae3semtl = [0.125]
  dof = 4

  delVTo_B = (c[2]-c[1])/2.

; Set defaults
  tdr = N_Elements(tdRes) ? tdRes : 0.5
  ext = N_Elements(exTip) ? exTip : 2.
  nt = N_Elements(nTracks) ? nTracks : 100000L
  md = N_Elements(maxDip) ? maxDip : 25.
  mtvtvb = N_Elements(maxTipVTVB) ? maxTipVTVB : 12.
  ueb = N_Elements(underetchBias) ? underetchBias : 3.
  sdl = N_Elements(sdLen) ? sdLen : 0.   ; or 0.5?

  GenerateTrackSet, 0.022, c[0], c[1], delVTo_B, STEPS = [15,20,25,30], STEP_RESULTS=se6mod, DEPTH_INC=tdr, TIME_INC=tdr, EX_TIP=ext, NUM_TRACKS=nt, MAX_DIP=md, UNDERETCH_BIAS=ueb, SD_LEN=sdl
  GenerateTrackSet, 0.022, c[0], c[1], delVTo_B, STEPS = [10], STEP_RESULTS=eae3mod, DEPTH_INC=tdr, TIME_INC=tdr, EX_TIP=ext, NUM_TRACKS=nt,MAX_DIP=md, UNDERETCH_BIAS=ueb, SD_LEN=sdl

  if (chisq eq 1) then sumsq = (Total(((se6mtl - se6mod)/se6semtl)^2) + Total(((eae3mtl - eae3mod)/eae3semtl)^2))/dof $
  else sumsq = Total((se6mtl - se6mod)^2) + Total((eae3mtl - eae3mod)^2)
  
  if (N_Elements(fileUnit) EQ 1) then PrintF, fileUnit, [c, sumSq] else print, [c, sumSq]
  return, sumsq
End

; Function VTXFit1SE5
;
; Evaluation function for fitting of constant-core model for experiments SE5 and EAE2
; See description of VTXFit1a for further information
Function VTXFit1SE5, c

  Common fitRecord, fileUnit, tdRes, nTracks, exTip, maxDip, maxTipVTVB, underetchBias, sdLen, chisq

  se5mtl = [12.58, 13.48, 13.87, 14.15]
  se5semtl = [0.094, 0.090, 0.089, 0.089]
  eae2mtl = [11.76]
  eae2semtl = [0.140]
  dof = 4.

  delVTo_B = (c[2]-c[1])/2.
  
; Initialize or use defaults
  tdr = N_Elements(tdRes) ? tdRes : 0.5
  ext = N_Elements(exTip) ? exTip : 2.
  nt = N_Elements(nTracks) ? nTracks : 100000L
  md = N_Elements(maxDip) ? maxDip : 25.
  mtvtvb = N_Elements(maxTipVTVB) ? maxTipVTVB : 12.
  ueb = N_Elements(underetchBias) ? underetchBias : 3.
  sdl = N_Elements(sdLen) ? sdLen : 0.   ; or 0.5?

  GenerateTrackSet, 0.022, c[0], c[1], delVTo_B, STEPS = [15,20,25,30], STEP_RESULTS=se5mod, DEPTH_INC=tdr, TIME_INC=tdr, EX_TIP=ext, NUM_TRACKS=nt, MAX_DIP=md, UNDERETCH_BIAS=ueb, SD_LEN=sdl
  GenerateTrackSet, 0.022, c[0], c[1], delVTo_B, STEPS = [10], STEP_RESULTS=eae2mod, DEPTH_INC=tdr, TIME_INC=tdr, EX_TIP=ext, NUM_TRACKS=nt, MAX_DIP=md, UNDERETCH_BIAS=ueb, SD_LEN=sdl

  if (chisq eq 1) then sumsq = (Total(((se5mtl - se5mod)/se5semtl)^2) + Total(((eae2mtl - eae2mod)/eae2semtl)^2))/dof $
  else sumsq = Total((se5mtl - se5mod)^2) + Total((eae2mtl - eae2mod)^2)
  
  if (N_Elements(fileUnit) EQ 1) then PrintF, fileUnit, [c, sumSq] else print, [c, sumSq]
  return, sumsq
End

; Function VTXFit1SE4
;
; Evaluation function for fitting of constant-core model for experiments SE4 and EAE1
; See description of VTXFit1a for further information
Function VTXFit1SE4, c

  Common fitRecord, fileUnit, tdRes, nTracks, exTip, maxDip, maxTipVTVB, underetchBias, sdLen, chisq

  se4mtl = [12.33,12.64,12.93]
  se4semtl = [0.072, 0.070, 0.069]
  eae1mtl = [11.25]
  eae1semtl = [0.123]
  dof = 3.

  delVTo_B = (c[2]-c[1])/2.

; Initialize or use defaults
  tdr = N_Elements(tdRes) ? tdRes : 0.5
  ext = N_Elements(exTip) ? exTip : 2.
  nt = N_Elements(nTracks) ? nTracks : 100000L
  md = N_Elements(maxDip) ? maxDip : 25.
  mtvtvb = N_Elements(maxTipVTVB) ? maxTipVTVB : 12.
  ueb = N_Elements(underetchBias) ? underetchBias : 3.
  sdl = N_Elements(sdLen) ? sdLen : 0.   

  GenerateTrackSet, 0.022, c[0], c[1], delVTo_B, STEPS = [20,25,30], STEP_RESULTS=se4mod, DEPTH_INC=tdr, TIME_INC=tdr, EX_TIP=ext, NUM_TRACKS=nt, MAX_DIP=md, UNDERETCH_BIAS=ueb, SD_LEN=sdl
  GenerateTrackSet, 0.022, c[0], c[1], delVTo_B, STEPS = [10], STEP_RESULTS=eae1mod, DEPTH_INC=tdr, TIME_INC=tdr, EX_TIP=ext, NUM_TRACKS=nt, MAX_DIP=md, UNDERETCH_BIAS=ueb, SD_LEN=sdl

  if (chisq eq 1) then sumsq = (Total(((se4mtl - se4mod)/se4semtl)^2) + Total(((eae1mtl - eae1mod)/eae1semtl)^2))/dof $ 
  else sumsq = Total((se4mtl - se4mod)^2) + Total((eae1mtl - eae1mod)^2)
  
  if (N_Elements(fileUnit) EQ 1) then PrintF, fileUnit, [c, sumSq] else print, [c, sumSq]
  return, sumsq
End


; Function VTXFit2
;
; Evaluation function for fitting of linear etching model for experiments SE1 and SE2.
;
; METHOD:
; Calls GenerateTrackSet for each data set, and compare results to measurements (hardwired).
;
; INPUT:
;   c: FltArr(2): vTmax (µm/s), mean latent length (µm)
;
; COMMON BLOCKS
;   fitRecord: An optional set of modifiers to vary fit from default parameters, including:
;     fileUnit: File unit number for capturing each function call
;     tdRes: Time and depth resolution for track revelation CDF; 0.5 means every 0.5 seconds in etch time and 0.5 microns in depth (default 0.5)
;     nTracks: How many random internal tracks to generate (default 100000L)
;     exTip: Length of tip region to exclude from intersection (µm, default 1)
;     maxDip: Maximum dip of internal tracks (degrees, default 25)
;     maxTipVTVB: Maximum etch rate ratio vT/vB for tip to be considered well etched (default 12)
;     underetchBias: Power law term for simulating under-etched track selection bias (default 3.)
;     sdLen: Standard deviation of generated latent track lengths (µm, default 0)
;     chiSq: Set to 1 to return reduced chi squared; otherwise return sum of squared differences (0-1, default 0)
;
; MODIFICATION HISTORY
;   WRITTEN BY: Richard Ketcham, June 2020

Function VTXFit2, c

  Common fitRecord, fileUnit, tdRes, nTracks, exTip, maxDip, maxTipVTVB, underetchBias, sdLen, chisq

  se2mtl = [9.89, 16.19, 16.99, 17.18]
  se2semtl = [0.175, 0.071, 0.073, 0.076]
  se1mtl = [15.77, 16.34, 16.92]
  se1semtl = [0.079, 0.084, 0.091]
  dof = 6.
; Data including 15-second etching step
;  se2mtl = [9.89, 15.31, 16.19, 16.99, 17.18]
;  se2semtl = [0.175, 0.081, 0.071, 0.073, 0.076]
;  dof = 7

  delVTo = 0.0
  delVTo_B = (c[1]-delVTo)/2.
  
; Initialize or use defaults
  tdr = N_Elements(tdRes) ? tdRes : 0.5
  ext = N_Elements(exTip) ? exTip : 2.
  nt = N_Elements(nTracks) ? nTracks : 100000L
  md = N_Elements(maxDip) ? maxDip : 25.
  mtvtvb = N_Elements(maxTipVTVB) ? maxTipVTVB : 12.
  ueb = N_Elements(underetchBias) ? underetchBias : 3.
  sdl = N_Elements(sdLen) ? sdLen : 0.   ; or 0.5?

  GenerateTrackSet, 0.022, c[0], delVTo, delVTo_B, STEPS = [10,20,25,30], STEP_RESULTS=se2mod, DEPTH_INC=tdr, TIME_INC=tdr, EX_TIP=ext, NUM_TRACKS=nt, MAX_DIP=md, UNDERETCH_BIAS=ueb, SD_LEN=sdl
  GenerateTrackSet, 0.022, c[0], delVTo, delVTo_B, STEPS = [20,25,30], STEP_RESULTS=se1mod, DEPTH_INC=tdr, TIME_INC=tdr, EX_TIP=ext, NUM_TRACKS=nt, MAX_DIP=md, MAX_TIP_VTVB=mtvtvb, SD_LEN=sdl
  
  if (chisq eq 1) then sumsq = (Total(((se1mtl - se1mod)/se1semtl)^2) + Total(((se2mtl - se2mod)/se2semtl)^2))/dof $
  else sumsq = Total((se1mtl - se1mod)^2) + Total((se2mtl - se2mod)^2)
  
  if (N_Elements(fileUnit) EQ 1) then PrintF, fileUnit, [c, sumSq] else print, [c, sumSq]
  return, sumsq
End

; Function VTXFit2
;
; Evaluation function for fitting of linear etching model for experiments SE3 and TK20.
; See description of VTXFit2 for further information
Function VTXFit2spon, c

  Common fitRecord, fileUnit, tdRes, nTracks, exTip, maxDip, maxTipVTVB, underetchBias, sdLen, chisq

  se3mtl = [9.11, 13.03, 14.89, 15.43, 15.69]
  se3semtl = [0.303, 0.284, 0.107, 0.110, 0.111]
  tk2020mtl = [14.43]
  tk2020semtl = [0.083]
  dof = 5.

  delVTo = 0.0
  delVTo_B = (c[1]-delVTo)/2.
  
; Initialize or use defaults
  tdr = N_Elements(tdRes) ? tdRes : 0.5
  ext = N_Elements(exTip) ? exTip : 2.
  nt = N_Elements(nTracks) ? nTracks : 100000L
  md = N_Elements(maxDip) ? maxDip : 25.
  mtvtvb = N_Elements(maxTipVTVB) ? maxTipVTVB : 12.
  ueb = N_Elements(underetchBias) ? underetchBias : 3.
  sdl = N_Elements(sdLen) ? sdLen : 0.   

; Both experiments Cf-irradiated, so include CF flag
  GenerateTrackSet, 0.022, c[0], delVTo, delVTo_B, STEPS = [10,15,20,25,30], STEP_RESULTS=se3mod, DEPTH_INC=tdr, TIME_INC=tdr, EX_TIP=ext, NUM_TRACKS=nt, MAX_DIP=md, UNDERETCH_BIAS=ueb, SD_LEN=sdl, /CF
  GenerateTrackSet, 0.022, c[0], delVTo, delVTo_B, STEPS = [20], STEP_RESULTS=tk2020mod, DEPTH_INC=tdr, TIME_INC=tdr, EX_TIP=ext, NUM_TRACKS=nt, MAX_DIP=md, MAX_TIP_VTVB=mtvtvb, SD_LEN=sdl, /CF
  
  if (chisq eq 1) then sumsq = (Total(((se3mtl - se3mod)/se3semtl)^2) + Total(((tk2020mtl - tk2020mod)/tk2020semtl)^2))/dof $
  else sumsq = Total((se3mtl - se3mod)^2) + Total((tk2020mtl - tk2020mod)^2)

  if (N_Elements(fileUnit) EQ 1) then PrintF, fileUnit, [c, sumSq] else print, [c, sumSq]
  return, sumsq
End


; Function VTXFit2SE6
;
; Evaluation function for fitting of linear etching model for experiments SE6 and EAE3.
; See description of VTXFit2 for further information
Function VTXFit2SE6, c

  Common fitRecord, fileUnit, tdRes, nTracks, exTip, maxDip, maxTipVTVB, underetchBias, sdLen, chisq

  se6mtl = [13.77, 15.13, 15.34, 15.65]
  se6semtl = [0.132, 0.078, 0.079, 0.078]
  eae3mtl = [13.38]
  eae3semtl = [0.125]
  dof = 4.

  delVTo = 0.
  delVTo_B = (c[1]-delVTo)/2.
  
; Initialize or use defaults
  tdr = N_Elements(tdRes) ? tdRes : 0.5
  ext = N_Elements(exTip) ? exTip : 2.
  nt = N_Elements(nTracks) ? nTracks : 100000L
  md = N_Elements(maxDip) ? maxDip : 25.
  mtvtvb = N_Elements(maxTipVTVB) ? maxTipVTVB : 12.
  ueb = N_Elements(underetchBias) ? underetchBias : 3.
  sdl = N_Elements(sdLen) ? sdLen : 0.   ; or 0.5?

  GenerateTrackSet, 0.022, c[0], delVTo, delVTo_B, STEPS = [15,20,25,30], STEP_RESULTS=se6mod, DEPTH_INC=tdr, TIME_INC=tdr, EX_TIP=ext, NUM_TRACKS=nt, MAX_DIP=md, UNDERETCH_BIAS=ueb, SD_LEN=sdl
  GenerateTrackSet, 0.022, c[0], delVTo, delVTo_B, STEPS = [10], STEP_RESULTS=eae3mod, DEPTH_INC=tdr, TIME_INC=tdr, EX_TIP=ext, NUM_TRACKS=nt, MAX_DIP=md, UNDERETCH_BIAS=ueb, SD_LEN=sdl

  if (chisq eq 1) then sumsq = (Total(((se6mtl - se6mod)/se6semtl)^2) + Total(((eae3mtl - eae3mod)/eae3semtl)^2))/dof $
  else sumsq = Total((se6mtl - se6mod)^2) + Total((eae3mtl - eae3mod)^2)
  
  if (N_Elements(fileUnit) EQ 1) then PrintF, fileUnit, [c, sumSq] else print, [c, sumSq]

  return, sumsq
End


; Function VTXFit2SE5
;
; Evaluation function for fitting of linear etching model for experiments SE5 and EAE2.
; See description of VTXFit2 for further information
Function VTXFit2SE5, c

  Common fitRecord, fileUnit, tdRes, nTracks, exTip, maxDip, maxTipVTVB, underetchBias, sdLen, chisq

  se5mtl = [12.58, 13.48, 13.87, 14.15]
  se5semtl = [0.094, 0.090, 0.089, 0.089]
  eae2mtl = [11.76]
  eae2semtl = [0.140]
  dof = 4.

  delVTo = 0.
  delVTo_B = (c[1]-delVTo)/2.
  
; Initialize or use defaults
  tdr = N_Elements(tdRes) ? tdRes : 0.5
  ext = N_Elements(exTip) ? exTip : 2.
  nt = N_Elements(nTracks) ? nTracks : 100000L
  md = N_Elements(maxDip) ? maxDip : 25.
  mtvtvb = N_Elements(maxTipVTVB) ? maxTipVTVB : 12.
  ueb = N_Elements(underetchBias) ? underetchBias : 3.
  sdl = N_Elements(sdLen) ? sdLen : 0.   ; or 0.5?

  GenerateTrackSet, 0.022, c[0], delVTo, delVTo_B, STEPS = [15,20,25,30], STEP_RESULTS=se5mod, DEPTH_INC=tdr, TIME_INC=tdr, EX_TIP=ext, NUM_TRACKS=nt, MAX_DIP=md, UNDERETCH_BIAS=ueb, SD_LEN=sdl
  GenerateTrackSet, 0.022, c[0], delVTo, delVTo_B, STEPS = [10], STEP_RESULTS=eae2mod, DEPTH_INC=tdr, TIME_INC=tdr, EX_TIP=ext, NUM_TRACKS=nt, MAX_DIP=md, UNDERETCH_BIAS=ueb, SD_LEN=sdl

  if (chisq eq 1) then sumsq = (Total(((se5mtl - se5mod)/se5semtl)^2) + Total(((eae2mtl - eae2mod)/eae2semtl)^2))/dof $
  else sumsq = Total((se5mtl - se5mod)^2) + Total((eae2mtl - eae2mod)^2)
  
  if (N_Elements(fileUnit) EQ 1) then PrintF, fileUnit, [c, sumSq] else print, [c, sumSq]
  return, sumsq
End


; Function VTXFit2SE4
;
; Evaluation function for fitting of linear etching model for experiments SE4 and EAE1.
; See description of VTXFit2 for further information
Function VTXFit2SE4, c

  Common fitRecord, fileUnit, tdRes, nTracks, exTip, maxDip, maxTipVTVB, underetchBias, sdLen, chisq

  se4mtl = [12.33,12.64,12.93]
  se4semtl = [0.072, 0.070, 0.069]
  eae1mtl = [11.25]
  eae1semtl = [0.123]
  dof = 3.

  delVTo = 0.
  delVTo_B = (c[1]-delVTo)/2.
  
; Initialize or use defaults
  tdr = N_Elements(tdRes) ? tdRes : 0.5
  ext = N_Elements(exTip) ? exTip : 2.
  nt = N_Elements(nTracks) ? nTracks : 100000L
  md = N_Elements(maxDip) ? maxDip : 25.
  mtvtvb = N_Elements(maxTipVTVB) ? maxTipVTVB : 12.
  ueb = N_Elements(underetchBias) ? underetchBias : 3.
  sdl = N_Elements(sdLen) ? sdLen : 0.  

  GenerateTrackSet, 0.022, c[0], delVTo, delVTo_B, STEPS = [20,25,30], STEP_RESULTS=se4mod, DEPTH_INC=tdr, TIME_INC=tdr, EX_TIP=ext, NUM_TRACKS=nt, MAX_DIP=md, UNDERETCH_BIAS=ueb, SD_LEN=sdl
  GenerateTrackSet, 0.022, c[0], delVTo, delVTo_B, STEPS = [10], STEP_RESULTS=eae1mod, DEPTH_INC=tdr, TIME_INC=tdr, EX_TIP=ext, NUM_TRACKS=nt, MAX_DIP=md, UNDERETCH_BIAS=ueb, SD_LEN=sdl

  if (chisq eq 1) then sumsq = (Total(((se4mtl - se4mod)/se4semtl)^2) + Total(((eae1mtl - eae1mod)/eae1semtl)^2))/dof $ 
  else sumsq = Total((se4mtl - se4mod)^2) + Total((eae1mtl - eae1mod)^2)
  
  if (N_Elements(fileUnit) EQ 1) then PrintF, fileUnit, [c, sumSq] else print, [c, sumSq]
  return, sumsq
End


; Pro FirstFitTryManual
;
; Tries a manual fit of Common-core model to the SE1+SE2 data.  
;
; METHOD:
; Sets up the common block and initial guess and range variables, and calls the IDL Amoeba (simplex method) routine.
;
; INPUT:
;   initTry: FltArr(3): vTmax (µm/s), delXTo (µm), mean latent length (µm)
;
; COMMON BLOCKS
;   fitRecord: An optional set of modifiers to vary fit from default parameters, including:
;     fileUnit: File unit number for capturing each function call
;     tdRes: Time and depth resolution for track revelation CDF; 0.5 means every 0.5 seconds in etch time and 0.5 microns in depth (default 0.5)
;     nTracks: How many random internal tracks to generate (default 100000L)
;     exTip: Length of tip region to exclude from intersection (µm, default 1)
;     maxDip: Maximum dip of internal tracks (degrees, default 25)
;     maxTipVTVB: Maximum etch rate ratio vT/vB for tip to be considered well etched (default 12)
;     underetchBias: Power law term for simulating under-etched track selection bias (default 3.)
;     sdLen: Standard deviation of generated latent track lengths (µm, default 0)
;     chiSq: Set to 1 to return reduced chi squared; otherwise return sum of squared differences (0-1, default 0)
;
; MODIFICATION HISTORY
;   WRITTEN BY: Richard Ketcham, June 2020

Pro FirstFitTryManual, initTry
  Common fitRecord, fileUnit, tdRes, nTracks, exTip, maxDip, maxTipVTVB, underetchBias, sdLen, chiSq
  maxTipVTVB = 12
  tdRes = 0.5
  exTip = 2.5
  chiSq = 1
  if (N_Params() EQ 0) then initTry = [1.0,9.,17.]
  scale = [0.1,1.5,0.1]
  result = Amoeba(1.e-4,SCALE=scale,P0=initTry,FUNCTION_NAME='VTXFit1a',FUNCTION_VALUE=fVal, NMAX=50)
  if N_Elements(result) EQ 1 then print, "Fail" $
  else print, result, fVal
End


; Pro FirstFitTry
;
; Tries a fit of Common-core model to the SE1+SE2 data.
; Uses all default parameters, as defined in function VTXFit1
;
; INPUT:
;   initTry: FltArr(3): vTmax (µm/s), delXTo (µm), mean latent length (µm)

Pro FirstFitTry, initTry
  if (N_Params() EQ 0) then initTry = [1.0,9.,17.]
  scale = [0.1,1.5,0.1]
  result = Amoeba(1.e-4,SCALE=scale,P0=initTry,FUNCTION_NAME='VTXFit1a',FUNCTION_VALUE=fVal, NMAX=50)
  if N_Elements(result) EQ 1 then print, "Fail" $
  else print, result, fVal
End

; Pro FirstFitTrySpon
;
; Tries a fit of Common-core model to the SE3+TK20 data.
; Uses all default parameters, as defined in function VTXFit1spon
Pro FirstFitTrySpon, initTry
  if (N_Params() EQ 0) then initTry = [0.8,9.,15.5]
  scale = [0.1,1.5,0.3]
  result = Amoeba(1.e-4,SCALE=scale,P0=initTry,FUNCTION_NAME='VTXFit1spon',FUNCTION_VALUE=fVal, NMAX=50)
  if N_Elements(result) EQ 1 then print, "Fail" $
  else print, result, fVal
End

Pro FirstFitTrySE6, initTry
  if (N_Params() EQ 0) then initTry = [2.0,11.6,15.0]
  scale = [0.1,1.,0.3]
  result = Amoeba(1.e-4,SCALE=scale,P0=initTry,FUNCTION_NAME='VTXFit1SE6',FUNCTION_VALUE=fVal, NMAX=50)
  if N_Elements(result) EQ 1 then print, "Fail" $
  else print, result, fVal
End

Pro FirstFitTrySE5, initTry
  if (N_Params() EQ 0) then initTry = [1.7,10.,13.8]
  scale = [0.1,1.5,0.3]
  result = Amoeba(1.e-4,SCALE=scale,P0=initTry,FUNCTION_NAME='VTXFit1SE5',FUNCTION_VALUE=fVal, NMAX=50)
  if N_Elements(result) EQ 1 then print, "Fail" $
  else print, result, fVal
End

Pro FirstFitTrySE4, initTry
  if (N_Params() EQ 0) then initTry = [1.7,10.,12.4]
  scale = [0.1,1.5,0.3]
  result = Amoeba(1.e-4,SCALE=scale,P0=initTry,FUNCTION_NAME='VTXFit1SE4',FUNCTION_VALUE=fVal, NMAX=50)
  if N_Elements(result) EQ 1 then print, "Fail" $
  else print, result, fVal
End

Pro SecondFitTry, initTry
  if (N_Params() EQ 0) then initTry = [1.5,17.2]
  scale = [0.1,0.1]
  result = Amoeba(1.e-4,SCALE=scale,P0=initTry,FUNCTION_NAME='VTXFit2',FUNCTION_VALUE=fVal, NMAX=50)
  if N_Elements(result) EQ 1 then print, "Fail" $
  else print, result, fVal
End

Pro SecondFitTrySpon, initTry
  if (N_Params() EQ 0) then initTry = [1.4,15.7]
  scale = [0.1,0.1]
  result = Amoeba(1.e-4,SCALE=scale,P0=initTry,FUNCTION_NAME='VTXFit2spon',FUNCTION_VALUE=fVal, NMAX=50)
  if N_Elements(result) EQ 1 then print, "Fail" $
  else print, result, fVal
End

Pro SecondFitTrySE6, initTry
  if (N_Params() EQ 0) then initTry = [3.9,15.1]
  scale = [0.1,0.1]
  result = Amoeba(1.e-4,SCALE=scale,P0=initTry,FUNCTION_NAME='VTXFit2SE6',FUNCTION_VALUE=fVal, NMAX=50)
  if N_Elements(result) EQ 1 then print, "Fail" $
  else print, result, fVal
End

Pro SecondFitTrySE5, initTry
  if (N_Params() EQ 0) then initTry = [3.5,13.8]
  scale = [0.1,0.1]
  result = Amoeba(1.e-4,SCALE=scale,P0=initTry,FUNCTION_NAME='VTXFit2SE5',FUNCTION_VALUE=fVal, NMAX=50)
  if N_Elements(result) EQ 1 then print, "Fail" $
  else print, result, fVal
End

Pro SecondFitTrySE4, initTry
  if (N_Params() EQ 0) then initTry = [3.5,12.5]
  scale = [0.1,0.1]
  result = Amoeba(1.e-4,SCALE=scale,P0=initTry,FUNCTION_NAME='VTXFit2SE4',FUNCTION_VALUE=fVal, NMAX=50)
  if N_Elements(result) EQ 1 then print, "Fail" $
  else print, result, fVal
End


; PRO RecordFits1
;
; Runs a series of fits from multiple starting points.

Pro RecordFits1
  Common fitRecord, fileUnit, tdRes, nTracks, exTip, maxDip, maxTipVTVB, underetchBias, sdLen, chisq

  fileUnit = 1  
  
  chisq = 1
  tdRes = 0.2
  
  maxTipVTVB = 12.
  underetchBias = 3.0
  sdLen = 0.5

  OpenW, fileUnit, "C:\Temp\VTXfitSE1+2_12-3e_x2.txt"
  FirstFitTry, [1.0,9.,16.8]
  FirstFitTry, [1.0,9.,17.2]
  FirstFitTry, [1.0,5.,16.8]
  FirstFitTry, [1.0,5.,17.2]
  FirstFitTry, [1.4,9.,16.8]
  FirstFitTry, [1.4,9.,17.2]
  FirstFitTry, [1.4,5.,16.8]
  FirstFitTry, [1.4,5.,17.2]
  FirstFitTry, [1.45,3.5,17.0]
  FirstFitTry, [1.6,1.,17.1]
  FirstFitTry, [1.35,3.5,16.9]
  FirstFitTry, [1.35,4.0,17.0]
  FirstFitTry, [1.55,1.1,17.0]
  FirstFitTry, [1.57,1.3,16.9]
  Close, fileUnit

  OpenW, fileUnit, "C:\Temp\VTXfitSE3_12-3e_x2.txt"
  FirstFitTrySpon, [0.7,10.,15.7]
  FirstFitTrySpon, [0.7,10.,15.2]
  FirstFitTrySpon, [0.7,5.,15.7]
  FirstFitTrySpon, [0.7,5.,15.2]
  FirstFitTrySpon, [1.4,10.,15.7]
  FirstFitTrySpon, [1.4,10.,15.2]
  FirstFitTrySpon, [1.4,5.,15.7]
  FirstFitTrySpon, [1.4,5.,15.2]
  FirstFitTrySpon, [1.2,3.,15.7]
  FirstFitTrySpon, [1.1,4.,15.7]
  FirstFitTrySpon, [1.05,5.,15.6]
  Close, fileUnit

  OpenW, fileUnit, "C:\Temp\VTXfitSE6_12-3e_x2.txt"
  FirstFitTrySE6, [3.2,8.,14.8]
  FirstFitTrySE6, [3.2,8.,15.3]
  FirstFitTrySE6, [3.2,2.,14.8]
  FirstFitTrySE6, [3.2,2.,15.3]
  FirstFitTrySE6, [2.4,8.,14.8]
  FirstFitTrySE6, [2.4,8.,15.3]
  FirstFitTrySE6, [2.4,2.,14.8]
  FirstFitTrySE6, [2.4,2.,15.3]
  FirstFitTrySE6, [3.0,5.,15.1]
  FirstFitTrySE6, [4.0,0.,15.1]
  FirstFitTrySE6, [2.7,3.35,15.1]
  FirstFitTrySE6, [2.6,3.3,15.2]
  FirstFitTrySE6, [3.7,1.0,15.1]
  FirstFitTrySE6, [3.6,1.0,15.2]
  FirstFitTrySE6, [2.8,5.0,15.2]
  FirstFitTrySE6, [2.7,5.1,15.1]
  FirstFitTrySE6, [3.9,0.5,15.2]
  FirstFitTrySE6, [3.8,0.6,15.1]
  FirstFitTrySE6, [3.5,2.0,15.2]
  FirstFitTrySE6, [3.55,1.3,15.1]
  FirstFitTrySE6, [3.4,2.7,15.2]
  FirstFitTrySE6, [3.15,3.75,15.1]
  FirstFitTrySE6, [2.9,4.7,15.1]
  FirstFitTrySE6, [3.35,2.75,15.1]
  FirstFitTrySE6, [3.3,2.7,15.0]
  FirstFitTrySE6, [3.15,3.5,15.1]
  FirstFitTrySE6, [3.1,3.7,15.0]
  FirstFitTrySE6, [3.1,4.,15.1]
  Close, fileUnit
End


Pro RecordFits1a
  Common fitRecord, fileUnit, tdRes, nTracks, exTip, maxDip, maxTipVTVB, underetchBias, sdLen, chisq

  fileUnit = 1

  chisq = 1
  tdRes = 0.2

  maxTipVTVB = 12.
  underetchBias = 3.0
  sdLen = 0.5

;  OpenW, fileUnit, "C:\Temp\VTXfitSE5_12-3e_x2.txt"
;  FirstFitTrySE5, [1.4,13.,13.8]
;  FirstFitTrySE5, [1.4,13.,13.0]
;  FirstFitTrySE5, [1.4,10.,13.8]
;  FirstFitTrySE5, [1.4,10.,13.0]
;  FirstFitTrySE5, [2.0,13.,13.8]
;  FirstFitTrySE5, [2.0,13.,13.0]
;  FirstFitTrySE5, [2.0,10.,13.8]
;  FirstFitTrySE5, [2.0,10.,13.0]
;  FirstFitTrySE5, [1.9,8.5,13.4]
;  FirstFitTrySE5, [1.8,8.3,13.6]
;  FirstFitTrySE5, [2.0,8.0,13.4]
;  FirstFitTrySE5, [2.1,8.1,13.6]
;  FirstFitTrySE5, [2.05,7.5,13.6]
;  FirstFitTrySE5, [2.1,7.4,13.4]
;  FirstFitTrySE5, [1.85,9.5,13.4]
;  FirstFitTrySE5, [1.95,8.7,13.5]
;  FirstFitTrySE5, [1.92,8.8,13.5]
;  FirstFitTrySE5, [1.9,8.9,13.5]
;  FirstFitTrySE5, [1.85,8.8,13.5]
;  FirstFitTrySE5, [1.95,8.2,13.5]
;  Close, fileUnit

  OpenW, fileUnit, "C:\Temp\VTXfitSE4_12-3e_x2x.txt"
;  FirstFitTrySE4, [1.4,11.,12.8]
;  FirstFitTrySE4, [1.4,11.,12.0]
;  FirstFitTrySE4, [1.4,7.,12.8]
;  FirstFitTrySE4, [1.4,7.,12.0]
;  FirstFitTrySE4, [2.0,11.,12.8]
;  FirstFitTrySE4, [2.0,11.,12.0]
;  FirstFitTrySE4, [2.0,7.,12.8]
;  FirstFitTrySE4, [2.0,7.,12.0]
;  FirstFitTrySE4, [1.9,8.,12.3]
;  FirstFitTrySE4, [2.45,4.,12.4]
;  FirstFitTrySE4, [2.7,3.,12.5]
;  FirstFitTrySE4, [3.,2.,12.5]
;  FirstFitTrySE4, [3.3,1.,12.5]
;  FirstFitTrySE4, [2.4,5.,12.4]
;  FirstFitTrySE4, [2.4,4.9,12.6]
;  FirstFitTrySE4, [2.6,3.5,12.4]
;  FirstFitTrySE4, [2.6,3.4,12.6]
;  FirstFitTrySE4, [2.85,2.5,12.4]
;  FirstFitTrySE4, [2.85,2.4,12.6]
;  FirstFitTrySE4, [3.2,1.5,12.4]
;  FirstFitTrySE4, [3.2,1.4,12.6]
;  FirstFitTrySE4, [2.3,5.0,12.4]
;  FirstFitTrySE4, [2.2,4.9,12.6]
;  FirstFitTrySE4, [2.5,3.5,12.4]
;  FirstFitTrySE4, [2.55,3.7,12.6]
;  FirstFitTrySE4, [3.1,1.5,12.4]
;  FirstFitTrySE4, [3.15,1.7,12.6]
;  FirstFitTrySE4, [2.8,2.8,12.4]
;  FirstFitTrySE4, [2.76,3.0,12.6]
;  FirstFitTrySE4, [1.8,9.5,12.3]
;  FirstFitTrySE4, [1.9,8.0,12.3]
;  FirstFitTrySE4, [1.85,9.0,12.3]
;  FirstFitTrySE4, [1.87,8.5,12.3]
;  FirstFitTrySE4, [2.95,2.25,12.4]
  FirstFitTrySE4, [3.4,1.,12.4]
  FirstFitTrySE4, [3.5,0.5,12.5]  
  FirstFitTrySE4, [3.6,0.0,12.4]
  Close, fileUnit
End


Pro RecordFits2
  Common fitRecord, fileUnit, tdRes, nTracks, exTip, maxDip, maxTipVTVB, underetchBias, sdLen, chisq

  fileUnit = 1

  tdRes = 0.2
  chisq = 1
  
  maxTipVTVB = 12.
  underetchBias = 3.0

  sdLen = 0.5
  OpenW, fileUnit, "C:\Temp\VTXfit2SE1+2_12-3e_x2.txt"
  SecondFitTry, [1.4,16.8]
  SecondFitTry, [1.4,17.2]
  SecondFitTry, [1.9,16.8]
  SecondFitTry, [1.9,17.2]
  Close, fileUnit
;

  OpenW, fileUnit, "C:\Temp\VTXfit2SE3_12-3e_x2.txt"
  SecondFitTrySpon, [1.3,16.1]
  SecondFitTrySpon, [1.3,15.2]
  SecondFitTrySpon, [1.8,16.1]
  SecondFitTrySpon, [1.8,15.2]
  Close, fileUnit

  OpenW, fileUnit, "C:\Temp\VTXfit2SE6_12-3e_x2.txt"
  SecondFitTrySE6, [3.5,14.8]
  SecondFitTrySE6, [3.5,15.3]
  SecondFitTrySE6, [4.5,14.8]
  SecondFitTrySE6, [4.5,15.3]
  Close, fileUnit

  OpenW, fileUnit, "C:\Temp\VTXfit2SE5_12-3e_x2.txt"
  SecondFitTrySE5, [3.0,13.9]
  SecondFitTrySE5, [3.0,13.2]
  SecondFitTrySE5, [4.1,13.9]
  SecondFitTrySE5, [4.1,13.2]
  Close, fileUnit

  OpenW, fileUnit, "C:\Temp\VTXfit2SE4_12-3e_x2.txt"
  SecondFitTrySE4, [3.0,12.9]
  SecondFitTrySE4, [3.0,12.1]
  SecondFitTrySE4, [4.1,12.9]
  SecondFitTrySE4, [4.1,12.1]
  Close, fileUnit
End


; Pro MakeModel1Plots
;
; Creates a set of plots showing teh result of one or more fits to the data, and reports
; the best-fit model parameters and their confidence limits.
;
; INPUT:
;   fileName: Full file path+name of the file with run results (from "RecordFits" routines)
;   confLimit: Merit function value estalishing confidence limits
;
; WRITTEN BY:
;   Richard Ketcham, July 2020

Pro MakeModel1Plots, fileName, confLimit
  
  fd2 = ReadModel1Fits(fileName)
  rgb = 13    ; RGB table to use
  defTickLen = 0.02  ; Tickmark setting for all graphs; percentage of graph width or height
  
; Calculate ranges to show, excluding bad results
  val = Where((fd2[1,*] LE fd2[2,*]) AND (fd2[1,*] GE 0) AND (fd2[2,*] GT 5), NCOMPLEMENT=numBad)
  minValidCore = Min(fd2[1,val],MAX=maxValidCore)
  coreRange = [Floor(minValidCore*10.-1)/10.,Ceil(maxValidCore*10.+1)/10.]
  minValidLen = Min(fd2[2,val],MAX=maxValidLen)
  lenRange = [Floor(minValidLen*10.-1)/10.,Ceil(maxValidLen*10.+1)/10.]
  
  print, numBad, " invalid results excluded."
  
  fd2g = fd2[*,val]
  
; Calculate color table 
  objVals = fd2g[3,*]
  objVals = 1/objVals
  objMaxLim = Max(objVals)
  objMinLim = 0

  confidenceLimit = N_Params() LT 2 ? objMaxLim*0.9 : 1/confLimit
  withinConf = Where(objVals GE confidenceLimit)
  confVals = fd2g[*,withinConf]
    
  minObj = Min(fd2g[3,*],minInd)
  print, "Best-fit model: ", fd2g[*,minInd]
  fitLabel = "Best: " + StrTrim(String(fd2g[0,minInd],FORMAT='(G10.5)'),2) 
  for i=1,3 do fitLabel += ", " + StrTrim(String(fd2g[i,minInd],FORMAT='(G10.5)'),2) 

  runLabel = 'Run: ' + File_Basename(fileName,'.txt')
  
  p1 = Plot(fd2g[0,*],fd2g[2,*],SYMBOL="circle",SYM_SIZE=0.5, /SYM_FILLED, LINESTYLE=6, RGB_TABLE=rgb, YRANGE=lenRange, $
    WINDOW_TITLE=fileName, XTITLE="Maximum Etching Velocity (µm/s)", YTITLE="Latent Track Length (µm)", $
    XTICKLEN=defTickLen, YTICKLEN=defTickLen, LAYOUT=[3,1,1], DIMENSION=[1200,420], MARGIN=[0.16,0.23,0.05,0.05], $
    VERT_COLORS=255*((objVals > objMinLim)-objMinLim)/(objMaxLim-objMinLim))
  p2 = Plot(fd2g[0,*],fd2g[1,*],SYMBOL="circle",SYM_SIZE=0.5, /SYM_FILLED, LINESTYLE=6, RGB_TABLE=rgb, YRANGE=coreRange, $
    WINDOW_TITLE=fileName, XTITLE="Maximum Etching Velocity (µm/s)", YTITLE="Core Zone Length (µm)", $
    XTICKLEN=defTickLen, YTICKLEN=defTickLen, /CURRENT, LAYOUT=[3,1,2], MARGIN=[0.16,0.23,0.05,0.05], $
    VERT_COLORS=255*((objVals > objMinLim)-objMinLim)/(objMaxLim-objMinLim))
  p3 = Plot(fd2g[1,*],fd2g[2,*],SYMBOL="circle",SYM_SIZE=0.5, /SYM_FILLED, LINESTYLE=6, RGB_TABLE=rgb, XRANGE=coreRange, $
    WINDOW_TITLE=fileName, YTITLE="Latent Track Length (µm)", XTITLE="Core Zone Length (µm)", YRANGE = lenRange, $
    XTICKLEN=defTickLen, YTICKLEN=defTickLen, /CURRENT, LAYOUT=[3,1,3], MARGIN=[0.16,0.23,0.05,0.05], $
    VERT_COLORS=255*((objVals > objMinLim)-objMinLim)/(objMaxLim-objMinLim))
  cb = ColorBar(TITLE='1/c2n',RANGE=[0,objMaxLim],POSITION=[0.4,0.08,0.6,0.12], FONT_NAME="Symbol")
  t1 = Text(0.01,0.08,runLabel,/CURRENT,FONT_SIZE=11) ;, FONT_STYLE='Bold')
  t2 = Text(0.01,0.03,fitLabel,/CURRENT,FONT_SIZE=11) ;, FONT_STYLE='Bold')
  p3d = Plot3D(Reform(fd2g[0,*]),Reform(fd2g[1,*]),Reform(fd2g[2,*]),SYMBOL="circle",SYM_SIZE=0.5, /SYM_FILLED, $
    LINESTYLE=6, RGB_TABLE=rgb, YRANGE=coreRange, ZRANGE=lenRange, DIMENSION = [500,550], $ 
    WINDOW_TITLE=fileName, XTITLE="Maximum Etching Velocity (µm/s)", YTITLE="Core Zone Length (µm)", ZTITLE="Latent Track Length (µm)", $
    XTICKLEN=defTickLen, YTICKLEN=defTickLen, VERT_COLORS=255*((objVals > objMinLim)-objMinLim)/(objMaxLim-objMinLim))
  cb3d = ColorBar(TITLE='1/c2n',RANGE=[0,objMaxLim],POSITION=[0.25,0.06,0.75,0.09], FONT_NAME="Symbol")

  tab = String(9B)
  print, "Confidence intervals"
  for i=0,2 do print, "c", StrTrim(String(i),2), ":", tab, min(confVals[i,*]), tab, max(confVals[i,*]), tab, (max(confVals[i,*])-min(confVals[i,*]))/2.

End

Pro MakeModel2Plot, fileName, confLimit

  fd2 = ReadModel2Fits(fileName)
  rgb = 13    ; RGB table to use
  defTickLen = 0.02  ; Tickmark setting for all graphs; percentage of graph width or height

  ; Calculate ranges to show, excluding bad results
  minValidLen = Min(fd2[1,*],MAX=maxValidLen)
  lenRange = [Floor(minValidLen*10.-1)/10.,Ceil(maxValidLen*10.+1)/10.]

  ; Calculate color table
  objVals = fd2[2,*]
  objVals = 1/objVals
  objMaxLim = Max(objVals)
  objMinLim = 0

  confidenceLimit = N_Params() LT 2 ? objMaxLim*0.9 : 1/confLimit
  withinConf = Where(objVals GE confidenceLimit)
  confVals = fd2[*,withinConf]
  
  minObj = Min(fd2[2,*],minInd)
  print, "Best-fit model: ", fd2[*,minInd]
  fitLabel = "Best: " + StrTrim(String(fd2[0,minInd],FORMAT='(G10.5)'),2)
  for i=1,2 do fitLabel += ", " + StrTrim(String(fd2[i,minInd],FORMAT='(G10.5)'),2)

  runLabel = 'Run: ' + File_Basename(fileName,'.txt')

  p1 = Plot(fd2[0,*],fd2[1,*],SYMBOL="circle",SYM_SIZE=0.5, /SYM_FILLED, LINESTYLE=6, RGB_TABLE=rgb, $ ;YRANGE=lenRange, $
    WINDOW_TITLE=fileName, XTITLE="Maximum Etching Velocity (µm/s)", YTITLE="Latent Track Length (µm)", $
    XTICKLEN=defTickLen, YTICKLEN=defTickLen, DIMENSIONS=[500,540], POSITION=[0.15,0.2,0.96,0.95], $
    VERT_COLORS=255*((objVals > objMinLim)-objMinLim)/(objMaxLim-objMinLim))
  t1 = Text(0.01,0.08,runLabel,/CURRENT,FONT_SIZE=10) ;, FONT_STYLE='Bold')
  t2 = Text(0.01,0.03,fitLabel,/CURRENT,FONT_SIZE=10) ;, FONT_STYLE='Bold')
  cb = ColorBar(TITLE='1/c2n',RANGE=[0,objMaxLim],POSITION=[0.5,0.07,0.96,0.11],FONT_NAME="Symbol")
  
  tab = String(9B)
  print, "Confidence intervals"
  for i=0,1 do print, "c", StrTrim(String(i),2), ":", tab, min(confVals[i,*]), tab, max(confVals[i,*]), tab, (max(confVals[i,*])-min(confVals[i,*]))/2.
End

Function ReadModel1Fits, fileName
  if (N_Params() EQ 0) then fileName =  "C:\Temp\VTXfit1.txt"
  OpenR, 1, fileName
  ReadF, 1, val1, val2, val3, val4
  data = [val1, val2, val3, val4]
  while not Eof(1) do begin
    ReadF, 1, val1, val2, val3, val4
    data = [data, val1, val2, val3, val4]
  endwhile
  Close, 1
  data = Reform(data, 4, N_Elements(data)/4)
  return, data
End

Function ReadModel2Fits, fileName
  if (N_Params() EQ 0) then fileName =  "C:\Temp\VTXfit2.txt"
  OpenR, 1, fileName
  ReadF, 1, val1, val2, val3
  data = [val1, val2, val3]
  while not Eof(1) do begin
    ReadF, 1, val1, val2, val3
    data = [data, val1, val2, val3]
  endwhile
  Close, 1
  data = Reform(data, 3, N_Elements(data)/3)
  return, data
End


; Pro GetSEData
;
; Returns a vector with the individual track lengths for any of the step-etch experiments.  Mainly used for generating
; histograms that can be compared to model outputs.
;
; NOTE: There have been some corrections (as of 7/17/2020) of a few of the lengths that have not been incorporated into
; this code. 
;
; WRITTEN BY:
;   Richard Ketcham, July 2020

Pro GetSEData, SE1_STEP1=se1, SE2_STEP1=se2, SE3_STEP1=se3, SE4_STEP1=se4, SE5_STEP1=se5, SE6_STEP1=se6, $
    EAE1_STEP1=eae1, EAE2_STEP2=eae2, EAE3_STEP3 = eae3
  se1 = [16.56713573, 16.12177484, 15.31196424, 16.92365283, 14.83430981, 14.6887074, 16.35433582, 17.07934905, $
    15.14346793, 16.84663539, 15.9589066, 16.51546536, 16.30647491, 16.03995771, 15.3217373, 15.27079729, 15.66648544, $
    15.65987582, 17.02070063, 16.0450187, 15.44875444, 15.66146532, 16.28037162, 15.18539875, 15.88822271, 15.52590365, $
    15.98835905, 14.92323343, 14.9394289, 15.00984081, 15.24592057, 16.66180425, 16.2927966, 14.62330388, 15.72957471, $
    16.15219864, 15.53901287, 16.65866093, 17.03844538, 14.85683964, 16.14256563, 15.22995565, 15.59724444, 14.89736624, $
    15.05013452, 14.85384869, 16.21435, 15.32341832, 15.44111343, 16.26296756, 15.06228681, 15.58267929, 16.53680441, $
    14.93107662, 15.33089854, 15.51430766, 15.7090436, 14.9961173, 15.12387917, 15.55788356, 16.47052962, 16.29466784, $
    16.48929437, 16.01361324, 14.94555533, 15.15365108, 16.42261032, 16.7749732, 15.4500128, 15.70630832, 17.07837489, 16.20716573]
  se2 = [8.004911783, 5.939631682, 8.43128564, 5.5755, 6.166371563, 6.332089072, 7.676233165, 7.478101343, 9.728042765, $
    5.931231475, 5.098779838, 5.49218758, 7.781710516, 5.947173111, 9.646573527, 12.27011881, 11.65902103, 6.265564074, $
    10.83653014, 11.10956841, 10.0244658, 12.8751484, 12.68449841, 11.34267879, 11.09302882, 9.888167131, 12.15134216, $
    9.704503336, 6.844, 8.605605481, 7.968519263, 11.77912196, 12.17623661, 7.692506286, 7.795207679, 8.624044532, $
    10.21963044, 10.14211513, 10.10376128, 9.494086813, 12.79946137, 9.224507298, 8.144418998, 9.581096911, 13.71692514, $
    11.14564433, 10.24336255, 8.452245373, 13.18825235, 9.545082818, 10.60442021, 9.797466293, 11.37942775, 11.72621076, $
    10.15471898, 11.32475646, 12.06268575, 8.828886164, 12.71650127, 11.44420234, 8.002673643, 10.66021651, 10.38264851, $
    11.78666483, 10.14555361, 7.847616986, 12.75614953, 9.735503553, 11.75987384, 10.15677925, 12.29391806, 10.75633323, $
    7.933596915, 8.685147553, 11.28727345, 9.321710063, 13.79725161, 10.37946697, 10.44739286, 8.584444837, 9.610370799, $
    9.185237408, 11.4532419, 12.92376277, 8.965779478, 9.449292694, 12.9895643, 9.003080504, 10.91140314, 10.46906395, $
    13.60804949, 11.19081041, 9.392709726, 7.365877697, 7.55138104, 7.950669154, 9.07586394, 11.45099755, 10.94095594, $
    8.846982379, 10.44887699, 10.37495512, 10.34772988, 9.221166399, 10.70580535, 9.801835338, 13.50916528, 8.112697761, $
    12.53105988, 12.68082997, 9.201920341, 9.818346521, 7.305895211, 9.431096077, 9.740691064, 10.78733361, 10.06159114, $
    9.849801695, 11.63517649, 10.60534006, 10.52340393, 10.58880114, 10.60669723, 11.12788664, 5.532796504, 9.476261694, 8.739105769]
  se3 = [6.052568983, 10.86044755, 9.567418683, 8.918799314, 9.133796573, 9.567593514, 8.190016132, 8.420028922, $
    9.575129585, 9.026907369, 11.83880003, 7.161704895, 9.895197001, 12.24989994, 5.549122382, 7.596327024, 4.963147115, $
    5.72490299, 7.259072724, 6.837871306, 12.4673407, 10.4785661, 8.162086601, 8.1347846, 10.4520931, 10.44483982, $
    9.851223114, 7.103385289, 13.81794881, 8.668484057, 13.17402494, 9.031713583, 9.500176523, 6.224464334, 11.71198054, $
    9.72389264, 9.178682409, 6.856358497, 5.582488908, 10.80968504, 12.93835786, 11.56657391, 11.99785626, 7.4004, $
    9.959178482, 10.53223687, 9.206125692, 10.14141876, 8.708431969] 
  se4 = [12.37192269, 12.21055248, 12.81697984, 13.85194539, 11.44937704, 12.59086298, 12.42183192, 11.34641009, $
    14.34135275, 13.20754664, 10.53675, 12.85003002, 11.67674845, 12.67590642, 11.37880378, 13.75501141, 11.14288889, $
    12.29904112, 12.62429318, 12.06027746, 11.55377904, 11.24682705, 11.45892799, 11.76679797, 12.96990266, 13.76885782, $
    12.56889501, 12.20249531, 11.46826305, 14.59787732, 12.43857216, 11.63659842, 12.46235768, 12.47022804, 14.48076324, $
    11.82461689, 13.04749525, 12.688293, 12.19963182, 12.87187203, 11.98798397, 12.03778334, 13.12093606, 10.50864505, $
    13.07185012, 12.5694805, 11.35830982, 13.01856936, 11.73495199, 12.61644109, 12.4003146, 12.49163577, 13.25138933, $
    13.99224292, 11.91047875, 11.66270048, 11.60908637, 11.81068681, 13.2976006, 12.02134545, 11.61358667, 13.74848521, $
    12.1968031, 12.38895482, 11.99551228, 12.07603551, 12.14599654, 11.19917555, 12.99476646, 10.43437502, 12.68246117, $
    11.74841019, 12.2451059, 12.62455172, 12.43729593, 12.25082481, 12.05694672, 13.13167516, 11.77902174, 11.63795444, $
    12.26304731, 12.09197321, 13.26753073, 11.52504354, 12.89201702, 12.68233548, 14.05608894, 13.25132519, 11.61840587, $
    12.96536655, 11.26935823, 11.65700965, 11.6916471, 11.79924005, 14.53792301, 12.02436976, 11.56382435, 12.19492223, $
    10.99718533, 11.19343751, 13.20566231, 11.97888129, 11.8080747, 10.47036721, 12.55864379, 13.20515112, 11.95680649, $
    11.14306813, 13.70014896, 11.02735091, 11.62707098, 12.40989708, 12.00535829, 12.45162933, 12.53695586, 11.77629423, $
    14.0630239, 11.48512533, 11.66296224, 13.76369555, 11.36460613, 12.93630833, 11.57714652, 12.71619901, 11.95079066, $
    13.79263884, 12.38629993, 13.00124059, 12.16591691, 12.77292863, 11.58111855, 11.84763479, 12.4847208, 12.2967661, $
    11.89295383, 11.99945066, 13.15567284, 11.02876752, 11.84058031, 12.65520845, 13.17106208, 12.18820995, 12.65521511, $
    14.38692032, 11.83167571, 14.07632536]
  se5 = [11.2592172, 12.60324193, 12.11383555, 13.15887155, 13.98726358, 12.99778887, 14.25690357, 10.35468732, 11.69400339, $
    12.51284284, 12.92090461, 13.23513706, 13.11730152, 14.27809663, 11.72941162, 10.73546104, 13.26812796, 13.18628342, $
    12.98599475, 12.03670591, 13.57195556, 12.51610005, 13.47892151, 12.10701507, 12.35993467, 11.18827463, 12.53555599, $
    14.76492986, 12.25454007, 12.20357258, 13.51204868, 11.71567842, 10.77369831, 13.04410794, 12.68339433, 12.0246233, $
    11.83625054, 11.93931386, 13.64417724, 11.12535747, 13.37452748, 13.03231892, 12.51767505, 12.91370138, 12.48938276, $
    13.08234182, 12.36025039, 11.37352562, 13.26482443, 13.3043938, 13.19240622, 14.59870898, 11.79926669, 12.45671158, $
    11.98533441, 12.28126272, 13.02684456, 13.90229125, 14.49525719, 14.80935394, 11.59293475, 14.71202815, 12.34005003, $
    13.13357175, 12.89731615, 12.98258052, 13.72819197, 11.15418173, 12.26669153, 11.49452115, 10.8324543, 11.7739252, $
    12.62499535, 12.66697641, 13.76840709, 14.4345553, 12.6548428, 12.78967321, 12.39809882, 13.01447026, 12.34316562, $
    12.91759309, 13.10693378, 11.90368477, 12.8230397, 13.08548303, 12.63507048, 12.70389186, 12.91303616, 12.39023407, $
    11.39388588, 12.33150326, 13.37150267, 11.77360992, 12.18183062, 13.46248064, 13.04752158, 12.90302086, 12.14340173, $
    12.15947235, 10.14054533, 10.87665924, 9.92249324, 12.95557655, 12.2061252, 12.30878464, 10.0084217, 12.61724886, $
    12.21911764, 12.34693749, 12.64457421, 14.12144976, 12.64357218]
  se6 = [15.08737925, 13.95096966, 12.38014606, 10.41002858, 11.07723349, 13.67861787, 14.08710531, 13.60801151, $
    13.87532508, 12.00954438, 14.60391721, 13.65426717, 15.38470589, 13.46592988, 12.65203809, 13.07333658, 13.4520145, $
    14.95366073, 13.88503218, 11.21831102, 15.00136686, 15.27191386, 14.84898474, 10.93267785, 11.31996286, 11.89343449, $
    15.54526444, 12.67729347, 13.8889, 14.5457818, 9.688794319, 11.62953132, 15.49492134, 15.20074172, 15.34166379, $
    14.70147642, 13.99503489, 13.81669535, 12.95074319, 11.86163305, 10.67042709, 14.87495466, 14.98584019, 13.74896005, $
    13.64594159, 14.94561081, 14.35064757, 14.27069676, 13.8423309, 12.05793288, 14.9829107, 14.90103055, 11.87217322, $
    13.54226117, 12.12257292, 14.57439406, 15.03139981, 13.89028067, 14.93188645, 14.54255304, 14.86151871, 13.31007296, $
    13.71480515, 14.17290396, 14.98942971, 13.82463093, 14.43098296, 14.80573889, 13.85481098, 13.55533459, 14.05212873, $
    14.03023824, 12.51902037, 13.94150858, 13.50634558, 12.48607896, 15.05521728, 13.32616636, 13.98334853, 11.77973214, $
    15.0675884, 15.4581895, 14.52433203, 11.46000959, 14.06723123, 15.2243597, 16.2475156, 13.2935449, 14.84104035, $
    15.49742532, 14.60384717, 13.47462324, 12.46065571, 14.3118632, 15.30079902, 15.46158628, 10.67471861, 12.76068504, $
    13.93443119, 12.59708338, 14.26348668, 13.22806937, 15.72505688, 14.47509666, 14.50000177]
  eae1 = [11.22583494, 11.51866456, 10.77599661, 12.60233506, 12.51088742, 11.19797218, 10.95524113, 10.94110816, $
    11.80345719, 11.25301325, 11.18151183, 9.063598747, 9.673037008, 11.57954533, 12.22807853, 8.418787469, 10.44798228, $
    9.01191154, 10.24929949, 12.12522268, 12.20339332, 11.42269653, 12.44828486, 12.47305981, 11.34744183, 11.69846598, $
    12.60907166, 10.68879022, 11.36765725, 11.02527231, 12.20202544, 11.28159371, 11.88327889, 11.57973252, 11.75596444, $
    11.44023875, 10.79650898, 11.04800565, 11.59562918, 10.78841116, 11.74074552, 9.964362523, 11.26881597, 12.98799212, $
    10.32492891, 11.21367023, 10.19489452, 11.31560648, 10.49486422, 11.84173112, 9.935062255, 10.83292943, 11.9780927, $
    11.20154522, 12.5203307, 11.61434146, 11.59337157]
  eae2 = [6.061391365, 13.79762383, 12.19364726, 11.8734865, 10.78768733, 12.56803955, 12.88222396, 11.00018484, $
    12.77100516, 11.6215188, 13.70317992, 11.74682561, 8.853468365, 8.247318552, 11.91958674, 13.05511534, 12.01399181, $
    9.33767639, 10.67107825, 11.98529982, 12.01985989, 12.28298273, 11.68971675, 9.291812689, 11.87627723, 8.545440212, $
    10.41094649, 11.12842515, 13.73882222, 10.91315778, 13.24169587, 11.45464906, 10.84476912, 11.85959104, 11.74525348, $
    11.22983428, 12.55255214, 12.38439013, 11.65265038, 12.78260453, 10.10304553, 11.2681242, 13.65798755, 12.57763266, $
    11.57523042, 10.92255935, 10.55220879, 12.15013667, 12.55043075, 11.40132927, 11.42024611, 12.74355081, 12.72193733, $
    10.82396603, 13.15735366, 13.13883809, 13.10583249, 11.76803788, 12.2515114, 11.92181329, 11.28990096, 7.565319045, $
    11.77678503, 11.23173572, 11.66276751, 13.20542187, 13.38786214, 13.24294475, 13.43323268, 12.07997399, 12.66987635, $
    11.5647382, 9.728396178, 11.96419114, 12.56114677, 12.72649206, 11.71707113, 12.56230761, 11.53021559, 12.53355311, $
    11.26192827, 12.50109499, 12.52267537, 12.81759805, 11.6140809, 11.26077881, 10.46183953, 12.59965851, 12.50240938, $
    12.32949005, 11.57912155, 12.92653799, 12.52803727, 12.72883171]
  eae3 = [13.68391492, 12.31453767, 13.88340572, 13.48578904, 15.02741441, 11.86613786, 12.00259266, 14.16100181, $
    14.4664123, 14.26698079, 12.16839083, 13.56812582, 14.50424675, 13.82387473, 14.04195707, 14.3853138, 14.26225235, $
    13.33986783, 13.33779958, 14.57018532, 15.18713937, 13.69100952, 14.08477625, 12.95321913, 12.80999399, 13.30795124, $
    13.13679459, 11.9306957, 15.06801663, 13.80190357, 12.61888517, 13.96585351, 14.13062401, 13.09389022, 12.75295561, $
    14.14517812, 13.91287876, 13.7644559, 12.80994779, 12.59299149, 12.21231274, 15.1363208, 13.22699271, 13.17180188, $
    10.84884673, 13.42796719, 13.0812527, 15.57486031, 12.91139774, 15.01392298, 12.51927017, 13.754318, 13.7028854, $
    13.62980825, 13.37092269, 14.66474886, 12.40271183, 12.24620934, 15.3170726, 15.0496707, 11.35586185, 12.98987168, $
    12.2854746, 13.98439406, 12.16461219, 14.27062408, 12.11003086, 13.40408537, 13.41070837, 15.51652648, 13.5497785, $
    13.74591609, 13.73327643, 13.03829969, 13.13691655, 11.49140133, 11.04820912, 10.1597526, 12.82177393, 14.10286453, $
    11.23579672]
End


Pro TestVTx
; 0.022, 1.64, 0.0, 8.48
; 0.022, 1.25, 4.874, 6.028
  vB = 0.022
  vTmax = 1.25
  delXTo = 4.874
  delXTo_B = 6.028
  
  tsBegin = 2.
  tsEnd = 12.
  tsStep = 2.
  numSteps = Round(abs((tsEnd-tsBegin)/tsStep))+1
  
  numPts = 100
  trueLen = delXTo + 2.*delXTo_B
  xInts = trueLen * IndGen(numPts)/(numPts-1.)
  
  lengths = FltArr(numSteps,numPts)
  vtvbs = lengths
  
  for step = 0, numSteps-1 do begin
    ts = tsBegin + step*tsStep
    for i=0,N_Elements(xInts)-1 do begin
      halfLen1 = VTX_ConstantCore(20., xInts[i], ts, vB, vTmax, delXTo, delXTo_B, END_VTVB=vtvb1)
      halfLen2 = VTX_ConstantCore(20., trueLen - xInts[i], ts, vB, vTmax, delXTo, delXTo_B, END_VTVB=vtvb2)
      lengths[step,i] = halfLen1 + halfLen2
      vtvbs[step,i] = vtvb1 > vtvb2
    endfor
    print, Moment(lengths[step,*]), Max(lengths[step,*]), Min(lengths[step,*])
  endfor
  lenPlot = plot(xInts, lengths[0,*], XTITLE="Impingement Point (µm)", YTITLE="Length (µm)")
  for i=0,numSteps-1 do ndPlot = plot(xInts,lengths[i,*],/OVERPLOT)

  vtvbPlot = plot(xInts, vtvbs[0,*], XTITLE="Impingement Point (µm)", YTITLE="Maximum tip vT/vB")
  for i=0,numSteps-1 do vtvbPlot = plot(xInts,vtvbs[i,*],/OVERPLOT)

  for i=0,N_Elements(xInts)-1 do print, [xInts[i],lengths[*,i],vtvbs[*,i]]

End
