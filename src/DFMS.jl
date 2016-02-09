module DFMS
using DataFrames
using LsqFit
using PyPlot
using HDF5
using Spice


export crop,
       clean,
       clean_h5,
       findPeaks,
       getGainFactor,
       getPixelGain,
       get_eph,
       get_ephemeris_spice,
       get_outliers,
       load_h5_data,
       load_h5_data!,
       load_pds,
       model,
       parseDataFileBothRows,
       parseDataFile,
       pds2h5,
       peakFit,
       polyFit,
       remove_outliers,
       remove_outliers!,
       show_me,
       show_time_series,
       show_single_spec,
       single_spec,
       time_series,
       time_series_autospec,
       update!,
       findSimilar,
       doCO2Correction

c_adc = 2.5 / (2.^12 - 1)
c_leda = 4.22 * 10.0^-12
q = 1.60217657 * 10.0^-19
const global t_integration = 19.8
const global C = c_adc * c_leda / q

const global _rowA = 1
const global _rowB = 2

type Spectrum
  gainStep::Int
  tStart::DateTime
  commandedMass::Int
  p_pt::Int
  row_A::Vector{Int}
  row_B::Vector{Int}
  pixelGain_A::Vector{Float64}
  pixelGain_B::Vector{Float64}
  baseline_A::Vector{Float64}
  baseline_B::Vector{Float64}
  ionsPerSpectrum_A::Vector{Float64}
  ionsPerSpectrum_B::Vector{Float64}
  peakIndices_A::Vector{Int}
  peakIndices_B::Vector{Int}
  peakAmplitudes_A::Vector{Float64}
  peakAmplitudes_B::Vector{Float64}
end

###############################################################################
# functions to load data
###############################################################################
function time_series(fileName)
  label = basename(fileName)
  df = readtable(fileName)
  df[:date] = DateTime[DateTime(tStr) for tStr in df[:date]];
  return df;
end

function time_series_autospec(fileName)
  label = basename(fileName)
  df = readtable(fileName, separator='\t')
  df[:date] = DateTime[DateTime(tStr, "yyyy-mm-dd HH:MM:SS") for tStr in df[:Time]];
  return df
end

function single_spec(fileName; dataOnly=false)
  y, gainStep, t, m0, p_pt = parseDataFileBothRows(fileName)
  row_A = vec(y[:,1])
  row_B = vec(y[:,2])
  gainFactor_A = getGainFactor(gainStep, _rowA)
  gainFactor_B = getGainFactor(gainStep, _rowB)
  pixelGain_A = getPixelGain(t, _rowA, gainStep)
  pixelGain_B = getPixelGain(t, _rowB, gainStep)
  baseline_A = polyFit(row_A, model)
  baseline_B = polyFit(row_B, model)

  ionsPerSpectrum_A = (row_A .- baseline_A) ./ pixelGain_A / gainFactor_A * C
  ionsPerSpectrum_B = (row_B .- baseline_B) ./ pixelGain_B / gainFactor_A * C
  yMin_A = 4 * std(ionsPerSpectrum_A[80:180])
  yMin_B = 4 * std(ionsPerSpectrum_B[80:180])


  @show(t)
  @show(gainFactor_A)
  @show(gainFactor_B)
  @show(maximum(ionsPerSpectrum_A))
  @show(yMin_A)

  pI_A, pA_A = findPeaks(ionsPerSpectrum_A, 4, yMin_A)
  pI_B, pA_B = findPeaks(ionsPerSpectrum_B, 4, yMin_B)

  s = Spectrum(gainStep, t, m0, p_pt, row_A, row_B, pixelGain_A,
               pixelGain_B, baseline_A, baseline_B, ionsPerSpectrum_A,
               ionsPerSpectrum_B, pI_A, pI_B, pA_A, pA_B)
  return s
end

################################################################################
# functions for plotting data
################################################################################
function show_me(df::DataFrame, var=:peakArea; logy=true)
  figure()
  if logy
    semilogy(df[:date], df[var], "ok")
  else
    plot(df[:date], df[var], "ok")
  end
  grid(true)
  ylabel(string(var))
end

function show_me(dfs::Vector{DataFrame}, var=:peakArea; logy=true)
  figure()
  for df in dfs
    if logy
      semilogy(df[:date], df[var], "o", lw=2)
    else
      plot(df[:date], df[var], "o", lw=2)
    end
  end
  ylabel(string(var))
  grid(true)
end

function show_me(s::Spectrum; logy=true)
  if logy
    figure()
    semilogy(s.baseline_A, "ok", markerfacecolor="white")
    semilogy(s.baseline_B, "or", markerfacecolor="white")
    semilogy(s.row_A, "-k", lw=2, label="Row A")
    semilogy(s.row_B, "-r", lw=2, label="Row B")
    legend()

    figure()
    semilogy(s.ionsPerSpectrum_A, lw=2, label="ions per spectrum A")
    semilogy(s.ionsPerSpectrum_B, lw=2, label="ions per spectrum B")
    semilogy(s.peakIndices_A, s.peakAmplitudes_A, "or")
    semilogy(s.peakIndices_B, s.peakAmplitudes_B, "sr")
  else
    figure()
    plot(s.baseline_A, "ok", markerfacecolor="white")
    plot(s.baseline_B, "or", markerfacecolor="white")
    plot(s.row_A, "-k", lw=2, label="Row A")
    plot(s.row_B, "-r", lw=2, label="Row B")
    legend()

    figure()
    plot(s.ionsPerSpectrum_A, lw=2, label="ions per spectrum A")
    plot(s.ionsPerSpectrum_B, lw=2, label="ions per spectrum B")
    plot(s.peakIndices_A, s.peakAmplitudes_A, "or")
    plot(s.peakIndices_B, s.peakAmplitudes_B, "sr")


  end
  grid(true)
  legend()
end

function show_me(fileName::ASCIIString; logy=false)
  y, gainStep, t, m0, p_pt = parseDataFileBothRows(fileName)
  @show(gainStep)
  @show(t)
  @show(m0)
  @show(p_pt)
  figure()
  if logy
    semilogy(y[:,1], "-ok", label="Row A")
    semilogy(y[:,2], "-sr", label="Row B")
  else
    plot(y[:,1], "-ok", label="Row A")
    plot(y[:,2], "-sr", label="Row B")
  end
  grid(true)
  legend()
end


function show_me(y::Vector{Real}; logy=true)
  figure()
  if logy
    semilogy(y, "-k", lw=2)
  else
    plot(y, "-k", lw=2)
  end
  grid(true)
  nothing
end


################################################################################
# functions to modify loaded time series
################################################################################

function crop!(df, tStart::DateTime, tStop::DateTime)
  flt1 = df[:date] .>= tStart
  flt2 = df[:date] .<= tStop
  flt = flt1 & flt2
  df = df[flt, :]
end

function crop!(df, tStart::ASCIIString, tStop::ASCIIString)
  crop!(df, DateTime(tStart), DateTime(tStop))
end

function crop(df, tStart::DateTime, tStop::DateTime)
  flt1 = df[:date] .>= tStart
  flt2 = df[:date] .<= tStop
  flt = flt1 & flt2
  return df[flt, :]
end

function crop(df, tStart::ASCIIString, tStop::ASCIIString)
  crop(df, DateTime(tStart), DateTime(tStop))
end

function clean!(df)
  flt = ones(Bool, size(df, 1))
  for var in [:peakArea, :peakCenter, :peakWidth]
    flt &= get_outliers(df, var=var)
  end
  df = df[flt, :]
end

function clean(df)
  flt = ones(Bool, size(df, 1))
  for var in [:peakArea, :peakCenter, :peakWidth]
    flt &= get_outliers(df, var=var)
  end
  try
    flt_offNadir = df[:offNadirAngle_deg] .<= 25.0
    flt &= flt_offNadir
  end
  return df[flt, :]
end

function remove_outliers(df; N=50, var=:peakArea)
  flt = get_outliers(df, N=N, var=var)
  return df[flt,:]
end

function remove_outliers!(df; N=50, var=:peakArea)
  flt = get_outliers(df, N=N, var=var)
  df = df[flt,:]
end

function get_outliers(df; N=50, var=:peakArea)
  y = convert(Array, df[var])
  nDataPoints = length(y)
  flt = ones(Bool, nDataPoints)
  for i = N+1:nDataPoints-N
    data = append!(y[i-N:i-1], y[i+1:i+N])
    if abs(median(data) - y[i]) > (3*std(data))
      flt[i] = false
    end
  end
  return flt
end

################################################################################
# functions for spice
################################################################################
function get_ephemeris_spice(t::ASCIIString)
  get_ephemeris_spice([DateTime(t)])
end

function get_ephemeris_spice(t::DateTime)
  get_ephemeris_spice([t])
end

function get_eph(tt::Vector{DateTime})
  gc_enable(false)
  furnsh("/home/abieler/rosetta/spiceKernels/metafiles/operationalKernels.tm")
  for t in tt
    tStr = string(t)
    et = utc2et(tStr)
  end
  gc_enable(true)
end

function get_ephemeris_spice(tt::Vector{DateTime})
  furnsh("/home/abieler/rosetta/spiceKernels/metafiles/operationalKernels.tm")
  N = length(tt)
  objectName = "ROSETTA"
  CG_ID = 1000012
  frame = "67P/C-G_CK"
  observer = "CHURYUMOV-GERASIMENKO"
  km2AU = 6.68458712 * 1e-9
  rotAxis = [0.0,0.0,1.0]
  pCenter = zeros(Float64, 3)
  lon = zeros(Float64, N)
  lat = zeros(Float64, N)
  localHourFloat = zeros(Float64, N)
  xSC = zeros(Float64, N)
  ySC = zeros(Float64, N)
  zSC = zeros(Float64, N)
  xSUN = zeros(Float64, N)
  ySUN = zeros(Float64, N)
  zSUN = zeros(Float64, N)
  lonSunArr = zeros(Float64, N)
  latSunArr = zeros(Float64, N)
  lonSC_cso = zeros(Float64, N)
  latSC_cso = zeros(Float64, N)
  heliocentricDistance = zeros(Float64, N)
  heliocentricDistance_AU = zeros(Float64, N)
  cometocentricDistance = zeros(Float64, N)
  phaseAngle = zeros(Float64, N)
  offNadirAngle = zeros(Float64, N)
  et = 0.0
  gc_enable(false)
  for (i,t) in enumerate(tt)
    timeStamp = string(t)
    et = utc2et(timeStamp)
    r_SC, lt = spkpos("ROSETTA", et, "67P/C-G_CK", "NONE", observer)
    r_SUN, lt = spkpos("SUN", et, "67P/C-G_CK", "NONE", observer)
    r_CG, lt = spkpos("67P/C-G", et, "ROS_VIRTIS-H", "NONE", "ROSETTA")
    r_SC_cso, lt = spkpos("ROSETTA", et, "67P/C-G_CSO", "NONE", observer)

    cometocentricDistance[i] = norm(r_SC)
    heliocentricDistance[i] = norm(r_SUN)
    heliocentricDistance_AU[i] = heliocentricDistance[i] * km2AU

    r_SC_hat = r_SC ./ cometocentricDistance[i]
    r_SUN_hat = r_SUN ./ heliocentricDistance[i]

    pA = vsep(r_SC, r_SUN)
    oNA = vsep(r_CG, rotAxis)
    phaseAngle[i] = pA / pi * 180.0
    offNadirAngle[i] = oNA / pi * 180.0

    xSC[i] = r_SC[1]
    ySC[i] = r_SC[2]
    zSC[i] = r_SC[3]
    xSUN[i] = r_SUN[1]
    ySUN[i] = r_SUN[2]
    zSUN[i] = r_SUN[3]

    r, llon, llat = reclat(r_SC)
    llon = llon / pi * 180.0
    llat = llat / pi * 180.0
    lon[i] = llon
    lat[i] = llat

    rSun, lonSun, latSun = reclat(r_SUN)
    lonSun = lonSun / pi * 180
    latSun = latSun / pi * 180
    lonSunArr[i] = lonSun
    latSunArr[i] = latSun

    rr, lon_cso, lat_cso = reclat(r_SC_cso)
    lon_cso = lon_cso / pi * 180
    lat_cso = lat_cso / pi * 180
    lonSC_cso[i] = lon_cso
    latSC_cso[i] = lat_cso


    deltaLongitude = (((llon - lonSun) + 180) % 360) - 180
    localHour = 12 + (deltaLongitude / 180.0 * 12.0)
    localHourFloat[i] = localHour
    localMin = round(Int, floor(localHour % 1 * 60))
    localSec = round(Int, (localHour % 1 * 60) % 1 * 60)
    localHour = round(Int, floor(localHour))
  end
  gc_enable(true)

  df = DataFrame()
  df[:date] = tt
  df[:localHourFloat] = localHourFloat
  df[:rSC_km] = cometocentricDistance
  df[:rSUN_km] = heliocentricDistance
  df[:rSUN_AU] = heliocentricDistance_AU
  df[:lonSC_deg] = lon
  df[:latSC_deg] = lat
  df[:lonSUN_deg] = lonSunArr
  df[:latSUN_deg] = latSunArr
  df[:lonSC_cso_deg] = lonSC_cso
  df[:latSC_cso_deg] = latSC_cso
  df[:phaseAngle_deg] = phaseAngle
  df[:offNadirAngle_deg] = offNadirAngle
  df[:xSC_km] = xSC
  df[:ySC_km] = ySC
  df[:zSC_km] = zSC
  df[:xSUN_km] = xSUN
  df[:ySUN_km] = ySUN
  df[:zSUN_km] = zSUN

  return df
end


function pds2h5(path)
  homeDir = "/home/abieler/rosetta/data/dfms/h5/"
  fileNames = readdir(path)

  for fileName in sort(fileNames)
    # remove the .TAB of fileName
    bName = fileName[1:end-4]
    yrStr, tStr, modeStr = matchall(r"(\d+)", bName)
    if (modeStr[end] == '2') & (modeStr != "0572") & (modeStr[1:2] != "06")
      h5Str = joinpath(homeDir, yrStr[1:6] * ".h5")
      y, gainStep, t, m0, p_pt = parseDataFileBothRows(joinpath(path, fileName))
      mstr = "m" * string(round(Int,m0))
      dataSetName = mstr * "/" * bName
      try
        h5write(h5Str, joinpath(dataSetName, "gainStep"), gainStep)
        h5write(h5Str, joinpath(dataSetName, "commandedMass"), m0)
        h5write(h5Str, joinpath(dataSetName, "p_pt"), p_pt)
        h5write(h5Str, joinpath(dataSetName, "y"), y)
        h5write(h5Str, joinpath(dataSetName, "tStart"), string(t))
      catch
      end
    end
  end
end


function pds2h5(path)
  homeDir = "/home/abieler/rosetta/data/dfms/h5/"
  fileNames = readdir(path)

  for fileName in sort(fileNames)
    # remove the .TAB of fileName
    bName = fileName[1:end-4]
    yrStr, tStr, modeStr = matchall(r"(\d+)", bName)
    if (modeStr[end] == '2') & (modeStr != "0572") & (modeStr[1:2] != "06")
      h5Str = joinpath(homeDir, yrStr[1:6] * ".h5")
      y, gainStep, t, m0, p_pt = parseDataFileBothRows(joinpath(path, fileName))
      mstr = "m" * string(round(Int,m0))
      dataSetName = mstr * "/" * bName
      try
        h5write(h5Str, joinpath(dataSetName, "gainStep"), gainStep)
        h5write(h5Str, joinpath(dataSetName, "commandedMass"), m0)
        h5write(h5Str, joinpath(dataSetName, "p_pt"), p_pt)
        h5write(h5Str, joinpath(dataSetName, "y"), y)
        h5write(h5Str, joinpath(dataSetName, "tStart"), string(t))
      catch
      end
    end
  end
end

################################################################################
# helper functions for the automatic DFMS data analysis of HR spectra
################################################################################
function getGainFactor(gainStep,row)
  fileName = ""
  gainFactor = 0.0
  if (row == _rowA)
    fileName = "gainFactor_SPACE_A.txt"
  elseif (row == _rowB)
    fileName = "gainFactor_SPACE_B.txt"
  end

  iFile = open(joinpath("/home/abieler/.julia/v0.4/DFMS/InstrumentData/", fileName), "r")
  i=0
  while !eof(iFile)
    line = readline(iFile)
    if (i == gainStep)
      gainFactor = parse(Float64, matchall(r"(\d+.\d+[eE][+-]\d+)", line)[1])
    end
    i += 1
  end
  close(iFile)
  return gainFactor
end

function getGainFactor(gainStep)
  fileName = ""
  gainFactor = [0.0, 0.0]
  counter = 1
  gF = 0.0
  for myrow in [_rowA, _rowB]
    if (myrow == _rowA)
      fileName = "gainFactor_SPACE_A.txt"
    elseif (myrow == _rowB)
      fileName = "gainFactor_SPACE_B.txt"
    end

    iFile = open(joinpath("/home/abieler/.julia/v0.4/DFMS/InstrumentData/", fileName), "r")
    i=0
    while !eof(iFile)
      line = readline(iFile)
      if (i == gainStep)
        gF = parse(Float64, matchall(r"(\d+.\d+[eE][+-]\d+)", line)[1])
      end
      i += 1
    end
    close(iFile)
    gainFactor[counter] = gF
    counter += 1
  end
  return gainFactor
end
function parseVar(fileName, varName, iSkip=329)
  iFile = open(fileName, "r")
  t = DateTime(2000,1,1)
  isStartTimeFound = false
  while !eof(iFile)
    line = readline(iFile)
    if (contains(line, "START_TIME") & (isStartTimeFound == false))
      tStr = matchall(r"\d+-\d+-\d+T\d+:\d+:\d+", line)[1]
      t = DateTime(tStr)
      isStartTimeFound = true
    elseif contains(line, varName)
      value = parse(Float64, matchall(r"(-?\d.\d+[eE][+-]\d+)", line)[1])
    end

    if (i >= iSkip)
      break
    end
  end
  close(iFile)
  return t, value
end

function parseDataFile(fileName, row=2)
  iFile = open(fileName, "r")
  i = 1
  y = zeros(Float64, 512)
  t = DateTime(2000)
  gainStep = 0
  m0 = 0
  p_pt = 0.0
  iSkip = 329
  isStartTimeFound = false
  while !eof(iFile)
    line = readline(iFile)
    if (contains(line, "START_TIME") & (isStartTimeFound == false))
      tStr = matchall(r"\d+-\d+-\d+T\d+:\d+:\d+", line)[1]
      t = DateTime(tStr)
      isStartTimeFound = true
    elseif contains(line, "ROSINA_DFMS_SCI_GAIN_OF_SPECTRUM")
      gainStep = parse(Int, matchall(r"(\d+)", line)[1])
    elseif contains(line , "ROSINA_DFMS_SCI_MASS")
      m0 = parse(Int, matchall(r"(\d+)", line)[1])
    elseif contains(line, "ROSINA_DFMS_SCI_P_PT_AT_PEAK")
      p_pt = parse(Int, matchall(r"(\d+)", line)[1])
    elseif contains(line, "SPICE_FILE_NAME")
      iSkip = 341
    end

    if (i >= iSkip)
      y[i-(iSkip-1)] = parse(Float64, matchall(r"(\d+)", line)[row])
    end
    i += 1
  end

  close(iFile)
  return gainStep, t, y, m0, p_pt
end

function load_pds(fileName)
  iFile = open(fileName, "r")
  i = 1
  y = zeros(Float64, 512, 2)
  t = DateTime(2000)
  gainStep = 0
  m0 = 0
  p_pt = 0.0
  iSkip = 329
  isStartTimeFound = false
  while !eof(iFile)
    line = readline(iFile)
    if (contains(line, "START_TIME") & (isStartTimeFound == false))
      tStr = matchall(r"\d+-\d+-\d+T\d+:\d+:\d+", line)[1]
      t = DateTime(tStr)
      isStartTimeFound = true
    elseif contains(line, "ROSINA_DFMS_SCI_GAIN_OF_SPECTRUM")
      gainStep = parse(Int, matchall(r"(\d+)", line)[1])
    elseif contains(line , "ROSINA_DFMS_SCI_MASS")
      m0 = parse(Int, matchall(r"(\d+)", line)[1])
    elseif contains(line, "ROSINA_DFMS_SCI_P_PT_AT_PEAK")
      p_pt = parse(Int, matchall(r"(\d+)", line)[1])
    elseif contains(line, "SPICE_FILE_NAME")
      iSkip = 341
    end

    if (i >= iSkip)
      y[i-(iSkip-1), 1] = parse(Float64, matchall(r"(\d+)", line)[2])
      y[i-(iSkip-1), 2] = parse(Float64, matchall(r"(\d+)", line)[3])
    end
    i += 1
  end

  close(iFile)
  return y, gainStep, t, m0, p_pt
end

function load_h5_data(dataset, y, row)
    gainStep = read(dataset, "gainStep")
    m0 = read(dataset, "commandedMass")
    p_pt = read(dataset, "p_pt")
    rows = read(dataset, "y")
    tStr = read(dataset, "tStart")
    t = DateTime(tStr)
    for k in 1:512
      y[k] = rows[k, row]
    end
    return gainStep, t, m0, p_pt
end

function update!(dataset, varName, var)
  if !has(dataset, varName)
    d_write(dataset, varName, var)
  else
    dataset[varName][:,:] = var
  end
end

function load_h5_data!(dataset, y)
    gainStep = read(dataset, "gainStep")
    m0 = read(dataset, "commandedMass")
    p_pt = read(dataset, "p_pt")
    rows = read(dataset, "y")
    tStr = read(dataset, "tStart")
    t = DateTime(tStr)
    for myRow = 1:2
      for pixel in 1:512
        y[pixel, myRow] = rows[pixel, myRow]
      end
    end
    return gainStep, t, m0, p_pt
end


function parseDataFileBothRows(fileName)
  iFile = open(fileName, "r")
  i = 1
  y = zeros(Int, 512, 2)
  t = DateTime(2000)
  gainStep = 0
  m0 = 0
  p_pt = 0.0
  iSkip = 329
  isStartTimeFound = false
  while !eof(iFile)
    line = readline(iFile)
    if (contains(line, "START_TIME") & (isStartTimeFound == false))
      tStr = matchall(r"\d+-\d+-\d+T\d+:\d+:\d+", line)[1]
      t = DateTime(tStr)
      isStartTimeFound = true
    elseif contains(line, "ROSINA_DFMS_SCI_GAIN_OF_SPECTRUM")
      gainStep = parse(Int, matchall(r"(\d+)", line)[1])
    elseif contains(line , "ROSINA_DFMS_SCI_MASS")
      m0 = parse(Int, matchall(r"(\d+)", line)[1])
    elseif contains(line, "ROSINA_DFMS_SCI_P_PT_AT_PEAK")
      p_pt = parse(Int, matchall(r"(\d+)", line)[1])
    elseif contains(line, "SPICE_FILE_NAME")
      iSkip = 341
    end

    if (i >= iSkip)
      y[i-(iSkip-1), 1] = parse(Int, matchall(r"(\d+)", line)[2])
      y[i-(iSkip-1), 2] = parse(Int, matchall(r"(\d+)", line)[3])
    end
    i += 1
  end

  close(iFile)
  return y, gainStep, t, m0, p_pt
end

function getPixelGain(t, row, gainStep)
  fileName = ""
  i = 0
  path = "/home/abieler/.julia/v0.4/DFMS/InstrumentData/pixelgain_SPACE"
  # pixel gains for july to 1 october
  if DateTime(2014,2,1) < t < DateTime(2014,10,1)
    fileName = joinpath(path, "pg_space_Aug2014.csv")
    i=2
  # pixel gains from 1 october to 31 december
  elseif DateTime(2014,10,1) <= t < DateTime(2015,1,1)
    if gainStep < 12
      fileName = joinpath(path, "pg_space_Nov2014_GS9_new.csv")
    elseif 12 <= gainStep < 15
      fileName = joinpath(path, "pg_space_Nov2014_GS14_new.csv")
    elseif gainStep >= 15
      fileName = joinpath(path, "pg_space_Nov2014_GS16_new.csv")
    end
    i=1
  # pixel gains from 1 january to april
  elseif DateTime(2015,1,1) <= t < DateTime(2016,4,1)
    if gainStep < 12
      fileName = joinpath(path, "pg_space_Feb2015_GS9_new.csv")
    elseif 12 <= gainStep < 15
      fileName = joinpath(path, "pg_space_Feb2015_GS14_new.csv")
    elseif gainStep >= 15
      fileName = joinpath(path, "pg_space_Feb2015_GS16_new.csv")
    end
    i=1
  end
  iFile = open(fileName, "r")
  pixelGain = ones(Float64, 512)
  while !eof(iFile)
    line = readline(iFile)
    pixelGain[i] = parse(Float64, matchall(r"(\d+.\d+)", line)[row+2])
    i += 1
  end
  return pixelGain
end

function getPixelGain(t, gainStep)
  fileName = ""
  i = 0
  path = "/home/abieler/.julia/v0.4/DFMS/InstrumentData/pixelgain_SPACE"
  pixelGain = zeros(Float64, 512, 2)
  for row in [_rowA, _rowB]
    # pixel gains for july to 1 october
    if DateTime(2014,2,1) < t < DateTime(2014,10,1)
      fileName = joinpath(path, "pg_space_Aug2014.csv")
      i=2
    # pixel gains from 1 october to 31 december
    elseif DateTime(2014,10,1) <= t < DateTime(2015,1,1)
      if gainStep < 12
        fileName = joinpath(path, "pg_space_Nov2014_GS9_new.csv")
      elseif 12 <= gainStep < 15
        fileName = joinpath(path, "pg_space_Nov2014_GS14_new.csv")
      elseif gainStep >= 15
        fileName = joinpath(path, "pg_space_Nov2014_GS16_new.csv")
      end
      i=1
    # pixel gains from 1 january to april
    elseif DateTime(2015,1,1) <= t < DateTime(2016,4,1)
      if gainStep < 12
        fileName = joinpath(path, "pg_space_Feb2015_GS9_new.csv")
      elseif 12 <= gainStep < 15
        fileName = joinpath(path, "pg_space_Feb2015_GS14_new.csv")
      elseif gainStep >= 15
        fileName = joinpath(path, "pg_space_Feb2015_GS16_new.csv")
      end
      i=1
    end
    iFile = open(fileName, "r")
    while !eof(iFile)
      line = readline(iFile)
      pixelGain[i, row] = parse(Float64, matchall(r"(\d+.\d+)", line)[row+2])
      i += 1
    end
  end
  return pixelGain
end

function multiDoubleGauss(x, p)
  nPeaks = round(Int, (length(p)-1) / 5)
  B = p[end]
  yFit = zeros(Float64, length(x))
  for i=1:nPeaks
    x0 = p[i]
    A = p[i+nPeaks]
    w = p[i+2*nPeaks]
    A2 = abs(p[i+3*nPeaks])
    w2 = p[i+4*nPeaks]
    yFit += (A * exp(-(x-x0).^2 / w^2) + A2 * exp(-(x-x0).^2 / w2^2))
  end
  yFit += B
  return yFit
end

function multiGauss(x, p)
  # x = bins
  # p = fitting parameters
  # B = baseline height
  w0 = 4.0
  dw = 1.25

  nPeaks = round(Int, (length(p)-1) / 3)
  B = 1.0
  #B = p[end]
  yFit = zeros(Float64, length(x))
  for i=1:nPeaks
    x0 = p[i]
    A = p[i+nPeaks]
    w = p[i+2*nPeaks]
    th = tanh((w-w0)/dw)
    w_th = w0 + th * dw
    yFit += A * exp(-(x-x0).^2 / w_th^2)
  end
  yFit += B
  return yFit
end

function model(x, p)
  f = zeros(Float64, length(x))
  for i=1:length(x)
    @inbounds f[i] = p[1] + p[2]*x[i] + p[3]*x[i]^2 + p[4]*x[i]^3
  end
  return f
end

function peakFit(y, pI, pA, fitMethod)
  nPeaks = length(pI)
  peakArea = zeros(Float64, nPeaks)
  peakIndex = Float64[value for value in pI]
  peakWidth = zeros(Float64, nPeaks)

  if fitMethod == "sum"
    peakArea = sumPixels(y,[peakIndex])
    return y, peakArea, pI[iPeak], 0.0

  elseif fitMethod == "singleGauss"
    # use single gauss to fit peak shape
    fitParams = Float64[]
    append!(fitParams, pI)
    append!(fitParams, pA)
    append!(fitParams, ones(Float64, length(pI)))
    #push!(fitParams, median(y))

    xmin = 100
    xmax = 400
    x = collect(1:512)

    # compute weight factors w
    weights = zeros(Float64, 512)
    sigma = 12.0
    if nPeaks == 1
      sigma = 50.0
    end
    for i=1:length(pI)
      weights += exp(-0.5 * ((x - pI[i]) / sigma).^2) / pA[i]
    end

    fit = curve_fit(multiGauss, x[xmin:xmax], y[xmin:xmax], weights[xmin:xmax],
                    fitParams)

    for i=1:nPeaks
      peakArea[i] = fit.param[i+nPeaks] * abs(fit.param[i+2*nPeaks]) * sqrt(pi)
      peakIndex[i] = fit.param[i]
      peakWidth[i] = fit.param[i+2*nPeaks]
    end

    w0 = 4.0
    dw = 1.25
    th = tanh((peakWidth-w0)/dw)
    peakWidth = w0 + th * dw

    return multiGauss(x, fit.param), peakArea, peakIndex, peakWidth
  end
end

function peakFit(y, pI, pA, fitMethod, iPeak)
  nPeaks = length(pI)
  peakArea = 0.0
  peakIndex = pI[iPeak]
  peakWidth = 0.0

  if fitMethod == "sum"
    peakArea = sumPixels(y,[peakIndex])
    return y, peakArea, pI[iPeak], 0.0

  elseif fitMethod == "singleGauss"
    # use single gauss to fit peak shape
    fitParams = Float64[]
    append!(fitParams, pI)
    append!(fitParams, pA)
    append!(fitParams, ones(Float64, length(pI)))
    push!(fitParams, median(y))
    x = collect(1:512)
    fit = curve_fit(multiGauss, x[20:500], y[20:500], fitParams)
    #fit = curve_fit(multiGauss, x[LHS:RHS], y[LHS:RHS], fitParams)
    peakArea = fit.param[iPeak+nPeaks] * abs(fit.param[iPeak+2*nPeaks]) * sqrt(pi)
    peakIndex = fit.param[iPeak]
    peakWidth = fit.param[iPeak+2*nPeaks]
    return multiGauss(x, fit.param), peakArea, peakIndex, peakWidth
  end
end

function polyFit(y, model)
  # select only subset of the 512 data bins for fitting
  # of baseline signal
  # --> center of the mass spectrum is ignored as it
  # contains all the peaks
  bgRangeLeft = 10:180
  bgRangeRight = 380:500
  fullRange = 1:512
  xBG = collect(bgRangeLeft)
  yBG = y[bgRangeLeft]
  append!(xBG, collect(bgRangeRight))
  append!(yBG, y[bgRangeRight])

  fit = curve_fit(model, xBG, yBG, [median(y), 0., 0., 0.])

  return model(collect(fullRange), fit.param)
end

function sumPixels(y, pI)
  lhs = 238
  rhs = 285
  if length(pI) > 1
    lhs = minimum(pI)-20
    rhs = maximum(pI)+20
  elseif length(pI) == 1
    lhs = pI[1] - 20
    rhs = pI[1] + 20
  else
    return 0.0
  end
   totalArea = sum(y[lhs:rhs])
end

function findPeaks(y, NN=3, yMin=0.01, pkLHS=180, pkRHS=380)
  peakIndexes = Int64[]
  peakAmpl = Float64[]
  for i=pkLHS:pkRHS
    score = 0
    for j=i-NN:i+NN
      if ((y[i] >= y[j]) & (y[i] >= yMin))
        score += 1
      end
    end
    if (score == 2*NN+1)
      push!(peakIndexes, i-1)
      push!(peakAmpl, y[i])
    end
  end

  return peakIndexes, peakAmpl
end


function findSimilar(t1, n1, t2, n2)
  const dtMax = DateTime(2000,1,1,0,35) - DateTime(2000,1,1,0)
  t_sync = DateTime[]
  n1_sync = Float64[]
  n2_sync = Float64[]

  for i=1:length(t2)
    for j=1:length(t1)
      if ((t2[i] > t1[j]) & ((t2[i] - t1[j]) < dtMax))
        push!(n1_sync, n1[j])
        push!(n2_sync, n2[i])
        push!(t_sync, t2[i])
        break
      end
    end
  end
  return t_sync, n1_sync, n2_sync

end

function findSimilar(df1::DataFrame, df2::DataFrame)
  const dtMax = DateTime(2000,1,1,0,35) - DateTime(2000,1,1,0)
  t_sync = DateTime[]
  n1_sync = Float64[]
  n2_sync = Float64[]
  dxPeaks_sync = Float64[]
  peakCenter1_sync = Float64[]
  peakCenter2_sync = Float64[]

  df1[:date] = DateTime[DateTime(tStr) for tStr in df1[:date]]
  df2[:date] = DateTime[DateTime(tStr) for tStr in df2[:date]]


  t1 = convert(Array, df1[:date])
  t2 = convert(Array, df2[:date])
  n1 = convert(Array, df1[:peakArea])
  n2 = convert(Array, df2[:peakArea])
  peakCenter1 = convert(Array, df1[:peakCenter])
  peakCenter2 = convert(Array, df2[:peakCenter])

  for i=1:length(t2)
    for j=1:length(t1)
      if ((t2[i] > t1[j]) & ((t2[i] - t1[j]) < dtMax))
        push!(n1_sync, n1[j])
        push!(n2_sync, n2[i])
        push!(t_sync, t2[i])
        push!(dxPeaks_sync, (peakCenter1[j] - peakCenter2[i]))
        push!(peakCenter1_sync, peakCenter1[j])
        push!(peakCenter2_sync, peakCenter2[i])
        break
      end
    end
  end
  return t_sync, n1_sync, n2_sync, dxPeaks_sync, peakCenter1_sync, peakCenter2_sync

end

function doCO2Correction(tCO, nCO, tCO2, nCO2)

  dtMax = DateTime(2000,1,1,0,10) - DateTime(2000,1,1,0,0)
  nCOCorrected = Float64[]
  tCOCorrected = DateTime[]
  wasCorrected = zeros(Int64, length(nCO))
  for i=1:length(nCO)
    for j=1:length(nCO2)
      if ((tCO2[j] > tCO[i]) & ((tCO2[j] - tCO[i]) < dtMax))
        newValue = nCO[i] - 0.099 * nCO2[j]
        if newValue > 0
          push!(nCOCorrected, newValue)
          push!(tCOCorrected, tCO[i])
          wasCorrected[i] = 1
          break
        end
      end
    end
  end
  return tCOCorrected, nCOCorrected, wasCorrected
end

function clean_h5(path)
  fileNames = readdir(path)
  for fname in fileNames
    fid = h5open(joinpath(path, fname),"r+")
    for grp in fid
      for dataset in grp
        if ismatch(r"M06.2", name(dataset))
          println(name(dataset))
          o_delete(grp, name(dataset))
        end
      end
    end
    close(fid)
  end
end



end
