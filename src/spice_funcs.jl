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
