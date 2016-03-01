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

function moving_avg(y, window_size)
  myAvg = copy(y)
  for i=window_size+1:length(y)-window_size
    myAvg[i] = median(y[i-window_size:i+window_size])
  end
  return myAvg
end

function pds2h5(path, homeDir="/home/abieler/rosetta/data/dfms/h5/")
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
