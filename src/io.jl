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
  baseline_A = polyFit(row_A, baseline_model)
  baseline_B = polyFit(row_B, baseline_model)

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
function getGainFactor(gainStep,row)
  fileName = ""
  gainFactor = 0.0
  if (row == _rowA)
    fileName = "gainFactor_SPACE_A.txt"
  elseif (row == _rowB)
    fileName = "gainFactor_SPACE_B.txt"
  end

  iFile = open(joinpath(homedir(), ".julia/v0.4/DFMS/InstrumentData", fileName), "r")
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


function get_gain_factor(gainstep, isFMdata)
  if isFMdata
    gain_factor_fm(gainstep)
  else
    gain_factor_pds(gainstep)
  end
end

function gain_factor_pds(gainstep)
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

    iFile = open(joinpath(homedir(), ".julia/v0.4/DFMS/InstrumentData", fileName), "r")
    #iFile = open(joinpath("C:\\Users\\leroy\\.julia\\v0.4\\DFMS\\InstrumentData", fileName), "r")
    i=0
    while !eof(iFile)
      line = readline(iFile)
      if (i == gainstep)
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

function gain_factor_fm(gainstep)
  gainFactor = [0.0, 0.0]
  fileName = "gainFactor_FM.txt"
  iFile = open(joinpath(homedir(), ".julia/v0.4/DFMS/InstrumentData", fileName), "r")
  i=1
  while !eof(iFile)
    line = readline(iFile)
    if i == gainstep
      px, volt, gf_a, gf_b = split(line, '\t')
      gf_a = parse(Float64, gf_a)
      gf_b = parse(Float64, gf_b)
      gainFactor[1] = gf_a
      gainFactor[2] = gf_b
    end
    i+=1
  end
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

    iFile = open(joinpath(homedir(), ".julia/v0.4/DFMS/InstrumentData", fileName), "r")
    #iFile = open(joinpath("C:\\Users\\leroy\\.julia\\v0.4\\DFMS\\InstrumentData", fileName), "r")
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

function get_dfms_data(fileName)
  fileType = fileName[end-2:end]
  if fileType == "DAT"
    return load_fm(fileName)
  elseif fileType == "TAB"
    return load_pds(fileName)
  else
    println("File type unknown:")
    @show(fileType)
  end
end

function load_fm(fileName)
  isFMdata = true
  i = 1
  y = zeros(Float64, 512, 2)
  t = DateTime(2000)
  gainStep = 0
  m0 = 0
  p_pt = 0.0
  iSkip = 130
  isStartTimeFound = false

  iFile = open(fileName, "r")
  while !eof(iFile)
    line = readline(iFile)
    if (contains(line, "Packet time") & (isStartTimeFound == false))
      tStr = matchall(r"\d+/\d+/\d+ \d+:\d+:\d+", line)[1]
      t = DateTime(tStr, "dd/mm/yyyy HH:MM:SS")
      isStartTimeFound = true
    elseif (contains(line, "Gain") & !contains(line, "(pre") & !contains(line, "@"))
      gainStep = round(Int, parse(Float64, matchall(r"(\d+.\d+e\+\d+)", line)[1]))
    elseif contains(line , "Mass     ")
      m0 = round(Int, parse(Float64, matchall(r"(\d+.\d+e\+\d+)", line)[1]))
    elseif contains(line, "p/p_T")
      p_pt = round(Int, parse(Float64, matchall(r"(\d+.\d+e\+\d+)", line)[1]))
    end

    if (i >= iSkip)
      y[i-(iSkip-1), 1] = parse(Float64, matchall(r"(\d+)", line)[2])
      y[i-(iSkip-1), 2] = parse(Float64, matchall(r"(\d+)", line)[3])
    end
    i += 1
  end
  close(iFile)
  return y, gainStep, t, m0, p_pt, isFMdata
end

function load_pds(fileName)
  isFMdata = false
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
  return y, gainStep, t, m0, p_pt, isFMdata
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

function get_pixel_gain(t, gainstep, isFMdata)
  if isFMdata
    pixel_gain_fm(gainstep)
  else
    pixel_gain_pds(t, gainstep)
  end
end

function pixel_gain_fm(gainstep)
  pixelGain = zeros(Float64, 512, 2)
  fileName = joinpath(homedir(), ".julia/v0.4/DFMS/InstrumentData/pixelgain_FM/pixelgain_FM_2015_06.txt")
  iCol = (gainstep - 1) * 6 + 5
  data = readdlm(fileName, '\t', Float64)
  for i=1:size(data,1)
    pixelGain[i,1] = data[i,iCol]
    pixelGain[i,2] = data[i,iCol+1]
  end
  return pixelGain
end

function pixel_gain_pds(t, gainstep)
  fileName = ""
  i = 0
  path = joinpath(homedir(), ".julia/v0.4/DFMS/InstrumentData/pixelgain_SPACE")
  pixelGain = zeros(Float64, 512, 2)
  for row in [_rowA, _rowB]
    # pixel gains for july to 1 october
    if DateTime(2014,2,1) < t < DateTime(2014,10,1)
      fileName = joinpath(path, "pg_space_Aug2014.csv")
      i=2
    # pixel gains from 1 october to 31 december
    elseif DateTime(2014,10,1) <= t < DateTime(2015,1,1)
      if gainstep < 12
        fileName = joinpath(path, "pg_space_Nov2014_GS9_new.csv")
      elseif 12 <= gainstep < 15
        fileName = joinpath(path, "pg_space_Nov2014_GS14_new.csv")
      elseif gainstep >= 15
        fileName = joinpath(path, "pg_space_Nov2014_GS16_new.csv")
      end
      i=1
    # pixel gains from 1 january to april
    elseif DateTime(2015,1,1) <= t < DateTime(2016,4,1)
      if gainstep < 12
        fileName = joinpath(path, "pg_space_Feb2015_GS9_new.csv")
      elseif 12 <= gainstep < 15
        fileName = joinpath(path, "pg_space_Feb2015_GS14_new.csv")
      elseif gainstep >= 15
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


function getPixelGain(t, gainStep)
  fileName = ""
  i = 0
  path = joinpath(homedir(), ".julia/v0.4/DFMS/InstrumentData/pixelgain_SPACE")
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
