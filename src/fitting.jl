function doubleGauss(x, p)
  nPeaks = round(Int, length(p) / 5)
  B = 0.1

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

function singleGauss(x, p)
  # x = bins
  # p = fitting parameters
  # B = baseline height
  #w0 = 4.0
  #dw = 1.25

  nPeaks = round(Int, (length(p)-1) / 2)
  B = 0.1
  w = p[end]
  th = tanh((w-w0)/dw)
  w_th = w0 + th * dw

  yFit = zeros(Float64, length(x))
  for i=1:nPeaks
    A = p[i+nPeaks]

    #########################################
    # limit the fitting parameters
    x0 = p[i]
    th_x = tanh((x0-pI_0[i])/dx)
    x0_th = pI_0[i] + th_x * dx

    #w = p[i+2*nPeaks]
    #th = tanh((w-w0)/dw)
    #w_th = w0 + th * dw
    ##########################################

    yFit += A * exp(-(x-x0_th).^2 / w_th^2)
  end
  yFit += B
  return yFit
end

function baseline_model(x, p)
  f = zeros(Float64, length(x))
  for i=1:length(x)
    @inbounds f[i] = p[1] + p[2]*x[i] + p[3]*x[i]^2 + p[4]*x[i]^3
  end
  return f
end

function peakFit(y, pI, pA, fitMethod)
  # global variables needed to artificially limit the
  # search space for the fitting routine by a tanh
  # transformation.
  # w0 is the base value for the width (sigma) of the
  # peaks, the fit routine can vary them by +/- dw.
  # pI_0 and dx consider the position of the peak
  # centroid.
  global w0 = 4.0
  global dw = 1.50
  global pI_0 = Float64[value for value in pI]
  global dx = 1.25

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
    push!(fitParams, 1.0)
    #append!(fitParams, ones(Float64, length(pI)))
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

    fit = curve_fit(singleGauss, x[xmin:xmax], y[xmin:xmax], weights[xmin:xmax],
                    fitParams)

    for i=1:nPeaks
      peakIndex[i] = fit.param[i]
      peakArea[i] = fit.param[i+nPeaks] * abs(fit.param[end]) * sqrt(pi)
      peakWidth[i] = fit.param[end]
    end

    # transform fitted parameter again
    th = tanh((peakWidth-w0)/dw)
    peakWidth = w0 + th * dw

    th_x = tanh((peakIndex-pI_0)/dx)
    peakIndex = pI_0 + th_x * dx

    return singleGauss(x, fit.param), peakArea, peakIndex, peakWidth
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
    fit = curve_fit(singleGauss, x[20:500], y[20:500], fitParams)
    #fit = curve_fit(singleGauss, x[LHS:RHS], y[LHS:RHS], fitParams)
    peakArea = fit.param[iPeak+nPeaks] * abs(fit.param[iPeak+2*nPeaks]) * sqrt(pi)
    peakIndex = fit.param[iPeak]
    peakWidth = fit.param[iPeak+2*nPeaks]
    return singleGauss(x, fit.param), peakArea, peakIndex, peakWidth
  end
end

function polyFit(y, baseline_model)
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

  fit = curve_fit(baseline_model, xBG, yBG, [median(y), 0., 0., 0.])

  return baseline_model(collect(fullRange), fit.param)
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
      push!(peakIndexes, i)
      push!(peakAmpl, y[i])
    end
  end

  return peakIndexes, peakAmpl
end
