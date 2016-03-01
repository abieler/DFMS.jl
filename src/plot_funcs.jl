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
