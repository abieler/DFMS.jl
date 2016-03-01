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
