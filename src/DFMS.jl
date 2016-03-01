module DFMS
using DataFrames
using LsqFit
using PyPlot
using HDF5
using Spice

include("typedefs.jl")
include("io.jl")
include("spice_funcs.jl")
include("fitting.jl")
include("plot_funcs.jl")
include("data_transform.jl")


export crop,
       clean,
       clean_h5,
       findPeaks,
       get_dfms_data,
       getGainFactor,
       get_gain_factor,
       getPixelGain,
       get_pixel_gain,
       get_eph,
       get_ephemeris_spice,
       get_outliers,
       load_h5_data,
       load_h5_data!,
       baseline_model,
       moving_avg,
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

#=
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
=#

end
