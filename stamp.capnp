# To regenerate the files that are based on this schema, run:
# capnp compile -oc++ stamp.capnp

# unique id generated by capnp
@0xb811e7262df2bb01;

using Hit = import "hit.capnp";

# A "postage stamp" of data extracted from a larger set of raw data.
# This data has been upchannelized with an FFT but is still multi-antenna
# complex voltages.
struct Stamp {
  # The seticore version that generated this data.     
  seticoreVersion @13 :Text;

  # source, ra, and dec refer to the boresight target.
  sourceName @0 :Text;
  ra @1 :Float64;   # hours
  dec @2 :Float64;  # degrees

  # Other standard metadata found in FBH5 files.
  # This metadata applies specifically to the postage stamp itself, not the larger
  # file we extracted it from.
  fch1 @3 :Float64;  # MHz
  foff @4 :Float64;  # MHz
  tstart @5 :Float64;
  tsamp @6 :Float64;
  telescopeId @7 :Int32;

  # Dimensions of the data
  numTimesteps @8 :Int32;
  numChannels @9 :Int32;
  numPolarizations @10 :Int32;
  numAntennas @11 :Int32;

  # An array of complex voltages. Interpret this as row-major:
  #   data[timestep][channel][polarization][antenna][real vs imag]
  data @12 :List(Float32);

  # Metadata describing how exactly we extracted this stamp.
  # The coarse channel in the original file that we extracted data from.
  # Matches coarseChannel in the hit.
  coarseChannel @14 :Int32;

  # The size of FFT we used to create fine channels.
  fftSize @15 :Int32;

  # The first post-FFT channel that we extracted.
  # So column zero in `data` corresponds to this column in the original
  # post-FFT data.
  # This will not exactly match startChannel in a hit, because we combine
  # adjacent hits and may use different window sizes. But if you consider the
  # intervals [startChannel, startChannel + numChannels), the interval for
  # the stamp should overlap with the interval for any relevant hits.
  startChannel @16 :Int32;

  # Metadata for the best hit we found for this stamp.
  # Not populated for stamps extracted with the CLI tool.
  signal @17 :Hit.Signal;

  # Metadata copied from the input raw file.
  # We need these so that we can match up this stamp to the beamforming recipe file.
  # "schan" is where the raw file starts in the beamforming recipe.
  schan @18 :Int32;
  obsid @19 :Text;
}
