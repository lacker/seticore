#pragma once

#include <string>
#include <vector>

using namespace std;

/*
  The FilterbankMetadata stores the metadata that is typically stored in conjunction
  with filterbank data.
*/
class FilterbankMetadata {
 public:
  FilterbankMetadata();

  bool has_dc_spike;

  // Represents the boresight source when there are multiple beams
  string source_name;

  // These vectors are for the case where this file represents a number of different
  // beams pointed at different things.
  // It's a bit of a hack and if there is more per-beam metadata like this, we might
  // want a different interface.
  vector<string> source_names;
  vector<double> ras;  // hours
  vector<double> decs; // degrees

  double fch1, foff, tstart, tsamp;

  // src_raj should be in hours (i.e. radians * 12/pi)
  double src_raj;

  // src_dej should be in degrees (i.e. radians * 180/pi)
  double src_dej;

  long num_timesteps, num_channels, coarse_channel_size, num_coarse_channels;
  int telescope_id;
  
  virtual ~FilterbankMetadata() {}

  FilterbankMetadata getSubsetMetadata(int beam, int band, int num_bands);

  // Providing information about beams when metadata is split among beams
  bool isCoherentBeam(int beam) const;
  string getBeamSourceName(int beam) const;

  // Hours
  double getBeamRA(int beam) const;

  // Degrees
  double getBeamDec(int beam) const;
  
 protected:
  void inferMetadata();
};
