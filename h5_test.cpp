#include "catch/catch.hpp"

#include <boost/filesystem.hpp>

#include "h5_reader.h"
#include "h5_writer.h"

TEST_CASE("h5 write then read", "[h5]") {
  string dir = boost::filesystem::temp_directory_path().c_str();
  string filename = dir + "/testing.h5";
  boost::filesystem::remove(filename);
    
  FilterbankMetadata m;
  m.source_name = "bob";
  m.fch1 = 1.0;
  m.foff = 2.0;
  m.tstart = 3.0;
  m.tsamp = 4.0;
  m.src_dej = 5.0;
  m.src_raj = 6.0;
  m.num_timesteps = 10;
  m.num_channels = 100;
  m.coarse_channel_size = 10;
  m.num_coarse_channels = 10;
  m.telescope_id = MEERKAT;

  vector<float> data;
  for (int time = 0; time < m.num_timesteps; ++time) {
    for (int chan = 0; chan < m.num_channels; ++chan) {
      data.push_back(100.0 * time + 1.0 * chan);
    }
  }
  
  H5Writer writer(filename, m);
  writer.setData(&data[0]);
  writer.close();

  H5Reader f(filename);

  REQUIRE(f.source_name == m.source_name);
  REQUIRE(f.fch1 == m.fch1);
  REQUIRE(f.foff == m.foff);
  REQUIRE(f.tstart == m.tstart);
  REQUIRE(f.tsamp == m.tsamp);
  REQUIRE(f.src_dej == m.src_dej);
  REQUIRE(f.src_raj == m.src_raj);
  REQUIRE(f.num_timesteps == m.num_timesteps);
  REQUIRE(f.num_channels == m.num_channels);
  REQUIRE(f.coarse_channel_size == m.coarse_channel_size);
  REQUIRE(f.num_coarse_channels == m.num_coarse_channels);
  REQUIRE(f.telescope_id == m.telescope_id);

  // Spot check data
  FilterbankBuffer buffer(f.num_timesteps, f.coarse_channel_size);
  f.loadCoarseChannel(2, &buffer);
  REQUIRE(buffer.get(1, 3) == 123.0);
  
  boost::filesystem::remove(filename);
}
