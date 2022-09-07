#include <boost/program_options.hpp>
#include <capnp/message.h>
#include <capnp/serialize-packed.h>
#include <errno.h>
#include <fmt/core.h>
#include <iostream>
#include "raw_file_group.h"
#include "raw_file_group_reader.h"
#include "stamp.capnp.h"
#include "upchannelizer.h"
#include "util.h"

using namespace std;

namespace po = boost::program_options;

/*
  Extract a single "postage stamp" from a group of raw files.
 */
int main(int argc, char* argv[]) {
  po::options_description desc("extract options");
  desc.add_options()

    ("raw_prefix", po::value<string>(),
     "prefix for the raw files to extract a stamp from")

    ("output", po::value<string>(),
     "location to write the stamp file")

    ("band", po::value<int>(),
     "parameter for raw-file-group reading")

    ("num_bands", po::value<int>(),
     "parameter for raw-file-group reading")

    ("coarse_channel", po::value<int>(),
     "coarse channel to extract data from (within the band)")

    ("fft_size", po::value<int>(),
     "size of the fft to use for upchannelization")

    ("start_channel", po::value<int>(),
     "channel to start extraction after fft")

    ("num_channels", po::value<int>(),
     "how many post-fft channels to extract")

    ("telescope_id", po::value<int>(),
     "telescope id")
    
    ;

  
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
  po::notify(vm);

  string raw_prefix = vm["raw_prefix"].as<string>();
  string output_filename = vm["output"].as<string>();
  int band = vm["band"].as<int>();
  int num_bands = vm["num_bands"].as<int>();
  int coarse_channel = vm["coarse_channel"].as<int>();
  int fft_size = vm["fft_size"].as<int>();
  int start_channel = vm["start_channel"].as<int>();
  int num_channels = vm["num_channels"].as<int>();
  int telescope_id = vm["telescope_id"].as<int>();
  
  cout << fmt::format("extracting stamp from {}.*.raw\n", raw_prefix);
  vector<string> raw_files = getFilesMatchingPrefix(raw_prefix);
  cout << fmt::format("reading {} raw files:\n", raw_files.size());
  for (auto f : raw_files) {
    cout << "  " << f << endl;
  }

  RawFileGroup file_group(raw_files, num_bands);
  int fine_channel = coarse_channel * fft_size + start_channel;
  double global_fch1 = file_group.obsfreq - 0.5 * file_group.obsbw;

  // Output metadata
  double foff = file_group.obsbw / file_group.num_coarse_channels / fft_size;
  double fch1 = global_fch1 + fine_channel * foff;
  
  // One of the fft size and block duration should divide the other
  int blocks_per_batch;
  if (file_group.timesteps_per_block >= fft_size) {
    assert(file_group.timesteps_per_block % fft_size == 0);
    blocks_per_batch = 1;
  } else {
    assert(fft_size % file_group.timesteps_per_block == 0);
    blocks_per_batch = fft_size / file_group.timesteps_per_block;
  }
  int num_batches = file_group.num_blocks / blocks_per_batch;
  
  Upchannelizer upchannelizer(0, fft_size,
                              file_group.timesteps_per_block * blocks_per_batch,
                              file_group.num_coarse_channels / num_bands,
                              file_group.npol,
                              file_group.nants);

  ComplexBuffer internal(upchannelizer.requiredInternalBufferSize());

  MultiantennaBuffer fine(upchannelizer.numOutputTimesteps(),
                          upchannelizer.numOutputChannels(),
                          file_group.npol,
                          file_group.nants);

  MultiantennaBuffer output(fine.num_timesteps * num_batches,
                            num_channels,
                            file_group.npol,
                            file_group.nants);

  RawFileGroupReader reader(file_group, band, 1, num_batches, blocks_per_batch);

  // Track where in output we're writing to
  int output_time = 0;
  
  for (int batch = 0; batch < num_batches; ++batch) {
    shared_ptr<DeviceRawBuffer> device_raw_buffer = reader.readToDevice();

    upchannelizer.run(*device_raw_buffer, internal, fine);

    fine.copyRange(fine_channel, output, output_time);

    output_time += fine.num_timesteps;
  }

  cudaDeviceSynchronize();

  // Output the file
  int fd = open(output_filename.c_str(), O_WRONLY | O_CREAT, 0664);
  if (fd < 0) {
    int err = errno;
    cerr << "could not open " << output_filename << " for writing. errno = "
         << err << endl;
    exit(1);
  }

  ::capnp::MallocMessageBuilder message;
  Stamp::Builder stamp = message.initRoot<Stamp>();
  stamp.setSourceName(file_group.source_name);
  stamp.setRa(file_group.ra);
  stamp.setDec(file_group.dec);
  stamp.setFch1(fch1);
  stamp.setFoff(foff);
  stamp.setTstart(file_group.getStartTime(0));
  stamp.setTsamp(file_group.tbin * fft_size);
  stamp.setTelescopeId(telescope_id);
  stamp.setNumTimesteps(output.num_timesteps);
  stamp.setNumChannels(output.num_channels);
  stamp.setNumPolarities(output.num_polarities);
  stamp.setNumAntennas(output.num_antennas);
  stamp.initData(2 * output.size);
  auto data = stamp.getData();

  for (int i = 0; i < (int) output.size; ++i) {
    thrust::complex<float> value = output.get(i);
    data.set(2 * i, value.real());
    data.set(2 * i + 1, value.imag());
  }

  writeMessageToFd(fd, message);
  
  // Clean up
  close(fd);
}
