#include <capnp/serialize-packed.h>
#include <errno.h>
#include <fmt/core.h>
#include <iostream>
#include "raw_file_group_reader.h"
#include "stamp.capnp.h"
#include "stamp_extractor.h"
#include "upchannelizer.h"
#include "util.h"

using namespace std;

StampExtractor::StampExtractor(RawFileGroup& file_group, int fft_size, int telescope_id,
                               const string& output_filename)
  : file_group(file_group), fft_size(fft_size), telescope_id(telescope_id),
    output_filename(output_filename), opened(false) {
}

StampExtractor::~StampExtractor() {
  if (opened) {
    close(fd);
  }
}

void StampExtractor::openOutputFile() {
  if (opened) {
    return;
  }
  fd = open(output_filename.c_str(), O_WRONLY | O_CREAT, 0664);
  if (fd < 0) {
    int err = errno;
    cerr << "could not open " << output_filename << " for writing. errno = "
         << err << endl;
    exit(1);
  }
  opened = true;
}

void StampExtractor::extract(int coarse_channel, int start_channel, int num_channels) {
  double global_fch1 = file_group.getFch1(fft_size);

  // Output metadata
  double foff = file_group.obsbw / file_group.num_coarse_channels / fft_size;
  double fch1 = global_fch1 + (coarse_channel * fft_size + start_channel) * foff;
  
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
                              1, file_group.npol, file_group.nants);

  ComplexBuffer internal(upchannelizer.requiredInternalBufferSize());

  MultiantennaBuffer fine(upchannelizer.numOutputTimesteps(),
                          upchannelizer.numOutputChannels(),
                          file_group.npol,
                          file_group.nants);

  MultiantennaBuffer output(fine.num_timesteps * num_batches,
                            num_channels,
                            file_group.npol,
                            file_group.nants);

  // Read just a single coarse channel, so one band equals one coarse channel
  RawFileGroupReader reader(file_group, file_group.num_coarse_channels,
                            coarse_channel, coarse_channel,
                            num_batches, blocks_per_batch);

  // Track where in output we're writing to
  int output_time = 0;
  
  for (int batch = 0; batch < num_batches; ++batch) {
    shared_ptr<DeviceRawBuffer> device_raw_buffer = reader.readToDevice();

    upchannelizer.run(*device_raw_buffer, internal, fine);

    fine.copyRange(start_channel, output, output_time);

    output_time += fine.num_timesteps;
  }

  cudaDeviceSynchronize();

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

  openOutputFile();
  writeMessageToFd(fd, message);
}
