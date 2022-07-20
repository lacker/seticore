#include <boost/algorithm/string.hpp>
#include <fmt/core.h>
#include <string>
#include <vector>

#include "dat_file_writer.h"
#include "h5_reader.h"

using namespace std;

/*
  Opens a dat file for writing, to log hits that we find.
  Some metadata for headers is copied out from the h5 file.
 */
DatFileWriter::DatFileWriter(const string& filename,
                             const FilterbankFileReader& metadata,
                             double max_drift) : metadata(metadata) {
  hit_count = 0;
  file.open(filename);

  vector<string> parts;
  boost::split(parts, filename, boost::is_any_of("/"));
  file << "# -------------------------- o --------------------------\n";
  file << "# File ID: " << parts[parts.size() - 1] << " \n";
  file << "# -------------------------- o --------------------------\n";
  file << "# Source:" << metadata.source_name << "\n";
  file << fmt::format("# MJD: {:18.12f}\tRA: {:.6f}\tDEC: {:.6f}\n",
                      metadata.tstart, metadata.src_raj, metadata.src_dej);
  file << fmt::format("# DELTAT: {:10.6f}\tDELTAF(Hz): {:10.6f}\t"
                      "max_drift_rate: {:10.6f}\tobs_length: {:10.6f}\n",
                      metadata.tsamp, metadata.foff * 1000000.0, max_drift,
                      metadata.num_timesteps * metadata.tsamp);
  file << "# --------------------------\n"
    "# Top_Hit_# \t"
    "Drift_Rate \t"
    "SNR \t"
    "Uncorrected_Frequency \t"
    "Corrected_Frequency \t"
    "Index \t"
    "freq_start \t"
    "freq_end \t"
    "SEFD \t"
    "SEFD_freq \t"
    "Coarse_Channel_Number \t"
    "Full_number_of_hits \t"
    "\n"
    "# --------------------------\n";
  file << flush;
}


DatFileWriter::~DatFileWriter() {
  file.close();
}


// Ignores the beam
void DatFileWriter::recordHit(DedopplerHit hit, int beam, int coarse_channel,
                              const float* input) {
  ++hit_count;

  int global_index = coarse_channel * metadata.coarse_channel_size + hit.index;
  double frequency = metadata.fch1 + global_index * metadata.foff;

  // Currently we just output one frequency for all frequency-type columns.
  // We also just output 1 instead of counting up pre-deduping hits because I suspect
  // that nobody wants it, and it's bothersome to calculate with our current algorithm.
  // TODO: see what the astronomers actually want here
  file << fmt::format("{}\t{:10.6f}\t{:10.6f}\t{:14.6f}\t{:14.6f}\t{}\t{:14.6f}\t{:14.6f}\t",
                      hit_count, hit.drift_rate, hit.snr, frequency,
                      frequency, hit.index, frequency, frequency);
  file << fmt::format("0.0\t0.000000\t{}\t1\t\n", coarse_channel);
  file << flush;
}
