#include "raw_file_group.h"

#include <assert.h>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/filesystem.hpp>
#include <fmt/core.h>
#include <iostream>
#include "util.h"
#include <vector> 

using namespace std;

// Extracts the prefix from a file name of the format:
//  <prefix>.<sequence-identifier>.raw
// This will include the directory.
// Returns the empty string if the filename doesn't match the format.
string getRawFilePrefix(const string& filename) {
  if (!boost::algorithm::ends_with(filename, ".raw")) {
    return "";
  }
  string no_suffix = filename.substr(0, filename.size() - 4);
  auto index = no_suffix.find_last_of(".");
  if (index >= no_suffix.size()) {
    // We can't find a prefix
    return "";
  }
  return no_suffix.substr(0, index);
}

// Given /foo/bar/baz, returns baz
string getBasename(const string& filename) {
  auto index = filename.find_last_of("/");
  if (index > filename.size()) {
    return filename;
  }
  return filename.substr(index + 1);
}

string getDirectory(const string& filename) {
  auto index = filename.find_last_of("/");
  assert(0 <= index && index < filename.size());
  return filename.substr(0, index);
}

vector<string> getRawFilesMatchingPrefix(const string& prefix) {
  string directory = getDirectory(prefix);
  boost::filesystem::path dir{directory};
  vector<string> filenames;
  boost::filesystem::directory_iterator it{dir};
  while (it != boost::filesystem::directory_iterator{}) {
    string filename(it->path().c_str());
    if (boost::starts_with(filename, prefix) &&
        boost::ends_with(filename, ".raw")) {
      filenames.push_back(filename);
    }
    ++it;
  }
  sort(filenames.begin(), filenames.end());
  return filenames;
}

RawFileGroup::RawFileGroup(const vector<string>& filenames)
  : current_file(-1), band(-1), num_bands(-1), read_size(-1), filenames(filenames) {
  assert(!filenames.empty());
  prefix = getBasename(getRawFilePrefix(filenames[0]));

  // Get metadata from the first file
  const RawFile& file(openFile(filenames[0]));
  const raw::Header& header(file.headers().front());

  nants = header.nants;
  num_coarse_channels = header.num_channels;
  npol = header.npol;
  obsbw = header.obsbw;
  obsfreq = header.obsfreq;
  ra = header.ra;
  dec = header.dec;
  obsid = header.getString("OBSID");
  source_name = header.src_name;
  telescope = header.telescop;
  start_pktidx = header.pktidx;
  next_pktidx = start_pktidx;
  tbin = header.tbin;
  timesteps_per_block = header.num_timesteps;

  schan = header.getInt("SCHAN", -1);
  assert(schan >= 0);
  start_time = header.getStartTime();

  if (header.mjd != 0) {
    // Sanity check
    double calculated_mjd = unixTimeToMJD(start_time);
    double read_mjd = header.mjd;
    if (fabs(calculated_mjd - read_mjd) > 0.001) {
      fatal(fmt::format("MJD mismatch. calculated {} but header has {} in {}",
                        calculated_mjd, read_mjd, filenames[0]));
    }
  }
  
  piperblk = header.getUnsignedInt("PIPERBLK", 0);
  assert(piperblk > 0);

  // Find the last block in the last file
  const RawFile& last_file(openFile(filenames[filenames.size() - 1]));
  int pktidx_diff = last_file.headers().back().pktidx - start_pktidx;
  assert(0 == pktidx_diff % piperblk);
  num_blocks = (pktidx_diff / piperblk) + 1;
}

RawFileGroup::~RawFileGroup() {}

/*
  Group up raw files of the form:
    <prefix>.<sequence-identifier>.raw

  Each prefix defines one group of raw files.
  The groups are sorted by prefix.
  Within the group, they are sorted by sequence identifier.
  A valid group must start with the sequence identifier "0000".
  Each of these sorts are *string* sorting, so if you provide raw files
  with numerical ids, be sure that they are the same length.
 */
vector<vector<string> > scanForRawFileGroups(const string& directory) {
  if (!boost::filesystem::is_directory(directory)) {
    fatal(directory, "is not a directory");
  }
  boost::filesystem::path dir{directory};
  vector<string> filenames;
  boost::filesystem::directory_iterator it{dir};
  while (it != boost::filesystem::directory_iterator{}) {
    filenames.push_back(it->path().c_str());
    ++it;
  }
  sort(filenames.begin(), filenames.end());

  // Gather filenames that share the group_prefix in group
  vector<string> group;
  string group_prefix("");
  vector<vector<string> > groups;

  for (string filename : filenames) {
    string prefix = getRawFilePrefix(filename);
    if (prefix.empty()) {
      continue;
    }

    if (group.empty()) {
      // This file is the first one overall
      group.push_back(filename);
      group_prefix = prefix;
      continue;
    }

    if (prefix == group_prefix) {
      // This file is a continuation of the current group
      group.push_back(filename);
      continue;
    }

    // This file represents a new group
    groups.push_back(group);
    group.clear();
    group.push_back(filename);
    group_prefix = prefix;
  }

  if (!group.empty()) {
    // We have to add the last group to our groups
    groups.push_back(group);
  }

  // Accept only the ones starting with a .0000.raw file
  vector<vector<string> > answer;
  for (auto group : groups) {
    if (boost::algorithm::ends_with(group[0], ".0000.raw")) {
      answer.push_back(group);
    }
  }
  
  return answer;
}

void RawFileGroup::resetBand(int new_band, int new_num_bands) {
  assert(0 <= new_band && new_band < new_num_bands);
  assert(num_coarse_channels % new_num_bands == 0);
  band = new_band;
  num_bands = new_num_bands;

  // Calculate the size of each read
  int data_size = oneBlockDataSize();
  assert(data_size % num_bands == 0);
  read_size = data_size / num_bands;
  
  // Prepare for iteration
  current_file = -1;
  next_pktidx = start_pktidx;
}

const RawFile& RawFileGroup::getFile() {
  return openFile(filenames[current_file]);
}

const raw::Reader& RawFileGroup::getReader() {
  return getFile().reader();
}

const raw::Header& RawFileGroup::getHeader() {
  return getFile().headers()[header_index];
}

const RawFile& RawFileGroup::openFile(const string& filename) {
  if (files.find(filename) == files.end()) {
    files[filename] = make_unique<RawFile>(filename);
  }
  return *files[filename];
}

void RawFileGroup::openNextFile() {
  ++current_file;
  header_index = 0;
  assert(current_file < (int) filenames.size());

  const raw::Header& header(getHeader());
  
  // Sanity check some values in this header
  assert(schan == header.getInt("SCHAN", -1));
  assert(nants == header.nants);
  assert(num_coarse_channels == header.num_channels);
  assert(npol == (int) header.npol);  
}

void RawFileGroup::readTasks(char* buffer, vector<function<bool()> >* tasks) {
  if (current_file == -1) {
    openNextFile();
  }

  const raw::Header& header(getHeader());
  if (header.pktidx < next_pktidx) {
    // Advance the header index
    ++header_index;    

    // Check if the header index is still valid
    const RawFile& file = getFile();
    if (header_index >= (int) file.headers().size()) {
      // We need to move on to the next file
      openNextFile();
      const raw::Header& h(getHeader());
      if (h.pktidx < next_pktidx) {
        fatal(fmt::format("first pktidx in {} is only {} when we expected at least {}",
                          getFile().filename, h.pktidx, next_pktidx));
      }
    }
  }

  const raw::Header& h(getHeader());
  if (h.pktidx == next_pktidx) {
    // We're pointing at the data we want to return
    getReader().readBandTasks(h, band, num_bands, buffer, tasks);
    next_pktidx += piperblk;
    return;
  }

  // We missed some blocks, so we'll have to return some zeros
  assert(h.pktidx > next_pktidx);
  cout << "missing block with pktidx = " << next_pktidx << endl;
  memset(buffer, 0, read_size);
  next_pktidx += piperblk;  
}

// Threadsafe.
// We can't calculate the time from the actual block, because we
// might have missed that block.
double RawFileGroup::getStartTime(int block) const {
  assert(block < num_blocks);
  double time_per_block = tbin * timesteps_per_block;
  return start_time + block * time_per_block;
}

float RawFileGroup::totalTime() const {
  return tbin * timesteps_per_block * num_blocks;
}

float RawFileGroup::totalDataGB() const {
  float giga = 1024.0 * 1024.0 * 1024.0;
  float data_per_block = oneBlockDataSize() / giga;
  return num_blocks * data_per_block;
}

float RawFileGroup::coarseChannelBandwidth() const {
  return obsbw / num_coarse_channels;
}

/*
  Gets the frequency for the first channel of this file after it's FFT'd.

  This formula is weird because FFTs aren't just dividing up buckets into
  smaller buckets. See:
    https://github.com/UCBerkeleySETI/bl_docs/wiki/calculating-fch1
*/
double RawFileGroup::getFch1(int fft_size) const {
  double fcchan0 = obsfreq - obsbw * (num_coarse_channels - 1) / (2 * num_coarse_channels);
  double fch1 = fcchan0 - floor(fft_size / 2) * obsbw / (num_coarse_channels * fft_size);
  return fch1;
}

int RawFileGroup::oneBlockDataSize() const {
  return nants * num_coarse_channels * timesteps_per_block * npol * 2;
}

int RawFileGroup::getTelescopeID() const {
  return telescopeID(telescope);
}
