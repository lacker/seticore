#include "raw_file_group.h"

#include <assert.h>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/filesystem.hpp>
#include <fmt/core.h>
#include <iostream>
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

RawFileGroup::RawFileGroup(const vector<string>& filenames, int num_bands)
  : current_file(-1), filenames(filenames), band(0),
    num_bands(num_bands) {
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
  obsid = header.getString("OBSID");
  source_name = header.src_name;
  start_pktidx = header.pktidx;
  next_pktidx = start_pktidx;
  tbin = header.tbin;
  timesteps_per_block = header.num_timesteps;

  schan = header.getInt("SCHAN", -1);
  assert(schan > 0);
  synctime = header.getUnsignedInt("SYNCTIME", -1);
  assert(synctime > 0);
  piperblk = header.getUnsignedInt("PIPERBLK", -1);
  assert(piperblk > 0);

  // Calculate the size of each read
  assert(0 == num_coarse_channels % num_bands);
  int channels_per_band = num_coarse_channels / num_bands;
  read_size = nants * channels_per_band * timesteps_per_block * npol * 2;
  
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
  Each of these sorts are *string* sorting, so if you provide raw files
  with numerical ids, be sure that they are the same length.
 */
vector<vector<string> > scanForRawFileGroups(const string& directory) {
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
  vector<vector<string> > answer;

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
    answer.push_back(group);
    group.clear();
    group.push_back(filename);
    group_prefix = prefix;
  }

  if (!group.empty()) {
    // We have to add the last group to our answer
    answer.push_back(group);
  }
  
  return answer;
}

void RawFileGroup::resetBand(int new_band) {
  assert(new_band < num_bands);
  band = new_band;
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
    files[filename] = make_unique<RawFile>(filename, num_bands);
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
        cerr << "first pktidx in " << getFile().filename << " is only "
             << h.pktidx << " when we expected at least " << next_pktidx
             << endl;
        exit(1);
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
  return synctime + block * time_per_block;
}

float RawFileGroup::totalTime() const {
  return tbin * timesteps_per_block * num_blocks;
}
