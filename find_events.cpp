#include <assert.h>
#include "dedoppler.h"
#include "event_file_writer.h"
#include "filterbank_buffer.h"
#include "filterbank_file_reader.h"
#include "find_events.h"
#include <fmt/core.h>
#include <map>
#include "util.h"

#include "cuda_util.h"

using namespace std;

struct FindHitResult {
  int count;

  // Returns a hit iff there is precisely one result.
  DedopplerHit* hit;
};

// low_index-high_index is inclusive.
// This method only cares about finding unique hits.
// If we find more than one result, just returns a count of 2 rather
// than counting them all.
FindHitResult findHit(const map<int, DedopplerHit*>& hitmap,
                      int low_index, int high_index) {
  auto iter = hitmap.lower_bound(low_index);

  if (iter == hitmap.end()) {
    // There's nothing in this range
    return FindHitResult{ 0, NULL };
  }

  FindHitResult answer{ 1, iter->second };

  ++iter;

  if (iter == hitmap.end() || iter->first > high_index) {
    return answer;
  }

  // We found at least two hits in the range.
  return FindHitResult{ 2, NULL };
}

/*
  Runs a dedoppler algorithm across files in a cadence, assuming they follow an
  "ABACAD" pattern.
  The results that appear in "on" but not "off" data are written to an .events file.

  max_drift is the maximum drift we are looking for, in Hz/sec
  snr_on is the minimum snr we require in the "on" data.
  snr_off is the minimum snr that makes us throw out an event when it
  occurs in the "off" data.
 */
void findEvents(const vector<string>& input_filenames, const string& output_filename,
                double max_drift, double snr_on, double snr_off) {
  vector<shared_ptr<FilterbankFileReader> > files;
  vector<shared_ptr<Dedopplerer> > dedopplerers;
  vector<shared_ptr<FilterbankBuffer> > buffers;

  for (auto& filename : input_filenames) {
    files.push_back(move(loadFilterbankFile(filename)));
    const auto& file = files.back();

    shared_ptr<Dedopplerer> dedopplerer(new Dedopplerer(file->num_timesteps,
                                                        file->coarse_channel_size,
                                                        file->foff, file->tsamp,
                                                        file->has_dc_spike));
    dedopplerers.push_back(dedopplerer);

    shared_ptr<FilterbankBuffer>
      buffer(new FilterbankBuffer(roundUpToPowerOfTwo(file->num_timesteps),
                                  file->coarse_channel_size));
    buffers.push_back(buffer);
  }

  EventFileWriter writer(output_filename, files);
  
  // Check the metadata lines up
  int num_timesteps = files[0]->num_timesteps;
  int coarse_channel_size = files[0]->coarse_channel_size;
  int num_coarse_channels = files[0]->num_coarse_channels;
  double foff = files[0]->foff;
  double tsamp = files[0]->tsamp;  
  bool has_dc_spike = files[0]->has_dc_spike;
  string source_name = files[0]->source_name;

  for (int i = 0; i < (int) files.size(); ++i) {
    assert(files[i]->num_timesteps == num_timesteps);
    assert(files[i]->coarse_channel_size == coarse_channel_size);
    assert(files[i]->num_coarse_channels == num_coarse_channels);
    assertFloatEq(files[i]->foff, foff);
    assertFloatEq(files[i]->tsamp, tsamp);
    assert(files[i]->has_dc_spike == has_dc_spike);

    // The targets should work like ABACAD
    if (i % 2 == 0) {
      if (files[i]->source_name != source_name) {
        fatal(fmt::format("file {} has source {} when we expected source {}",
                          files[i]->filename, files[i]->source_name, source_name));
      }
    } else {
      if (files[i]->source_name == source_name && source_name != "VOYAGER-1") {
        // The voyager cadence is bad like this
        fatal(fmt::format("file {} has source {} but it is supposed to be an 'off'",
                          files[i]->filename, files[i]->source_name));
      }
    }
  }

  // Handle one coarse channel at a time
  for (int coarse_channel = 0; coarse_channel < num_coarse_channels; ++coarse_channel) {
    // Each hit list in hit_lists corresponds to a single file. 
    vector<vector<DedopplerHit>> hit_lists(files.size());

    // Always scan the first file
    files[0]->loadCoarseChannel(coarse_channel, buffers[0].get());
    dedopplerers[0]->search(*buffers[0], *files[0], NO_BEAM, coarse_channel, max_drift,
                            0.0, snr_on, &hit_lists[0]);

    if (hit_lists[0].empty()) {
      // No hits in this coarse channel
      continue;
    }
    
    // Scan the rest of the files
    for (int i = 1; i < (int) dedopplerers.size(); ++i) {
      files[i]->loadCoarseChannel(coarse_channel, buffers[i].get());
      dedopplerers[i]->search(*buffers[i], *files[i], NO_BEAM, coarse_channel, max_drift,
                              0.0, (i % 2 == 0) ? snr_on : snr_off, &hit_lists[i]);
    }

    // For each input file, make a map keying each hit by their
    // starting index.
    // This should be unique because the dedopplerer will already only
    // report one hit per index.
    // This will let us match up the hits for events without doing a
    // linear scan for each candidate.
    vector<map<int, DedopplerHit*> > hitmaps(hit_lists.size());
    for (int i = 0; i < (int) hit_lists.size(); ++i) {
      for (int j = 0; j < (int) hit_lists[i].size(); ++j) {
        DedopplerHit* hit = &hit_lists[i][j];
        hitmaps[i][hit->index] = hit;
      }
    }

    // For each hit in the first file, we build a potential event
    // candidate
    double initial_tstart = files[0]->tstart;
    int num_events = 0;
    for (const auto& pair : hitmaps[0]) {
      auto initial_hit = pair.second;
      
      // First search the "ons" to make sure each "on" has a hit in
      // the right place and find the total frequency range we're looking over
      bool looks_ok = true;
      int low_index = INT_MAX;
      int high_index = -1;
      for (int i = 0; i < (int) hitmaps.size(); i += 2) {
        // Figure out where we expect to see a hit
        double delta_seconds = files[i]->tstart - initial_tstart;
        int timesteps = round(delta_seconds / tsamp);
        int expected_index = initial_hit->expectedIndex(timesteps);
        int wiggle = 10;
        FindHitResult result = findHit(hitmaps[i], expected_index - wiggle,
                                       expected_index + wiggle);
        if (result.count != 1) {
          looks_ok = false;
          break;
        }

        low_index = min(low_index, result.hit->lowIndex());
        high_index = max(high_index, result.hit->highIndex());
      }
      if (!looks_ok) {
        continue;
      }

      // If these fail, there's a bug in the code
      assert(low_index < coarse_channel_size);
      assert(high_index >= 0);
      
      // Now search every file to make sure it has the right number of
      // hits in the range. "on"s should have 1, "off"s should have 0.
      bool candidate_good = true;
      vector<DedopplerHit*> candidate;
      for (int i = 0; i < (int) hitmaps.size(); ++i) {
        int wiggle = 50;
        int ideal_count = i % 2 == 0 ? 1 : 0;
        FindHitResult result = findHit(hitmaps[i], low_index - wiggle,
                                       high_index - wiggle);
        if (result.count != ideal_count) {
          candidate_good = false;
          break;
        }
        candidate.push_back(result.hit);
      }

      if (!candidate_good) {
        continue;
      }

      // We actually have a good candidate. Write it out
      cout << "found event starting at " << initial_hit->toString() << endl;
      ++num_events;
      writer.write(candidate, buffers);
    }
    cout << pluralize(num_events, "event") << " found" << endl;
  }
}
