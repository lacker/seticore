#include <assert.h>
#include <cstdlib>
#include "hdf5.h"
#include <iostream>
#include <string.h>

#include "h5_file.h"

using namespace std;


/*
  Opens an h5 file for reading, assuming it contains standard radio
  telescope input from one of the known telescopes. If the data is an
  unexpected size or shape we should be conservative and exit.
 */
H5File::H5File(const string& filename) : FilterbankFile(filename) {
  file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file == H5I_INVALID_HID) {
    cerr << "could not open file: " << filename << endl;
    exit(1);
  }
  dataset = H5Dopen2(file, "data", H5P_DEFAULT);
  if (dataset == H5I_INVALID_HID) {
    cerr << "could not open dataset\n";
    exit(1);
  }
  if (!H5Tequal(H5Dget_type(dataset), H5T_NATIVE_FLOAT)) {
    cerr << "dataset is not float\n";
    exit(1);
  }
    
  fch1 = getDoubleAttr("fch1");
  foff = getDoubleAttr("foff");
  tstart = getDoubleAttr("tstart");
  tsamp = getDoubleAttr("tsamp");
  src_dej = getDoubleAttr("src_dej");
  src_raj = getDoubleAttr("src_raj");
  source_name = getStringAttr("source_name");
  
  dataspace = H5Dget_space(dataset);
  if (dataspace == H5I_INVALID_HID) {
    cerr << "could not open dataspace\n";
    exit(1);
  }
  if (H5Sget_simple_extent_ndims(dataspace) != 3) {
    cerr << "data is not three-dimensional\n";
    exit(1);
  }
  hsize_t dims[3];
  H5Sget_simple_extent_dims(dataspace, dims, NULL);
  if (dims[1] != 1) {
    cerr << "unexpected second dimension: " << dims[1] << endl;
    exit(1);
  }
  num_timesteps = dims[0];
  num_freqs = dims[2];

  telescope_id = getLongAttr("telescope_id");

  if (attrExists("nfpc")) {
    coarse_channel_size = getLongAttr("nfpc");
  }

  inferMetadata();
}

double H5File::getDoubleAttr(const string& name) const {
  double output;
  auto attr = H5Aopen(dataset, name.c_str(), H5P_DEFAULT);
  if (attr == H5I_INVALID_HID) {
    cerr << "could not access attr " << name << endl;
    exit(1);
  }
  if (H5Aread(attr, H5T_NATIVE_DOUBLE, &output) < 0) {
    cerr << "attr " << name << " could not be read as double\n";
    exit(1);
  }
  H5Aclose(attr);
  return output;
}

long H5File::getLongAttr(const string& name) const {
  long output;
  auto attr = H5Aopen(dataset, name.c_str(), H5P_DEFAULT);
  if (attr == H5I_INVALID_HID) {
    cerr << "could not access attr " << name << endl;
    exit(1);
  }
  if (H5Aread(attr, H5T_NATIVE_LONG, &output) < 0) {
    cerr << "attr " << name << " could not be read as long\n";
    exit(1);
  }
  H5Aclose(attr);
  return output;  
}

bool H5File::attrExists(const string& name) const {
  auto answer = H5Aexists(dataset, name.c_str());
  if (answer < 0) {
    cerr << "existence check for attr " << name << " failed\n";
    exit(1);
  }
  return answer > 0;
}

/*
  This assumes the string is stored as variable-length UTF8.
  I'm not sure what if any type conversion the library will do, and
  historically we do not store attributes with consistent string
  subtypes, so when we run into attributes with different formats we
  might have to improve this method.
 */
string H5File::getStringAttr(const string& name) const {
  auto attr = H5Aopen(dataset, name.c_str(), H5P_DEFAULT);
  if (attr == H5I_INVALID_HID) {
    cerr << "could not access attr " << name << endl;
    exit(1);
  }

  // Check the attribute's character type
  auto attr_type = H5Aget_type(attr);
  auto cset = H5Tget_cset(attr_type);
  if (cset < 0) {
    cerr << "H5Tget_cset failed\n";
    exit(1);
  }
  
  // Create mem_type for variable-length string of our character type
  auto mem_type = H5Tcopy(H5T_C_S1);
  if (H5Tset_size(mem_type, H5T_VARIABLE) < 0) {
    cerr << "H5Tset_size failed\n";
    exit(1);
  }
  if (H5Tset_strpad(mem_type, H5T_STR_NULLTERM) < 0) {
    cerr << "H5Tset_strpad failed\n";
    exit(1);
  }
  if (H5Tset_cset(mem_type, cset) < 0) {
    cerr << "H5Tset_cset failed\n";
    exit(1);
  }

  // We need to add one ourselves for a null
  auto storage_size = H5Aget_storage_size(attr) + 1;
  char* buffer = (char*)malloc(storage_size * sizeof(char));
  memset(buffer, 0, storage_size);

  // The API for variable-length and fixed-length attributes is
  // different, so first we determine which one we are reading
  bool variable_length = H5Tequal(mem_type, attr_type);
  if (variable_length) {
    if (H5Aread(attr, mem_type, &buffer) < 0) {
      cerr << "variable-length H5Aread failed for " << name << endl;
      exit(1);
    }
  } else {
    auto fixed_type = H5Tcopy(attr_type);
    if (H5Aread(attr, fixed_type, buffer) < 0) {
      cerr << "fixed-length H5Aread failed for " << name << endl;
      exit(1);
    }
    H5Tclose(fixed_type);
  }
  
  string output(buffer);
  
  free(buffer);
  H5Aclose(attr);
  H5Tclose(attr_type);
  H5Tclose(mem_type);

  return output;
}

/*
  Loads the data in row-major order.

  If the buffer has extra space beyond that needed to load the coarse channel, we zero it out.

  This also corrects for the DC spike, if needed.
*/
void H5File::loadCoarseChannel(int i, FilterbankBuffer* buffer) const {
  assert(num_timesteps <= buffer->num_timesteps);
  assert(coarse_channel_size == buffer->num_channels);
  
  // Select a hyperslab containing just the coarse channel we want
  const hsize_t offset[3] = {0, 0, unsigned(i * coarse_channel_size)};
  const hsize_t coarse_channel_dim[3] = {unsigned(num_timesteps), 1,
                                         unsigned(coarse_channel_size)};
  if (H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
                          offset, NULL, coarse_channel_dim, NULL) < 0) {
    cerr << "failed to select coarse channel hyperslab\n";
    exit(1);
  }

  // Define a memory dataspace
  hid_t memspace = H5Screate_simple(3, coarse_channel_dim, NULL);
  if (memspace == H5I_INVALID_HID) {
    cerr << "failed to create memspace\n";
    exit(1);
  }
    
  // Copy from dataspace to memspace
  if (H5Dread(dataset, H5T_NATIVE_FLOAT, memspace, dataspace, H5P_DEFAULT, buffer->data) < 0) {    
    cerr << "h5 read failed. make sure that plugin files are in the plugin directory: " << H5_DEFAULT_PLUGINDIR << "\n";
    exit(1);
  }
    
  H5Sclose(memspace);

  if (num_timesteps < buffer->num_timesteps) {
    // Zero out the extra buffer space
  }
}
  
H5File::~H5File() {
  H5Sclose(dataspace);
  H5Dclose(dataset);
  H5Fclose(file);
}

