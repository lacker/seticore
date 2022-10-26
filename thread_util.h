#pragma once

#include <boost/lockfree/queue.hpp>
#include <functional>
#include <vector>

using namespace std;

class TaskList {
 public:
  vector<function<bool()> > tasks;
  boost::lockfree::queue<int> indexes;
  atomic<bool> error;
  
 TaskList(vector<function<bool()> > tasks);
  ~TaskList() {};
  void run(int num_threads);

 private:
  void runTempWorkerThread();
};

// Returns whether all worker threads completed successfully.
bool runInParallel(vector<function<bool()> > tasks, int num_threads);

void setThreadName(const string& name);
