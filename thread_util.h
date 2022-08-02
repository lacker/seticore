#pragma once

#include <boost/lockfree/queue.hpp>
#include <functional>
#include <vector>

using namespace std;

class TaskList {
 public:
  vector<function<bool()> > tasks;
  boost::lockfree::queue<int> indexes;

 TaskList(vector<function<bool()> > tasks);
  ~TaskList() {};
  void run(int num_threads);

 private:
  void runWorkerThread();
};

void runInParallel(vector<function<bool()> > tasks, int num_threads);
