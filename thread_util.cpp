#include "thread_util.h"

#include <iostream>
#include <pthread.h>
#include <thread>

TaskList::TaskList(vector<function<bool()> > tasks)
  : tasks(tasks), indexes(tasks.size()) {
  for (int i = 0; i < (int) tasks.size(); ++i) {
    indexes.push(i);
  }
}

void TaskList::runTempWorkerThread() {
  setThreadName("tempworker");
  int index;
  while (indexes.pop(index)) {
    if (!tasks[index]()) {
      cerr << "worker task failed\n";
      exit(1);
    }
  }
}

void TaskList::run(int num_threads) {
  vector<thread> threads;
  for (int i = 0; i < num_threads; ++i) {
    threads.emplace_back(&TaskList::runTempWorkerThread, this);
  }
  for (auto& t : threads) {
    t.join();
  }
}

void runInParallel(vector<function<bool()> > tasks, int num_threads) {
  TaskList task_list(move(tasks));
  task_list.run(num_threads);
}

void setThreadName(const string& name) {
  auto p = pthread_self();
  pthread_setname_np(p, name.c_str());
}
