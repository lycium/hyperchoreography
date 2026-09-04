// A worker pool that outlives the loop it serves: one thread creation a proof, not one per
// integration step.  Free on a desktop; in wasm a thread is a Web Worker and the pool ran dry.
#pragma once
#include <condition_variable>
#include <cstdint>
#include <functional>
#include <mutex>
#include <thread>
#include <vector>

struct Pool {
  explicit Pool(int n) : n_(n < 1 ? 1 : n) {
    for (int i = 1; i < n_; i++) w_.emplace_back([this, i] { serve(i); });
  }
  ~Pool() {
    { std::lock_guard<std::mutex> lk(m_); quit_ = true; }
    cv_.notify_all();
    for (auto& t : w_) if (t.joinable()) t.join();
  }
  Pool(const Pool&) = delete;
  Pool& operator=(const Pool&) = delete;
  int size() const { return n_; }

  // body(k) for every k in [0, size()), the caller taking k = 0 rather than sitting idle.
  // Returns when all of them are done.
  void run(const std::function<void(int)>& body) {
    if (w_.empty()) { body(0); return; }
    { std::lock_guard<std::mutex> lk(m_); body_ = &body; left_ = (int)w_.size(); ++epoch_; }
    cv_.notify_all();
    body(0);
    std::unique_lock<std::mutex> lk(m_);
    done_.wait(lk, [this] { return left_ == 0; });
    body_ = nullptr;
  }

 private:
  void serve(int k) {
    std::uint64_t seen = 0;
    for (;;) {
      const std::function<void(int)>* job = nullptr;
      { std::unique_lock<std::mutex> lk(m_);
        cv_.wait(lk, [&] { return quit_ || epoch_ != seen; });
        if (quit_) return;
        seen = epoch_; job = body_; }
      if (job) (*job)(k);
      { std::lock_guard<std::mutex> lk(m_); if (--left_ == 0) done_.notify_one(); }
    }
  }

  const int n_;
  std::vector<std::thread> w_;
  std::mutex m_;
  std::condition_variable cv_, done_;
  const std::function<void(int)>* body_ = nullptr;
  std::uint64_t epoch_ = 0;
  int left_ = 0;
  bool quit_ = false;
};
