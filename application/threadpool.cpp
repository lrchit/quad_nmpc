#include "threadpool.h"
constinit thread_local int ThreadPool::worker_id_ = -1;

// ThreadPool constructor
ThreadPool::ThreadPool(int num_threads) : ctr_(0) {
    for (int i = 0; i < num_threads; i++) {
        threads_.push_back(std::thread(&ThreadPool::WorkerThread, this, i));
    }
}

// ThreadPool destructor
ThreadPool::~ThreadPool() {
    {
        std::unique_lock<std::mutex> lock(m_);
        for (int i = 0; i < threads_.size(); i++) {
            queue_.push(nullptr);
        }
        cv_in_.notify_all();
    }
    for (auto& thread : threads_) {
        thread.join();
    }
}

// ThreadPool scheduler
void ThreadPool::Schedule(std::function<void()> task) {
    std::unique_lock<std::mutex> lock(m_);
    queue_.push(std::move(task));
    cv_in_.notify_one();
}

// ThreadPool worker
void ThreadPool::WorkerThread(int i) {
    worker_id_ = i;
    while (true) {
        auto task = [&]() {
            std::unique_lock<std::mutex> lock(m_);
            cv_in_.wait(lock, [&]() { return !queue_.empty(); });
            std::function<void()> task = std::move(queue_.front());
            queue_.pop();
            cv_in_.notify_one();
            return task;
        }();
        if (task == nullptr) {
            {
                std::unique_lock<std::mutex> lock(m_);
                ++ctr_;
                cv_ext_.notify_one();
            }
            break;
        }
        task();

        {
            std::unique_lock<std::mutex> lock(m_);
            ++ctr_;
            cv_ext_.notify_one();
        }
    }
}
