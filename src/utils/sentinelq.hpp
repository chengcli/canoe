#ifndef SENTINELQ_HPP
#define SENTINELQ_HPP

#include <stdexcept>

template <typename T>
class SentinelQ {
 public:
  SentinelQ() : next(nullptr) {}

  virtual ~SentinelQ() {
    while (next != nullptr) pop();
  }

  void push(T const& other) {
    SentinelQ* q = this;
    while (q->next != nullptr) q = q->next;
    q->next = new SentinelQ;
    q->next->data = other;
    q->next->next = nullptr;
  }

  T& top() {
    if (next != nullptr) {
      return next->data;
    }
    std::runtime_error("SentinelQ is empty");
  }

  void pop() {
    if (next != nullptr) {
      SentinelQ* q = next;
      next = q->next;
      q->next = nullptr;
      delete q;
    }
  }

  bool empty() { return next == nullptr; }

  SentinelQ* getNext() { return next; }

  T& getData() { return data; }

 protected:
  T data;
  SentinelQ* next;
};

#endif
