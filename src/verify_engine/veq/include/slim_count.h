#ifndef VEQ_SLIM_COUNT_H
#define VEQ_SLIM_COUNT_H

#include <cstdint>
#include <cstdio>
#include <limits>

namespace veq_slim_count {
struct count_value {
  long long value;
};

inline long long add(long long a, long long b) {
  long long r;
  if (__builtin_add_overflow(a, b, &r))
    return (a < 0) ? std::numeric_limits<long long>::min()
                   : std::numeric_limits<long long>::max();
  return r;
}

inline long long sub(long long a, long long b) {
  long long r;
  if (__builtin_sub_overflow(a, b, &r))
    return (a < 0) ? std::numeric_limits<long long>::min()
                   : std::numeric_limits<long long>::max();
  return r;
}

inline long long mul(long long a, long long b) {
  long long r;
  if (__builtin_mul_overflow(a, b, &r))
    return ((a < 0) != (b < 0)) ? std::numeric_limits<long long>::min()
                                : std::numeric_limits<long long>::max();
  return r;
}
}

typedef veq_slim_count::count_value mpq_struct;
typedef mpq_struct mpq_t[1];
typedef veq_slim_count::count_value mpz_struct;
typedef mpz_struct mpz_t[1];

inline void mpq_init(mpq_t q) { q->value = 0; }
inline void mpq_clear(mpq_t) {}
inline void mpq_set(mpq_t out, mpq_t in) { out->value = in->value; }

inline void mpq_set_d(mpq_t out, double d) {
  if (d >= static_cast<double>(std::numeric_limits<long long>::max())) {
    out->value = std::numeric_limits<long long>::max();
    return;
  }
  if (d <= static_cast<double>(std::numeric_limits<long long>::min())) {
    out->value = std::numeric_limits<long long>::min();
    return;
  }
  out->value = static_cast<long long>(d);
}

inline void mpq_add(mpq_t out, mpq_t a, mpq_t b) {
  out->value = veq_slim_count::add(a->value, b->value);
}
inline void mpq_sub(mpq_t out, mpq_t a, mpq_t b) {
  out->value = veq_slim_count::sub(a->value, b->value);
}
inline void mpq_mul(mpq_t out, mpq_t a, mpq_t b) {
  out->value = veq_slim_count::mul(a->value, b->value);
}

inline int mpq_cmp(mpq_t a, mpq_t b) {
  return (a->value > b->value) ? 1 : ((a->value < b->value) ? -1 : 0);
}

inline int mpq_cmp_ui(mpq_t a, unsigned long num, unsigned long den) {
  if (den == 0) return (a->value > 0) ? 1 : ((a->value < 0) ? -1 : 0);
  long double value = static_cast<long double>(a->value) * static_cast<long double>(den);
  long double rhs = static_cast<long double>(num);
  return (value > rhs) ? 1 : ((value < rhs) ? -1 : 0);
}

inline void mpq_get_num(mpz_t num, mpq_t q) { num->value = q->value; }

inline void mpz_init(mpz_t v) { v->value = 0; }
inline long long mpz_get_si(mpz_t v) { return v->value; }

inline int mpq_out_str(std::FILE*, int, mpq_t) { return 0; }

template <typename... Args>
inline int gmp_printf(const char*, Args&&...) {
  return 0;
}

#endif
