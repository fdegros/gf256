#pragma once

#include <cassert>
#include <compare>
#include <cstddef>
#include <cstdint>
#include <ostream>
#include <span>
#include <stdexcept>
#include <vector>

// Implements the operations of the Galois Field GF(256) using the reducing
// polynomial x^8 + x^4 + x^3 + x + 1.
//
// In this field, each element is represented by a single byte. Therefore,
// sizeof(GF) == 1.
//
// The provided operations are: addition, subtraction, multiplication, division,
// logarithm and power.
class GF {
 public:
  // An element of GF(256) is a single byte.
  using Bits = std::uint8_t;
  Bits bits;

  // Constructors.
  GF() noexcept : bits(0) {}
  explicit GF(Bits v) noexcept : bits(v) {}
  explicit GF(std::byte v) noexcept : bits(Bits(v)) {}

  // Conversion operators.
  explicit operator Bits() const noexcept { return bits; }
  explicit operator std::byte() const noexcept { return std::byte(bits); }

  // Indicates whether this element is zero.
  bool is_zero() const noexcept { return bits == 0; }

  // Indicates whether this element is not zero.
  explicit operator bool() const noexcept { return bits != 0; }

  // Comparison operator.
  friend std::strong_ordering operator<=>(GF a, GF b) = default;

  // Stream insertion operator.
  // Prints the value in hexadecimal.
  friend std::ostream& operator<<(std::ostream& out, const GF a) {
    const char c[] = "0123456789ABCDEF";
    return out << c[a.bits >> 4] << c[a.bits & 0xF];
  }

  // The opposite of an element is this element itself.
  GF operator+() const noexcept { return *this; }
  GF operator-() const noexcept { return *this; }

  // Addition and subtraction are the same operation.
  friend GF operator+(const GF a, const GF b) noexcept {
    return GF(a.bits ^ b.bits);
  }

  friend GF operator-(const GF a, const GF b) noexcept {
    return GF(a.bits ^ b.bits);
  }

  GF& operator+=(const GF b) noexcept {
    bits ^= b.bits;
    return *this;
  }

  GF& operator-=(const GF b) noexcept {
    bits ^= b.bits;
    return *this;
  }

  // Multiplication.
  friend GF operator*(const GF a, const GF b) noexcept {
    if (!a || !b) return GF(0);

    int c = int(logs[a.bits - 1]) + int(logs[b.bits - 1]);
    if (c >= max) c -= max;

    assert(c >= 0);
    assert(c < max);
    return GF(ilogs[c]);
  }

  GF& operator*=(GF b) { return operator=((*this) * b); }

  // Division.
  // Throws: std::runtime_error if b == GF(0).
  friend GF operator/(GF a, GF b) {
    if (!b) throw std::runtime_error("Cannot divide by GF(0)");
    if (!a) return GF(0);

    int c = int(logs[a.bits - 1]) - int(logs[b.bits - 1]);
    if (c < 0) c += max;

    assert(c >= 0);
    assert(c < max);
    return GF(ilogs[c]);
  }

  GF& operator/=(GF b) { return operator=(*this / b); }

  // Returns the discrete logarithm of `a`. The base of the logarithm is the
  // generator element 3.
  // Throws: std::runtime_error if a == GF(0).
  friend int log(const GF a) {
    if (!a) throw std::runtime_error("Cannot compute log(GF(0))");
    return logs[a.bits - 1];
  }

  // Returns the generator element 3 raised to the power `b`. The exponent `b`
  // can be negative.
  static GF exp(int b) noexcept {
    b %= max;
    if (b < 0) b += max;
    assert(b >= 0);
    assert(b < max);
    return GF(ilogs[b]);
  }

  // Returns `a` raised to the power `b`. The exponent `b` can be negative if
  // `a` is not zero.
  // Throws: std::runtime_error if a == GF(0) and b <= 0.
  friend GF pow(const GF a, int b) {
    if (!a) {
      if (b <= 0)
        throw std::runtime_error("Cannot compute pow(GF(0), b) for b <= 0");
      assert(b > 0);
      return GF(0);
    }

    b %= max;
    if (b == 0) return GF(1);
    b *= int(logs[a.bits - 1]);
    return GF::exp(b);
  }

  // Maximum value. This is also the size of the multiplicative group, which
  // contains all the elements except zero.
  static constexpr int max = 255;

  // Logarithm table of generator element 3.
  static constexpr Bits logs[max] = {
      0,   25,  1,   50,  2,   26,  198, 75,  199, 27,  104, 51,  238, 223, 3,
      100, 4,   224, 14,  52,  141, 129, 239, 76,  113, 8,   200, 248, 105, 28,
      193, 125, 194, 29,  181, 249, 185, 39,  106, 77,  228, 166, 114, 154, 201,
      9,   120, 101, 47,  138, 5,   33,  15,  225, 36,  18,  240, 130, 69,  53,
      147, 218, 142, 150, 143, 219, 189, 54,  208, 206, 148, 19,  92,  210, 241,
      64,  70,  131, 56,  102, 221, 253, 48,  191, 6,   139, 98,  179, 37,  226,
      152, 34,  136, 145, 16,  126, 110, 72,  195, 163, 182, 30,  66,  58,  107,
      40,  84,  250, 133, 61,  186, 43,  121, 10,  21,  155, 159, 94,  202, 78,
      212, 172, 229, 243, 115, 167, 87,  175, 88,  168, 80,  244, 234, 214, 116,
      79,  174, 233, 213, 231, 230, 173, 232, 44,  215, 117, 122, 235, 22,  11,
      245, 89,  203, 95,  176, 156, 169, 81,  160, 127, 12,  246, 111, 23,  196,
      73,  236, 216, 67,  31,  45,  164, 118, 123, 183, 204, 187, 62,  90,  251,
      96,  177, 134, 59,  82,  161, 108, 170, 85,  41,  157, 151, 178, 135, 144,
      97,  190, 220, 252, 188, 149, 207, 205, 55,  63,  91,  209, 83,  57,  132,
      60,  65,  162, 109, 71,  20,  42,  158, 93,  86,  242, 211, 171, 68,  17,
      146, 217, 35,  32,  46,  137, 180, 124, 184, 38,  119, 153, 227, 165, 103,
      74,  237, 222, 197, 49,  254, 24,  13,  99,  140, 128, 192, 247, 112, 7};

  // Inverse logarithm table of generator element 3.
  static constexpr Bits ilogs[max] = {
      1,   3,   5,   15,  17,  51,  85,  255, 26,  46,  114, 150, 161, 248,
      19,  53,  95,  225, 56,  72,  216, 115, 149, 164, 247, 2,   6,   10,
      30,  34,  102, 170, 229, 52,  92,  228, 55,  89,  235, 38,  106, 190,
      217, 112, 144, 171, 230, 49,  83,  245, 4,   12,  20,  60,  68,  204,
      79,  209, 104, 184, 211, 110, 178, 205, 76,  212, 103, 169, 224, 59,
      77,  215, 98,  166, 241, 8,   24,  40,  120, 136, 131, 158, 185, 208,
      107, 189, 220, 127, 129, 152, 179, 206, 73,  219, 118, 154, 181, 196,
      87,  249, 16,  48,  80,  240, 11,  29,  39,  105, 187, 214, 97,  163,
      254, 25,  43,  125, 135, 146, 173, 236, 47,  113, 147, 174, 233, 32,
      96,  160, 251, 22,  58,  78,  210, 109, 183, 194, 93,  231, 50,  86,
      250, 21,  63,  65,  195, 94,  226, 61,  71,  201, 64,  192, 91,  237,
      44,  116, 156, 191, 218, 117, 159, 186, 213, 100, 172, 239, 42,  126,
      130, 157, 188, 223, 122, 142, 137, 128, 155, 182, 193, 88,  232, 35,
      101, 175, 234, 37,  111, 177, 200, 67,  197, 84,  252, 31,  33,  99,
      165, 244, 7,   9,   27,  45,  119, 153, 176, 203, 70,  202, 69,  207,
      74,  222, 121, 139, 134, 145, 168, 227, 62,  66,  198, 81,  243, 14,
      18,  54,  90,  238, 41,  123, 141, 140, 143, 138, 133, 148, 167, 242,
      13,  23,  57,  75,  221, 124, 132, 151, 162, 253, 28,  36,  108, 180,
      199, 82,  246};
};

// Struct used as input and output of the `interpolate` function.
struct Share {
  GF x;
  std::vector<GF> ys;

  friend bool operator==(const Share& a, const Share& b) = default;

  friend std::ostream& operator<<(std::ostream& out, const Share& s) {
    out << "{x: " << s.x << ", ys: [";
    const char* sep = "";
    for (const GF y : s.ys) {
      out << sep << y;
      sep = " ";
    }
    return out << "]}";
  }
};

// Interpolates polynomials using the Lagrange polynomial method.
//
// The given `shares` define the polynomials to interpolate. There must be at
// least two shares. All the shares must have distinct `x` values. All the
// shares must have the same number of `y` values.
//
// Let's `n` be the number of given shares, and `m` be the number of `y` values
// of each share:
// n == shares.size()
// m == s.ys.size() for each s in shares
//
// The degree of the interpolated polynomials is `n - 1`. A different polynomial
// is computed for each set of `y` values at each position `j` in [0..m). More
// precisely, for each `j` in [0..m), there is a unique polynomial `p[j]` of
// degree `n - 1` determined by the set of points {(s.x, s.ys[j]) for `s` in
// `shares`}. For each `j` in [0..m), and for each share `s` in `shares`, the
// polynomial p[j] of degree `n - 1` verifies the equation: p[j](s.x) ==
// s.ys[j].
//
// The result `r` of this `interpolate` function is a new share containing all
// these `m` polynomials evaluated at the destination value `dest_x`:
// r.x == dest_x
// r.ys.size() == m
// r.ys[j] == p[j](dest_x) for j in [0..m)
//
// The time complexity of this method is O(n*(n+m)).
// The space complexity is just O(m) for the resulting share.
//
// Precondition: shares.size() >= 2
// Precondition: shares[i].x != shares[j].x for i != j
// Precondition: shares[i].ys.size() == shares[j].ys.size() for i != j
inline Share interpolate(std::span<const Share> shares, GF dest_x) {
  if (shares.size() < 2) {
    throw std::runtime_error("Too few shares");
  }

  // Number of y values of each share.
  const size_t m = shares.front().ys.size();

  // Logarithm of the product of (s.x - dest_x) for s in shares.
  int a = 0;

  for (const Share& s : shares) {
    if (s.ys.size() != m) {
      throw std::runtime_error(
          "All the shares must have the same number of y values");
    }

    const GF d = s.x - dest_x;
    if (!d) {
      return s;
    }

    a += log(d);
  }

  Share r;
  r.x = dest_x;
  r.ys.resize(m);

  for (const Share& s : shares) {
    // Logarithm of the Lagrange basis polynomial evaluated at dest_x.
    int b = a - log(s.x - dest_x);
    for (const Share& t : shares) {
      if (&s != &t) {
        const GF d = s.x - t.x;
        if (!d) {
          throw std::runtime_error(
              "All the shares must have distinct x values");
        }
        b -= log(d);
      }
    }

    for (size_t i = 0; i < m; ++i) {
      if (const GF y = s.ys[i]) {
        r.ys[i] += GF::exp(b + log(y));
      }
    }
  }

  return r;
}
