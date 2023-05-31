#include "gf256.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <concepts>
#include <iomanip>
#include <iostream>

namespace {

// Reference implementation of multiplication in GF(256).
GF mult_slow(GF a, GF b) {
  if (!a || !b) return GF(0);

  GF p;
  while (true) {
    if (b.bits & 1) p += a;

    b.bits >>= 1;
    if (!b) return p;

    const bool carry = a.bits & (1 << 7);
    a.bits <<= 1;
    if (carry) a.bits ^= 0b11011;
  }
}

// Reference implementation of power operation in GF(256).
GF pow_slow(GF a, int b) {
  assert(a || b >= 0);

  b %= GF::max;
  if (b == 0) return GF(1);
  if (b < 0) b += GF::max;
  assert(0 < b && b < GF::max);

  while ((b & 1) == 0) {
    b >>= 1;
    a *= a;
  }

  GF p = a;
  while (b >>= 1) {
    a *= a;
    if (b & 1) p *= a;
  }

  return p;
}

}  // namespace

TEST(GF256, Zero) {
  EXPECT_TRUE(std::regular<GF>);
  EXPECT_TRUE(std::totally_ordered<GF>);
  EXPECT_TRUE(std::three_way_comparable<GF>);

  GF x;
  EXPECT_FALSE(x);
  EXPECT_TRUE(x.is_zero());
  EXPECT_EQ(x, GF());
  EXPECT_EQ(x, GF(0));
  EXPECT_EQ(x, GF(std::byte()));
  EXPECT_NE(x, GF(1));
  EXPECT_LT(x, GF(1));
}

TEST(GF256, Add) {
  for (GF a(GF::max); a; --a.bits) {
    EXPECT_TRUE(a);
    EXPECT_FALSE(a.is_zero());

    EXPECT_EQ(a + GF(0), a);
    EXPECT_EQ(a - GF(0), a);

    EXPECT_EQ(a, +a);
    EXPECT_EQ(a, -a);

    EXPECT_FALSE(a + a);
    EXPECT_FALSE(a - a);

    for (GF b(GF::max); b; --b.bits) {
      EXPECT_TRUE(b);

      const GF c = a + b;
      EXPECT_NE(c, a);
      EXPECT_NE(c, b);

      EXPECT_EQ(c, b + a);
      EXPECT_EQ(c, a - b);
      EXPECT_EQ(c, b - a);

      EXPECT_EQ(c + a, b);
      EXPECT_EQ(c + b, a);
      EXPECT_EQ(c - a, b);
      EXPECT_EQ(c - b, a);
    }
  }
}

TEST(GF256, Multiply) {
  EXPECT_EQ(GF(1) * GF(1), GF(1));
  EXPECT_EQ(GF(0) * GF(1), GF(0));
  EXPECT_EQ(GF(1) * GF(0), GF(0));
  EXPECT_EQ(GF(0) * GF(0), GF(0));

  for (GF a(GF::max); a.bits > 1; --a.bits) {
    EXPECT_TRUE(a);

    EXPECT_EQ(a * GF(1), a);
    EXPECT_EQ(GF(1) * a, a);
    EXPECT_EQ(a / a, GF(1));
    EXPECT_EQ(a / GF(1), a);

    EXPECT_FALSE(a * GF(0));
    EXPECT_FALSE(GF(0) * a);

    for (GF b(GF::max); b.bits > 1; --b.bits) {
      EXPECT_TRUE(b);

      const GF c = a * b;
      EXPECT_EQ(c, mult_slow(a, b));
      EXPECT_NE(c, a);
      EXPECT_NE(c, b);

      EXPECT_EQ(c, b * a);

      EXPECT_EQ(c / a, b);
      EXPECT_EQ(c / b, a);
    }
  }
}

TEST(GF256, Power) {
  for (GF a(GF::max); a; --a.bits) {
    EXPECT_EQ(pow(a, 0), GF(1));
    int i = 0;
    GF p(1);
    do {
      EXPECT_LT(i, GF::max);
      p *= a;
      i += 1;
      EXPECT_TRUE(p);
      EXPECT_EQ(p, pow(a, i));
      EXPECT_EQ(p * pow(a, -i), GF(1));
    } while (p != GF(1));
    EXPECT_EQ(GF::max % i, 0);
  }
}

TEST(GF256, Inverse) {
  for (GF a(GF::max); a; --a.bits) {
    const GF b = pow(a, -1);
    EXPECT_EQ(b, pow_slow(a, -1));
    EXPECT_EQ(b, GF(1) / a);
    EXPECT_EQ(a * b, GF(1));
  }
}

TEST(GF256, Division) {
  for (GF a;; ++a.bits) {
    for (GF b(GF::max); b; --b.bits) {
      const GF c = a / b;
      EXPECT_EQ(c * b, a);
    }
    if (a == GF(GF::max)) break;
  }
}

TEST(GF256, Distribution) {
  for (GF a(GF::max); a; --a.bits) {
    for (GF b(GF::max); b; --b.bits) {
      for (GF c = b; c; --c.bits) {
        EXPECT_EQ(a * (b + c), a * b + a * c);
      }
    }
  }
}
