#include "gf256.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <concepts>
#include <iomanip>
#include <iostream>
#include <random>
#include <type_traits>

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
  EXPECT_TRUE(std::is_nothrow_default_constructible_v<GF>);
  EXPECT_TRUE(std::is_trivially_copyable_v<GF>);
  EXPECT_TRUE(std::is_trivially_copy_constructible_v<GF>);
  EXPECT_TRUE(std::is_trivially_move_constructible_v<GF>);
  EXPECT_TRUE(std::is_trivially_copy_assignable_v<GF>);
  EXPECT_TRUE(std::is_trivially_move_assignable_v<GF>);
  EXPECT_TRUE(std::is_trivially_destructible_v<GF>);
  EXPECT_EQ(sizeof(GF), 1);

  // Default constructor should initialize to zero.
  GF x;
  EXPECT_FALSE(x);
  EXPECT_TRUE(x.is_zero());
  EXPECT_EQ(x, GF());
  EXPECT_EQ(x, GF(0));
  EXPECT_EQ(x, GF(std::byte()));
  EXPECT_EQ(x, GF(std::byte(0)));
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

TEST(GF256, Divide) {
  for (GF a;; ++a.bits) {
    for (GF b(GF::max); b; --b.bits) {
      const GF c = a / b;
      EXPECT_EQ(c * b, a);
    }
    if (a == GF(GF::max)) break;
  }
}

TEST(GF256, Distribute) {
  for (GF a(GF::max); a; --a.bits) {
    for (GF b(GF::max); b; --b.bits) {
      for (GF c = b; c; --c.bits) {
        EXPECT_EQ(a * (b + c), a * b + a * c);
      }
    }
  }
}

TEST(GF256, Interpolate) {
  std::random_device dev;
  std::mt19937 rng(dev());
  std::uniform_int_distribution<int> dist(0, 255);

  // Start with 3 random shares.
  std::vector<Share> shares;
  for (int i = 1; i <= 3; ++i) {
    Share& s = shares.emplace_back();
    s.x = GF(i);
    s.ys.resize(16);
    for (GF& y : s.ys) {
      y = GF(dist(rng));
    }

    // std::cout << "Input:  " << s << std::endl;
  }

  // Create a few other shares by interpolating the first three shares.
  for (int i = 240; i <= 255; ++i) {
    Share r = interpolate({shares.cbegin(), 3}, GF(i));
    // std::cout << "Output: " << r << std::endl;
    shares.push_back(std::move(r));
  }

  const int n = shares.size();

  // Interpolating any three shares should compute all the other ones.
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      for (int k = j + 1; k < n; ++k) {
        const Share in[3] = {shares[i], shares[j], shares[k]};
        for (const Share& s : shares) {
          EXPECT_EQ(interpolate(in, s.x), s);
        }
      }
    }
  }

  // Interpolating more than three shares should give the same result.
  for (int i = 4; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      EXPECT_EQ(interpolate({shares.cbegin(), size_t(i)}, shares[j].x),
                shares[j]);
    }
  }

  // Interpolating too few shares should throw an exception.
  EXPECT_THROW(interpolate({}, GF()), std::runtime_error);
  EXPECT_THROW(interpolate({shares.cbegin(), 1}, GF(255)), std::runtime_error);

  // Interpolating duplicated shares should throw an exception.
  {
    const Share in[3] = {shares[0], shares[1], shares[0]};
    EXPECT_THROW(interpolate(in, GF(255)), std::runtime_error);
  }
}
