# Galois Field GF(256)

C++20 implementation of the [Galois
Field](https://en.wikipedia.org/wiki/Finite_field_arithmetic) GF(256), aka
GF($2^8$), using the reducing polynomial $x^8 + x^4 + x^3 + x + 1$.

In this field, each element is represented by a single byte.

The provided operations are: addition, subtraction, multiplication, division,
logarithm, power, inverse and polynomial interpolation.
