#!/usr/bin/env python3
# AfAME - Developed by Dr. Zakaria Dahbi. See AUTHORS.md for development credits.
"""Independent, dependency-free verifier. Never imports or executes the C++ search."""
import argparse
import itertools
import json
import math
import sys


def require(ok, message):
    if not ok:
        raise ValueError(message)


def is_prime(p):
    return p >= 2 and all(p % v for v in range(2, math.isqrt(p) + 1))


def remainder(a, b, p):
    a = a[:]
    while a and not a[-1]:
        a.pop()
    while len(a) >= len(b):
        coef = a[-1] * pow(b[-1], -1, p) % p
        shift = len(a) - len(b)
        for j in range(len(b)):
            a[shift + j] = (a[shift + j] - coef * b[j]) % p
        while a and not a[-1]:
            a.pop()
    return a


def check_modulus(modulus, p, m):
    require(len(modulus) == m + 1 and modulus[-1] == 1, "Invalid modulus degree")
    require(all(type(x) is int and 0 <= x < p for x in modulus), "Invalid modulus coefficients")
    # Independent trial-division test, not the C++ Frobenius/GCD algorithm.
    for degree in range(1, m // 2 + 1):
        for lower in itertools.product(range(p), repeat=degree):
            require(remainder(modulus, list(lower) + [1], p), "Reducible field modulus")


class Field:
    def __init__(self, p, m, modulus):
        require(type(p) is int and type(m) is int and is_prime(p) and m >= 1, "Invalid field")
        self.p, self.m, self.q = p, m, p**m
        require(self.q <= 1024, "Verifier supports field sizes up to 1024")
        self.modulus = modulus
        if m > 1:
            check_modulus(modulus, p, m)
        else:
            require(modulus == [], "Prime fields have empty modulus metadata")
        self.inverses = {}

    def digits(self, value):
        result = []
        for _ in range(self.m):
            result.append(value % self.p)
            value //= self.p
        return result

    def encode(self, coeffs):
        return sum(c * self.p**i for i, c in enumerate(coeffs))

    def sub(self, a, b):
        if self.m == 1:
            return (a - b) % self.p
        return self.encode([(x - y) % self.p for x, y in zip(self.digits(a), self.digits(b))])

    def mul(self, a, b):
        if self.m == 1:
            return a * b % self.p
        aa, bb = self.digits(a), self.digits(b)
        product = [0] * (2 * self.m - 1)
        for i, x in enumerate(aa):
            for j, y in enumerate(bb):
                product[i+j] = (product[i+j] + x*y) % self.p
        return self.encode(remainder(product, self.modulus, self.p))

    def power(self, a, exponent):
        result = 1
        while exponent:
            if exponent & 1:
                result = self.mul(result, a)
            a = self.mul(a, a)
            exponent >>= 1
        return result

    def inverse(self, a):
        require(a != 0, "Zero has no inverse")
        if a not in self.inverses:
            value = pow(a, -1, self.p) if self.m == 1 else self.power(a, self.q - 2)
            require(self.mul(a, value) == 1, "Invalid inverse")
            self.inverses[a] = value
        return self.inverses[a]

    def rank(self, matrix):
        a = [row[:] for row in matrix]
        rows, cols, r = len(a), len(a[0]), 0
        for c in range(cols):
            pivot = next((i for i in range(r, rows) if a[i][c]), None)
            if pivot is None:
                continue
            a[r], a[pivot] = a[pivot], a[r]
            inv = self.inverse(a[r][c])
            a[r] = [self.mul(x, inv) for x in a[r]]
            for i in range(r+1, rows):
                v = a[i][c]
                if v:
                    a[i] = [self.sub(x, self.mul(v, y)) for x, y in zip(a[i], a[r])]
            r += 1
            if r == rows:
                break
        return r


def verify(data, progress=False):
    require(data.get("schema") == 1, "Unsupported certificate schema")
    n, d = data["N"], data["d"]
    require(type(n) is int and 2 <= n <= 24 and type(d) is int and d >= 2, "Invalid dimensions")
    fields, matrices, bases = [], [], set()
    for item in data["fields"]:
        f = Field(item["p"], item["m"], item["modulus"])
        require(item["q"] == f.q and f.p not in bases, "Inconsistent or repeated field factor")
        bases.add(f.p)
        matrix = item["matrix"]
        require(len(matrix) == n and all(len(row) == n for row in matrix), "Invalid matrix shape")
        require(all(type(x) is int and 0 <= x < f.q for row in matrix for x in row), "Invalid field label")
        require(all(matrix[i][j] == matrix[j][i] for i in range(n) for j in range(n)), "Nonsymmetric matrix")
        fields.append(f)
        matrices.append(matrix)
    require(fields and math.prod(f.q for f in fields) == d, "Field dimensions do not multiply to d")
    cuts = failures = cost = 0
    for k in range(1, n // 2 + 1):
        for subset in itertools.combinations(range(n), k):
            complement = [i for i in range(n) if i not in subset]
            deficits = []
            for f, matrix in zip(fields, matrices):
                r = f.rank([[matrix[i][j] for j in complement] for i in subset])
                deficits.append(k - r)
            local_cost = sum(v*v for v in deficits)
            cuts += 1
            cost += local_cost
            failures += local_cost != 0
        if progress:
            print(f"k={k}: cumulative cuts={cuts}, failures={failures}", flush=True)
    result = {"ame": failures == 0, "cuts": cuts, "failing_cuts": failures, "verified_cost": cost}
    for key, value in result.items():
        require(key in data and type(data[key]) is type(value) and data[key] == value,
                f"Certificate mismatch: {key}, declared={data.get(key)}, recomputed={value}")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("certificate")
    parser.add_argument("--progress", action="store_true")
    args = parser.parse_args()
    try:
        with open(args.certificate) as stream:
            result = verify(json.load(stream), args.progress)
        print(json.dumps(result, indent=2))
        sys.exit(0 if result["ame"] else 2)
    except (ValueError, KeyError, TypeError, IndexError, OSError) as exc:
        print(f"VERIFICATION ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
