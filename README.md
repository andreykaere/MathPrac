## About
This repository is mainly dedicated to solving the following problem: for 
$N \in \mathbb{N}$ find a natural solution, i.e. $a, b, c \in \mathbb{N}$,
such that the following equation holds:
```math
\frac{a}{b + c} + \frac{b}{a + c} + \frac{c}{a + b} = N.
```

This problem is equivalent to finding rational points with the same sign on
the cubic, defined by this equation:
```math
x^3 + y^3 + z^3 + (1 - N) (x^2 y + x^2 z + y^2 x + y^2 z + z^2 x + z^2 y) + (3 - 2 N) x y z = 0
```

To do that, we perform a projective transformation (with integer
coefficients), which maps this cubic onto cubic with equation:
```math
y^2 z = x^3 + a x z^2 + b z^3
```
also known as __Weierstrass form__.

This form is very useful when it comes to addition of the points on the cubic
curve. So, now we're looking for any rational point on this cubic and once we
find it, we add it to itself and look at the inverse image of the point. This
inverse image has a rational coordinates and if all signs are the same -- we
win and the problem is solved. If not, we continue adding point to itself.
There is a chance, however, that this point is from finite cyclic subgroup of
group of all rational points on our Weierstrass cubic. If this happens, we
have to search for point from another subgroup and start process from the top.


## Requirements
First you need to install `python3` (I think with any version of `python3` you 
will be good to go). And you will also have to install the following packages 
for `python3`, which you can install via `pip`: 
- `matplotlib`
- `sympy`
- `numpy`


## Usage
To find the solution for the problem run `Programs/main/main.py`. You can also
run the file `Programs/main/weierstrass_form/weierstrass_main.py` to find a
Weierstrass form for the given cubic (if possible).

