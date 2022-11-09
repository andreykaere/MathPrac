#ifndef __fraction__
#define __fraction__

int gcd(int a, int b);

struct Frac
{
	int num, den;
	
	Frac()
	{
		num = 1;
		den = 1;
	}
	Frac(int p, int q)
	{
		num = p;
		den = q;
	}
	Frac& operator=(Frac q)
	{
		num = q.num;
		den = q.den;
		return *this;
	}

	Frac& operator+=(Frac q)
	{
		num = q.num * den + num * q.den;
		den = q.den * den;
        int c = gcd(num, den);
        num = num / c;
        den = den / c;
		return *this;
	}
};

struct triple
{
    Frac x;
    Frac y;
    Frac z;
};

const Frac operator+(const Frac& q1, const Frac& q2);
const Frac operator*(const int &l, const Frac &q);
const Frac operator-(const Frac& q1);
const Frac operator-(const Frac& q1, const Frac& q2);
const Frac operator*(const Frac& q1, const Frac& q2);
const Frac operator/(const Frac& q1, const Frac& q2);

#endif
