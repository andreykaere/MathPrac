#ifndef __fraction__
#define __fraction__

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
		num += q.num;
		den += q.den;
		return *this;
	}
};

int gcd(int a, int b);
//const Frac operator+=(const Frac& q);
const Frac operator+(const Frac& q1, const Frac& q2);
const Frac operator*(const int &l, const Frac &q);
const Frac operator-(const Frac& q1);
const Frac operator-(const Frac& q1, const Frac& q2);
const Frac operator*(const Frac& q1, const Frac& q2);
const Frac operator/(const Frac& q1, const Frac& q2);

#endif
