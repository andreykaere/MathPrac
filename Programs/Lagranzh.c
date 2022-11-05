#include <stdio.h>
#include <math.h>
#include <stdlib.h>

struct Frac
{
	int num, den;
};

const Frac operator+();


//Умножение матрицы А справа на В: А1 = А*В
void mat_prod_r(double *a, double *b, int n)
{
	double buf[9];
	for(int i = 0; i < n; ++i)
		for(int j = 0; j < n; ++j)
			buf[i*n + j] = a[i*n + j];


	double h = 0;
	
	for(int i = 0; i < n; ++i)
		for(int j = 0; j < n; ++j)
		{
			h = 0;
			for(int k = 0; k < n; ++k)
				h += buf[i*n + k] * b[k*n + j];
			a[i*n + j] = h;
		}
	return ;
}


//Умножение матрицы А слева на В: А1 = В*А
void mat_prod_l(double *a, double *b, int n)
{
	double buf[9];
	for(int i = 0; i < n; ++i)
		for(int j = 0; j < n; ++j)
			buf[i*n + j] = a[i*n + j];
	
	double h = 0;
	
	for(int i = 0; i < n; ++i)
		for(int j = 0; j < n; ++j)
		{
			h = 0;
			for(int k = 0; k < n; ++k)
				h += b[i*n + k] * buf[k*n + j];
			a[i*n + j] = h;
		}

	
	return ;
}


//Преобразование матрицы квадратичной формы А1 = С^T * A * C
void transf(double *a, double *tr, int n)
{
	double trt[9];
	for(int i = 0; i < n; ++i)
		for(int j = 0; j < n; ++j)
			trt[i*n + j] = tr[j*n + i];
	
	mat_prod_l(a, trt, n);
	mat_prod_r(a, tr, n);
	
	
	
	return ;
}

//Замена р-bIX и q-bIX базиснbIX векторoB
void swap(double *a, int n, int p, int q, double *tr)
{

	double c[9];
	if(p == q)
		return ;
	for(int i = 0; i < n; ++i)
	{
		for(int j = 0; j < n; ++j)
		{
			if((i != p) && (i != q) && (j != p) && (j != q))
				c[i*n + j] = (i == j);
		}
		c[i*n + p] = 0;
		c[p*n + i] = 0;
		c[q*n + i] = 0;
		c[i*n + q] = 0;
	}

	c[p*n + p] = 0;
	c[p*n + q] = 1;
	c[q*n + p] = 1;
	c[q*n + q] = 0;

	mat_prod_r(tr, c, n);
	transf(a, c, n);

	
	return ;
}

int main()
{
	
	double a[9], tr[9], c[9], a1[9];
	int n = 3;
	char per[3] = {'x', 'y', 'z'};
	

	printf("Введите квадратичную форму:\n");
	
	for(int i = 0; i < n; ++i)
		for(int j = 0; j < n; ++j)
		{
			scanf("%lf", &a[i*n + j]);
			a1[i*n + j] = a[i*n + j];
			tr[i*n + j] = (i == j);
		}

	printf("\n");

	int nenul = -1;
	for(int i = 0; i < n; ++i)
		if(fabs(a[i*n + i]) > 0.0000000001)
		{
			nenul = i;
			break;
		}
	
	if(nenul > -1)
	{
		
		swap(a, n, 0, nenul, tr);
		
		c[0*n + 0] = 1;
		c[0*n + 1] = -a[0*n + 1] / a[0*n + 0];
		c[0*n + 2] = -a[0*n + 2] / a[0*n + 0];
		c[1*n + 0] = 0;
		c[1*n + 1] = 1;
		c[1*n + 2] = 0;
		c[2*n + 0] = 0;
		c[2*n + 1] = 0;
		c[2*n + 2] = 1;
		
		mat_prod_r(tr, c, n);
		transf(a, c, n);
	}
	else
	{
		int i_nen = -1, j_nen = -1;
		for(int i = 0; i < n; ++i)
		{
			for(int j = 0; j < n; ++j)
			{
				if(fabs(a[i*n + j]) > 0.0000000001)
				{
					i_nen = i;
					j_nen = j;
					break;
				}
			}
			if(i_nen >= 0)
				break;
		}
		
		if(i_nen < 0)
		{
			//В этом случае квадратичная форма нулевая
			for(int i = 0; i < n; ++i)
			{
				for(int j = 0; j < n; ++j)
				{
					printf("%lf ", a[i*n + j]);
				}
				printf("\n");
			}
			return 0;
		}

		

		swap(a, n, i_nen, 0, tr);
		swap(a, n, j_nen, 1, tr);
		
		c[0*n + 0] = 1;
		c[0*n + 1] = 1;
		c[0*n + 2] = 0;
		c[1*n + 0] = 1;
		c[1*n + 1] = -1;
		c[1*n + 2] = 0;
		c[2*n + 0] = 0;
		c[2*n + 1] = 0;
		c[2*n + 2] = 0;
		mat_prod_r(tr, c, n);
		transf(a, c, n);
	}
	
	
	nenul = -1;
	for(int i = 1; i < n*n; ++i)
		if(fabs(a[i*n + i]) > 0.0000000001)
		{
			nenul = i;
			break;
		}
	
	if(nenul >= 1)
	{
		swap(a, n, nenul, 1, tr);
		c[0*n + 0] = 1;
		c[0*n + 1] = 0;
		c[0*n + 2] = 0;
		c[1*n + 0] = 0;
		c[1*n + 1] = 1;
		c[1*n + 2] = -a[1*n + 2]/a[1*n + 1];
		c[2*n + 0] = 0;
		c[2*n + 1] = 0;
		c[2*n + 2] = 1;

		mat_prod_r(tr, c, n);
		transf(a, c, n);
	}
	else
	{
		if(fabs(a[1*n + 2]) < 0.0000000001)
			printf("Квадратичная форма вырожденная\n"); //В этом случае дополнительный к (1, 1) минор нулевой

		else
		{

			c[0*n + 0] = 1;
			c[0*n + 1] = 0;
			c[0*n + 2] = 0;
			c[1*n + 0] = 0;
			c[1*n + 1] = 1;
			c[1*n + 2] = 1;
			c[2*n + 0] = 0;
			c[2*n + 1] = 1;
			c[2*n + 2] = -1;
			
			mat_prod_r(tr, c, n);
			transf(a, c, n);
		}
		
	}
	

	printf("Канонический вид этой квадратичной формы:\n");
	for(int i = 0; i < n; ++i)
	{
		for(int j = 0; j < n; ++j)
			printf("%lf ", a[i*n + j]);

		printf("\n");
	}
	printf("\n");

	for(int i = 0; i < n; ++i)
	{
		for(int j = i; j < n; ++j)
		{
			if(i == j)
				printf("%lf*%c^2", a[i*n +j], per[i]);
			else
				printf("%lf*%c%c", 2*a[i*n + j], per[i], per[j]);
			if(i*n + j < n*n - 1)
				printf(" + ");

		}
	}


	printf("\n");
	printf("Матрица преобразования:\n");
	for(int i = 0; i < n; ++i)
	{
		for(int j = 0; j < n; ++j)
			printf("%lf ", tr[i*n + j]);
		printf("\n");
	}

	return 0;
}
