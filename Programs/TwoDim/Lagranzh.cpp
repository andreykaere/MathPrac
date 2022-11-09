#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include "fraction.h"
#include "Lagranzh.h"
    
//Умножение матрицы А справа на В: А1 = А*В
void mat_prod_r(Frac *a, Frac *b, int n)
{
	Frac buf[9];
	for(int i = 0; i < n; ++i)
		for(int j = 0; j < n; ++j)
			buf[i*n + j] = a[i*n + j];
    
	Frac h, zer(0, 1);
	
	for(int i = 0; i < n; ++i)
		for(int j = 0; j < n; ++j)
		{
			h = zer;
			for(int k = 0; k < n; ++k)
				h += buf[i*n + k] * b[k*n + j];
        
			a[i*n + j] = h;
		}
    
	return ;
}


//Умножение матрицы А слева на В: А1 = В*А
void mat_prod_l(Frac *a, Frac *b, int n)
{

	Frac buf[9];
	for(int i = 0; i < n; ++i)
		for(int j = 0; j < n; ++j)
			buf[i*n + j] = a[i*n + j];
	
	Frac h(0, 1), zer(0, 1);
	
	for(int i = 0; i < n; ++i)
		for(int j = 0; j < n; ++j)
		{
			h = zer;
			for(int k = 0; k < n; ++k)
				h += b[i*n + k] * buf[k*n + j];
			a[i*n + j] = h;
		}
	
	return ;
}


//Преобразование матрицы квадратичной формы А1 = С^T * A * C
void transf(Frac *a, Frac *tr, int n)
{
	Frac trt[9];
	for(int i = 0; i < n; ++i)
		for(int j = 0; j < n; ++j)
			trt[i*n + j] = tr[j*n + i];
	
    
    
	mat_prod_l(a, trt, n);
	mat_prod_r(a, tr, n);
	
	return ;
}

//Замена р-bIX и q-bIX базиснbIX векторoB
void swap(Frac *a, int n, int p, int q, Frac *tr)
{
	Frac c[9];
	Frac ed(1, 1), zer(0, 1);
	if(p == q)
		return ;
	for(int i = 0; i < n; ++i)
	{
		for(int j = 0; j < n; ++j)
			if((i != p) && (i != q) && (j != p) && (j != q))
            {
				int rav = (i == j);
                c[i*n + j] = Frac(rav, 1);
            }
        
		c[i*n + p] = zer;
		c[p*n + i] = zer;
		c[q*n + i] = zer;
		c[i*n + q] = zer;
        
	}	

	c[p*n + p] = zer;
	c[p*n + q] = ed;
	c[q*n + p] = ed;
	c[q*n + q] = zer;
    
	mat_prod_r(tr, c, n);
	transf(a, c, n);
	
	return ;
}

triple Lagranzh(Frac *a)
{
	Frac tr[9], c[9];
	int n = 3;
	char per[3] = {'x', 'y', 'z'};
	//int val = 0;
	
	Frac ed(1, 1), zer(0, 1);

	//printf("Введите квадратичную форму:\n");
	
	for(int i = 0; i < n; ++i)
		for(int j = 0; j < n; ++j)
		{
			//scanf("%d", &val);
			//a[i*n + j] = Frac(val, 1);
			tr[i*n + j] = Frac((i == j), 1);
		}
    
    
    
    
    printf("\n");

	int nenul = -1;
	for(int i = 0; i < n; ++i)
		if(a[i*n + i].num != 0)
		{
			nenul = i;
			break;
		}
	
	if(nenul > -1)
	{
		
		swap(a, n, 0, nenul, tr);
		c[0*n + 0] = ed;
		c[0*n + 1] = -a[0*n + 1] / a[0*n + 0];
		c[0*n + 2] = -a[0*n + 2] / a[0*n + 0];
		c[1*n + 0] = zer;
		c[1*n + 1] = ed;
		c[1*n + 2] = zer;
		c[2*n + 0] = zer;
		c[2*n + 1] = zer;
		c[2*n + 2] = ed;
		
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
				if(a[i*n + j].num != 0)
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
					printf("%d/%d ", a[i*n + j].num, a[i*n +j].den);
				}
				printf("\n");
			}
            triple ans;
            ans.x = a[0];
            ans.y = a[1*n + 1];
            ans.z = a[2*n + 2];
			return ans;
		}

		swap(a, n, i_nen, 0, tr);
		swap(a, n, j_nen, 1, tr);
		
		c[0*n + 0] = ed;
		c[0*n + 1] = ed;
		c[0*n + 2] = zer;
		c[1*n + 0] = ed;
		c[1*n + 1] = -ed;
		c[1*n + 2] = zer;
		c[2*n + 0] = zer;
		c[2*n + 1] = zer;
		c[2*n + 2] = ed;
		mat_prod_r(tr, c, n);
		transf(a, c, n);
	}
	
	nenul = -1;
	for(int i = 1; i < n; ++i)
		if(a[i*n + i].num != 0)
		{
			nenul = i;
			break;
		}
	
	if(nenul > 0)
	{

		swap(a, n, nenul, 1, tr);
		
		c[0*n + 0] = ed;
		c[0*n + 1] = zer;
		c[0*n + 2] = zer;
		c[1*n + 0] = zer;
		c[1*n + 1] = ed;
		c[1*n + 2] = -a[1*n + 2]/a[1*n + 1];
		c[2*n + 0] = zer;
		c[2*n + 1] = zer;
		c[2*n + 2] = ed;
		mat_prod_r(tr, c, n);
		transf(a, c, n);
	}
	else
	{
		if(a[1*n + 2].num == 0)
			printf("Квадратичная форма вырожденная\n"); //В этом случае дополнительный к (1, 1) минор нулевой

		else
		{
			c[0*n + 0] = ed;
			c[0*n + 1] = zer;
			c[0*n + 2] = zer;
			c[1*n + 0] = zer;
			c[1*n + 1] = ed;
			c[1*n + 2] = ed;
			c[2*n + 0] = zer;
			c[2*n + 1] = ed;
			c[2*n + 2] = -ed;
			
			mat_prod_r(tr, c, n);
			transf(a, c, n);
		}
		
	}
	

	printf("Канонический вид этой квадратичной формы:\n");
	for(int i = 0; i < n; ++i)
	{
		for(int j = 0; j < n; ++j)
			printf("%d/%d ", a[i*n + j].num, a[i*n+j].den);

		printf("\n");
	}
	printf("\n");

	for(int i = 0; i < n; ++i)
	{
		for(int j = i; j < n; ++j)
		{
			if(i == j)
				printf("%d/%d*%c^2", a[i*n +j].num, a[i*n + j].den, per[i]);
			else
			{
				Frac out = 2*a[i*n + j];
				printf("%d/%d*%c%c", out.num, out.den, per[i], per[j]);
			}
			if(i*n + j < n*n - 1)
				printf(" + ");

		}
	}


	printf("\n");
	printf("Матрица преобразования:\n");
	for(int i = 0; i < n; ++i)
	{
		for(int j = 0; j < n; ++j)
			printf("%d/%d ", tr[i*n + j].num, tr[i*n+j].den);
		printf("\n");
	}
    
    triple ans;
    ans.x = a[0*n + 0];
    ans.y = a[1*n + 1];
    ans.z = a[2*n + 2];
    
	return ans;
}
