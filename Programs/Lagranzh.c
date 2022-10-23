#include <stdio.h>
#include <math.h>
#include <stdlib.h>

void swap(double *a, int n, int p, int q)
{
	double *a2;
	a2 = (double*)malloc(n * n * sizeof(double));
	for(int i = 0; i < n; ++i)
	{
		if((i != p))
		{
			a2[i*n + p] = a[i*n + q];
			a2[p*n + i] = a[q*n + i];
		}
		if(i != q)
		{
			a2[q*n + i] = a[p*n + i];
			a2[i*n + q] = a[i*n + p];
		}
	}
	a2[p*n + p] = a[q*n + q];
	a2[q*n + q] = a[p*n + p];
	
	a2[p*n + q] = a[q*n + p];
	a2[q*n + p] = a[p*n + q];

	for(int i = 0; i < n; ++i)
	{
		for(int j = 0; j < n; ++j)
		{
			if(!((i!= p) && (i!= q) && (j != p) && (j != q)))
				a[i*n + j] = a2[i*n + j];

			//printf("%lf ", a[i*n + j]);
		}

		//printf("\n");
	}
	return ;
}

int main()
{
	
	double *a, *a1;
	int n = 3;
	char per[3] = {'x', 'y', 'z'};

	a = (double*)malloc(n * n * sizeof(double));
	a1 = (double*)malloc(n*n * sizeof(double));
	
	printf("Введите квадратичную форму:\n");
	
	for(int i = 0; i < n; ++i)
		for(int j = 0; j < n; ++j)
			scanf("%lf", &a[i*n + j]);

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
		//printf("%d\n", nenul);
		swap(a, n, 0, nenul);
		a1[0*n + 0] = a[0 * n + 0];
		a1[0*n + 1] = 0;
		a1[0*n + 2] = 0;
		

		a1[1*n + 1] = a[1*n + 1] - (a[0*n + 1]* a[0*n + 1])/a[0*n + 0];
		a1[1*n + 2] = a[1*n + 2] - a[0*n + 1] * a[0*n + 2] / a[0*n + 0];
		a1[2*n + 2] = a[2*n + 2] - (a[0*n + 2] * a[0*n + 2]) / a[0*n + 0];
		
		a1[1*n + 0] = a1[0*n + 1];
		a1[2*n + 0] = a1[0*n + 2];
		a1[2*n + 1] = a1[1*n + 2];
		
		//printf("\n");
		for(int i = 0; i < n; ++i)
			for(int j = 0; j < n; ++j)
				a[i*n + j] = a1[i*n + j];
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
			printf("Квадратичная форма вырожденная\n");
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

		//printf("%d %d\n", i_nen, j_nen);

		swap(a, n, i_nen, 0);
		swap(a, n, j_nen, 1);
		
		a1[0*n + 0] = 2 * a[0*n + 1];
		a1[1*n + 1] = -2 * a[0*n + 1];
		a1[0*n + 2] = (a[0*n + 2] + a[1*n + 2]);
		a1[1*n + 2] = (a[0*n + 2] - a[1*n +2]);
		
		a1[2*n + 2] = a[2*n + 2];
		a1[0*n + 1] = 0;
		a1[1*n + 0] = 0;
		
		a1[2*n + 0] = a1[0*n + 2];
		a1[2*n + 1] = a1[1*n + 2];

		
		for(int i = 0; i < n; ++i)
			for(int j = 0; j < n; ++j)
				a[i*n + j] = a1[i*n + j];

		a1[0*n + 0] = a[0 * n + 0];
		a1[0*n + 1] = 0;
		a1[0*n + 2] = 0;
		

		a1[1*n + 1] = a[1*n + 1] - (a[0*n + 1]* a[0*n + 1])/a[0*n + 0];
		a1[1*n + 2] = a[1*n + 2] - a[0*n + 1] * a[0*n + 2] / a[0*n + 0];
		a1[2*n + 2] = a[2*n + 2] - (a[0*n + 2] * a[0*n + 2]) / a[0*n + 0];
		
		a1[1*n + 0] = a1[0*n + 1];
		a1[2*n + 0] = a1[0*n + 2];
		a1[2*n + 1] = a1[1*n + 2];

		for(int i = 0; i < n; ++i)
			for(int j = 0; j < n; ++j)
				a[i*n + j] = a1[i*n + j];

	}
	
	/*for(int i = 0; i < n; ++i)
	{
		for(int j = 0; j < n; ++j)
			printf("%lf ", a[i*n + j]);

		printf("\n");
	}*/
	
	
	nenul = -1;
	for(int i = 1; i < n*n; ++i)
		if(fabs(a[i*n + i]) > 0.0000000001)
		{
			nenul = i;
			break;
		}
	
	if(nenul >= 1)
	{
		swap(a, n, nenul, 1);
		a1[1*n + 1] = a[1*n + 1];
		a1[1*n + 2] = 0;
		a1[2*n + 1] = 0;
		a1[2*n + 2] = a[2*n + 2] - a[1*n + 2] * a[1*n + 2] / a[1*n + 1];
		for(int i = 1; i < n; ++i)
			for(int j = 1; j < n; ++j)
				a[i*n + j] = a1[i*n + j];

	}
	else
	{
		if(fabs(a[1*n + 2]) < 0.0000000001)
			printf("Квадратичная форма вырожденная\n");

		else
		{
			a1[1*n + 1] = 2 * a[1*n + 2];
			a1[1*n + 2] = 0;
			a1[2*n + 1] = 0;
			a1[2*n + 2] = -2 * a[1*n + 2];
			for(int i = 0; i < n; ++i)
				for(int j = 0; j < n; ++j)
					a[i*n + j] = a1[i*n + j];
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

	free(a);
	return 0;
}
