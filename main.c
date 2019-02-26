/*
		Autores: Alaf do Nascimento Santos e Igor Batista Vieira.
		Algoritmos Numéricos - Engenharia Elétrica UFES	
																	*/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

void main()
{
	int n = 0, numOp = 0;
	double dA, dB, DP, dC, dD, e_max = 0, e_medio = 0, tol = pow(10,-10);
	printf("n: "); scanf("%ld", &n);		
	printf("dA: "); scanf("%lf", &dA);
	printf("dB ");  scanf("%lf", &dB);
	printf("DP: "); scanf("%lf", &DP);
	printf("dC: "); scanf("%lf", &dC);
	printf("dD: "); scanf("%lf", &dD);


	double  **A = calloc(n, sizeof(double*)); //aloca as linhas da matriz;
			for(int i = 0; i < n; i++)
				A[i] = (double*) calloc(n, sizeof(double)); //aloca as colunas da matriz
	double  **A2 = calloc(n, sizeof(double*)); //aloca as linhas da matriz;
			for(int i = 0; i < n; i++)
				A2[i] = (double*) calloc(n, sizeof(double)); //aloca as colunas da matriz

	double *b, *b2, *e, *x, *x_anterior;
	b = (double*) calloc(n, sizeof(double));
	b2 = (double*) calloc(n, sizeof(double));
	e = (double*) calloc(n, sizeof(double));
	x = (double*) calloc(n, sizeof(double));
	x_anterior = (double*) calloc(n, sizeof(double));

	/*CRIA A MATRIA "A"*/
	for(int i = 0; i < n; i++)
	{
		for(int j = 0; j < n; j++)
		{
			if((i-j) == 2) A[i][j] = dA;
			else if((i-j) == 1) A[i][j] = dB;
			else if((i-j) == 0) A[i][j] = DP;
			else if((i-j) == -1) A[i][j] = dC;
			else if((i-j) == -2) A[i][j] = dD;
		}
	}

	/*CRIA O VETOR "b"*/
	for(int i = 0; i < n; i++) b[i] = 0;
	for(int i = 0; i < n; i++)
		for(int j = 0; j < n; j++)
			if(fabs(i-j) <= 2) //considera a matriz como pentadiagonal
				b[i] += A[i][j];


	/*Faz uma cópia de A e b para usar no método de seidel*/
	for(int i = 0; i < n; i++)
	{
		b2[i] = b[i];
		for(int j = 0; j < n; j++)
			A2[i][j] = A[i][j];
	}

	/* IMPRIME O SISTEMA DADO*/
	printf("\n\n[A|b] = \n\n");
    for(int i = 0; i < n; i++)
    {
      for(int j = 0;j < n; j++) printf("  %.1f ", A[i][j]);
      printf(" | %.1lf\n", b[i]);
    }

	printf("\n\nELIMINAÇÃO DE GAUSS\n");
	for(int k = 0; k < (n-1); k++)
	{
		int maior = fabs(A[k][k]), linhaMaior = k;
		/*PIVOTEAMENTO*/
		for(int i = k+1; i < (k+3); i++)
		{
			if(fabs(A[i][k]) > maior)
			{
				maior = fabs(A[i][k]);
				linhaMaior = i;
			}
			if(k == n-2)
				break;
		}
		
		if(linhaMaior != k)
		{
			for(int j = 0; j < n; j++)
			{
				double auxiliar = A[k][j];
				A[k][j] = A[linhaMaior][j];
				A[linhaMaior][j] = auxiliar;
			}
			double bAuxiliar = b[k];
			b[k] = b[linhaMaior];
			b[linhaMaior] = bAuxiliar;
		}//FIM PIVOTEAMENTO
		
		/*ESCALONAMENTO*/
		unsigned int contador = 0;
		for(int i = k+1; i < k+3; i++)
		{
		    if(contador >= 2) //matriz pentadiagonal soh tem no máximo dois multiplicadores nao nulos para cada i
		        break;
		    else
		        contador++;
			double m = A[i][k]/A[k][k];
			numOp++;
			A[i][k] = 0;
			for(int j = k+1; j < k+2; j++) 
			{
					A[i][j] -= m*A[k][j];
					numOp+=2;
			}
			b[i] -= m*b[k];
			numOp+=2;

			if(k == n-2)
				break;
		}
		//FIM ESCALONAMENTO
	}

	/* RETRO-SUBSTITUIÇÃO*/
	x[n-1]= b[n-1]/A[n-1][n-1];
	numOp++;
    for(int i = (n-2); i >= 0; i--)
    {
      double soma = b[i];
      for(int j = i+1; j < n;j++)
      {
      	if(i == 0 && j > 2)
      		break;
      	else if(abs(i - j) != 1)
      		continue;
      	soma -= A[i][j]*x[j];
      	numOp+=2;
      }
      x[i]= soma/A[i][i];
      numOp++;
    }

	/*CALCULO DO VETOR DE ERROS*/
    for(int i = 0; i < n; i++)
    	e[i] = fabs(x[i] - 1);

    for(int i = 0; i < n; i++)
    {
    	e_medio += e[i];
    	if(e_max < e[i])
    		e_max = e[i];
    } e_medio = e_medio/n;
	printf("Número de Operações: %d\n", numOp);
	printf("Erros:\n\te_max = %e\n\te_medio = %e\n\n", e_max, e_medio);
	printf("Solução:\n");
		for(int i = 0; i < n; i++)
			printf("x[%d] = \t%lf\n", i+1, x[i]);
	numOp = 0; e_max = 0; e_medio = 0;

	printf("\n\nMÉTODO DE GAUSS-SEIDEL\n");
	/*Cria Vetor Chute inicial*/
	for(int i = 0; i < n; i++)
		x_anterior[i] = b2[i]/A2[i][i];

	unsigned int iteracao = 1;
	for(;;)
	{
        /*Seidel p/ Pentadiagonal*/
        for(int i = 0; i < n; i++)
        {
            x[i] = b2[i];
            for(int j = 0; j < n; j++)
            {
                if(j < (i-2))
                    continue; 
                else if(j > (i+2))
                    break;
                else if(j > i)
                    x[i] -= A2[i][j]*x_anterior[j];
                else if(i > j)
                    x[i] -= A2[i][j]*x[j];  
                if(i != j)
                    numOp+=2;          
            }
            x[i] = x[i]/A2[i][i];
            numOp++;
        }

       	/* VERIFICANDO CRITÉRIOS DE PARADA */
       	double maiorBaixo = 0, maiorCima = 0, difRel = 0;
       	for(int i = 0; i < n; i++)
       	{
			if(maiorBaixo < x[i])
				maiorBaixo = x[i];
			if(maiorCima < fabs(x_anterior[i]-x[i]))
				maiorCima = fabs(x_anterior[i]-x[i]);
		}

		difRel = maiorCima/maiorBaixo;
		if(difRel < tol) 
			break;
		else iteracao++;

		for(int i = 0; i < n; i++)
      		x_anterior[i] = x[i];
	 }

	/*CALCULO DO VETOR DE ERROS*/
    for(int i = 0; i < n; i++)
    	e[i] = fabs(x[i] - 1);

    for(int i = 0; i < n; i++)
    {
    	e_medio += e[i];
    	if(e_max < e[i])
    		e_max = e[i];
    } e_medio = e_medio/n;
	printf("Número de Operações: %d\n", numOp);
	printf("Quantidade de Iterações: %d\n", iteracao);
	printf("Erros:\n\te_max = %e\n\te_medio = %e\n\n", e_max, e_medio);
	printf("Solução:\n");
		for(int i = 0; i < n; i++)
			printf("x[%d] = \t%lf\n", i+1, x[i]);
}
