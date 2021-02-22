#include "class.h"

rho::rho(double e, double d, int j)
	{
		Emax = e;
		deltaE = d;
		Jmax = j;
		if (deltaE != 0) 
			size = int(Emax / deltaE) + 1;
		else size = int(Emax / 0.2) + 1;
		
		counter = new int* [Jmax]; //counter = new int[Jmax][size];
		for (int i = 0; i < Jmax; ++i) {
			counter[i] = new int[size];
		}
		for (int i = 0; i < Jmax; i++) 
		{
			for (int j=0; j < size; j++) 
			{
				counter[i][j] = 0;

			}
		}
	}
	void rho::Add(int J, double E)
	{
		if((J<Jmax)&&(E<Emax))counter[J][(int)(E/deltaE)]++;
	}
	int rho::read(int J, int E)
	{
		return counter[J][E];
	}
	double rho::SumJ(int index)
	{
		int res = 0;
		if (index < size)
		{
			for (int i = 0; i < Jmax; i++)
			{
				res += counter[i][index];
			}
			return res / (deltaE);
		}
		
		else
			for (int i = 0; i < Jmax; i++)
				for(int j=0;j<size;j++)
					res += counter[i][j];
		return res / deltaE;
	}

	int rho::total()
	{
		int res = 0;
		for (int j = 0; j < Jmax; j++)
			for (int E = 0; E < size; E++)
				if (counter[j][E] != 0)res+= counter[j][E];
		return res;

	}

	B::B(double e, double d, int j, double emin)
	{
		Emax = e;
		deltaE = d;
		Jmax = j;
		Emin = emin;
		if (deltaE != 0) size = int(Emax / deltaE) + 1;
		else size = int(Emax / 0.2) + 1;
		value = new double** [Jmax]; //counter = new int[Jmax][size];
		for (int i = 0; i < Jmax; ++i) {
			value[i] = new double*[size];
			for(int j=0;j<size;j++)
				value[i][j] = new double[size];
		}
		counter = new int** [Jmax]; //counter = new int[Jmax][size];
		for (int i = 0; i < Jmax; ++i) {
			counter[i] = new int* [size];
			for (int j = 0; j < size; j++)
				counter[i][j] = new int[size];
		}
		for (int i = 0; i < Jmax; i++) 
		{
			for (int j=0; j < size; j++)
			{
				for (int k = 0; k < size; k++)
				{
					value[i][j][k] = 0;
					counter[i][j][k] = 0;
				}
				

			}
		}
	}
	int B::Add(int J, double E,double v, double Ei)
	{
		if ((J>=0)&&(J < Jmax) && (E < Emax) && (E >= 0) && (Ei >= 0) && (Ei < Emax))  
		{
			value[J][(int)(E / deltaE)][(int)(Ei / deltaE)] += v;
			counter[J][(int)(E / deltaE)][(int)(Ei / deltaE)]++;
		}
		else
			return 1;
		return 0;

	}
	double B::read(int J,int E,int Ei)
	{
		return value[J][E][Ei];
	}
	double B::Avg(int index, double Eimin = 0, double Eimax = 0)
	{
		int index_max, index_min;
		if (Eimax == 0)index_max = size;
		else index_max = (int)(Eimax / deltaE);
		index_min = (int)(Eimin / deltaE);
		int sum = 0;
		double res = 0;
		for (int j = 0; j < Jmax; j++)
			for (int k = index_min; k < index_max; k++)
				if (counter[j][index][k] != 0)
				{
					res += value[j][index][k]/ counter[j][index][k];
					sum ++ ;
				}
		if (sum != 0)return res / sum;
		else
			return 0;
	}
	double B::f(int index, rho R, int& count ,double Eimin = 0, double Eimax = 0)
	{
		int index_max, index_min;
		if (Eimax == 0)index_max = size;
		else index_max = (int)(Eimax / deltaE);
		index_min = (int)(Eimin / deltaE);
		int sum = 0;
		double res = 0;
		for (int j = 0; j < Jmax; j++)
			for (int k = index_min; k < index_max; k++)
				if (counter[j][index][k] != 0)
				{
					res += value[j][index][k]* R.read(j, k) / deltaE;
					sum += counter[j][index][k];
				}
		count = sum;
		if (sum != 0)return res / sum;
		else
			return 0;
	}
