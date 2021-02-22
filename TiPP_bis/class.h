
class rho
{
	int size = 0;
	int** counter;
	double Emax;
	double deltaE;
	int Jmax;
public:
	rho(double, double, int );
	void Add(int J, double E);
	double SumJ(int J);
	int read(int J, int E);
	int total();
};
class B
{
	int size = 0;
	double*** value;
	int*** counter;
	double Emax;
	double deltaE;
	double Emin;
	int Jmax;
public:
	B(double, double, int, double);
	int Add(int J, double dE, double v,double E);
	double read(int ,int,int);
	double Avg(int J, double, double);
	double f(int J,rho,int&, double, double);
};