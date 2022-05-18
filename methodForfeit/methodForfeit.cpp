#include <iostream>
using namespace std;

double f(double* x){

    return x[0] * x[0] + 7 * x[1] * x[1] - x[0] * x[1] + x[0];
}

double g(double* x)
{
    return x[0] + x[1] - 2;
}

double P(double* x, double r)
{
    return  (r / 2) * pow(g(x),2);
}

double F(double* x, double r)
{
    return f(x) + P(x, r);
}

double Ftx(double* x, double r, double t, double* d, int sizeX) {
	double* xt = new double[sizeX];
	for (int i = 0; i < sizeX; i++) {
		xt[i] = x[i] - t * d[i];
	}
	return F(xt, r);
}

double findNorm(double* grad, double numOfValues) {
	double sumSquare = 0;
	for (int i = 0; i < numOfValues; i++) {
		sumSquare += pow(grad[i], 2);
	}
	sumSquare = pow(sumSquare, 0.5);
	return pow(sumSquare, 0.5);
}

double* findGradMassiv(double* x, double r, double numOfValues, double eps = 0.00001) {
	double* ans = new double[numOfValues];
	for (int i = 0; i < numOfValues; i++) {
		double* xWithEps = new double[numOfValues];
		for (int j = 0; j < numOfValues; j++) {
			xWithEps[j] = x[j];
		}
		xWithEps[i] += eps;
		ans[i] = (F(xWithEps, r) - F(x, r)) / eps;
	}
	return ans;
}

void printAnswer(double* x, double numOfValues, string reasonForStop) {
	cout << "Программа завершила работу из-за:" << endl << reasonForStop << endl;
	cout << "Функция достигает минимума в точке:" << endl;
	cout << "x(";
	for (int i = 0; i < numOfValues; i++) {
		if (i != numOfValues - 1) {
			cout << x[i] << ";";
		}
		else {
			cout << x[i] << ")" << endl;
		}
	}
	cout << "Функция имеет значение в точке:" << endl;
	cout << f(x) << endl;
}

double* findIntervalMethodSvenna(double* x, double r, double numOfValues, double t, double h) {
	double a, b;
	//Метод Свенна
	double xt = t;
	a = xt - h;//берем точку [0;0]
	b = xt + h;
	double fxPlusT = Ftx(x, r, b, findGradMassiv(x, r, numOfValues), numOfValues), fxMinusT = Ftx(x, r, a, findGradMassiv(x, r, numOfValues), numOfValues), fx = Ftx(x, r, xt, findGradMassiv(x, r, numOfValues), numOfValues);
	if (fxMinusT <= fx && fx >= fxPlusT) {
		cout << "Функция не унимодальна" << endl;
	}
	else {
		if (fxMinusT >= fx && fx <= fxPlusT) {
			//cout << "a0 = " << a << " ; b0 = " << b << endl;
		}
		else {
			if (fxMinusT >= fx && fx >= fxPlusT) {
				a = xt;
				xt = xt + h;
			}
			else {
				h = -h;
				b = xt;
				xt = xt + h;
			}
			int k = 1;
			double x1;
			while (true) {
				x1 = xt + pow(2, k) * h;
				double deltafx = Ftx(x, r, x1, findGradMassiv(x, r, numOfValues), numOfValues) - Ftx(x, r, xt, findGradMassiv(x, r, numOfValues), numOfValues);
				if (deltafx < 0) {
					if (h > 0) {
						a = xt;
					}
					if (h < 0) {
						b = xt;
					}
					xt = x1;
					k++;
				}
				else {
					if (h > 0) {
						b = x1;
					}
					if (h < 0) {
						a = x1;
					}
					break;
				}

			}
		}
	}
	double ans[2] = { a, b };
	return ans;
}

double findMinZolotoeSech(double* x, double r, double numOfValues, double* interval, double eps1) {
	double ak, zolotoeSceh = 0.38196, bk, delta = eps1 / 10 + (eps1 / 100) * 5, eps = eps1 / 10, xmin;
	ak = interval[0];
	bk = interval[1];
	/*cout << ak << "   " << bk << endl;
	cout << findGradMassiv(x, r, numOfValues)[0] << "   " << findGradMassiv(x, r, numOfValues)[1] << endl;*/
	int k = 0;
	double yk, zk;
	double l = 2 * eps;
	double value_fYK, value_fZK;
	yk = ak + zolotoeSceh * (bk - ak);
	zk = ak + bk - yk;
	value_fYK = Ftx(x, r, yk, findGradMassiv(x, r, numOfValues), numOfValues);
	value_fZK = Ftx(x, r, zk, findGradMassiv(x, r, numOfValues), numOfValues);
	while (true)
	{

		if (value_fYK <= value_fZK) {
			bk = zk;
			zk = yk;
			yk = ak + bk - yk;
			k++;
			value_fZK = value_fYK;
			value_fYK = Ftx(x, r, yk, findGradMassiv(x, r, numOfValues), numOfValues);

		}
		else {
			ak = yk;
			yk = zk;
			zk = ak + bk - zk;
			k++;
			value_fYK = value_fZK;
			value_fZK = Ftx(x, r, zk, findGradMassiv(x, r, numOfValues), numOfValues);
		}
		if (abs(bk - ak) <= l) {
			xmin = (ak + bk) / 2;
			break;
		}

	}
	//cout << xmin << endl;
	return xmin;
}

double* findMinXFastestDescend(double* x, double r, int numOfValues, double eps, int M) {
	double eps1 = eps / 10, eps2 = eps / 10 + (eps / 100) *5;
	double* deltaX = new double[numOfValues];
	double* x1 = new double[numOfValues];
	int k = 0;
	double t = 0;
	bool flagK = false;
	while (true) {
		double* grad = findGradMassiv(x, r, numOfValues);
		if (findNorm(grad, numOfValues) < eps1) {
			return x;
		}
		if (k >= M) {
			return x;
		}
		else {
			t = findMinZolotoeSech(x, r, numOfValues, findIntervalMethodSvenna(x, r, numOfValues, t, eps/1000), eps);
			for (int i = 0; i < numOfValues; i++) {
				x1[i] = x[i] - t * grad[i];
				deltaX[i] = x1[i] - x[i];
			}
			if (findNorm(deltaX, numOfValues) < eps2 && abs(F(x1, r) - F(x, r)) < eps2) {
				if (flagK) {
					for (int i = 0; i < numOfValues; i++) {
						x[i] = x1[i];
					}
					break;
				}
				flagK = true;
			}
			else {
				flagK = false;
			}
			for (int i = 0; i < numOfValues; i++) {
				x[i] = x1[i];
			}
			k++;
		}

	}
	return x;
}

int main()
{
	setlocale(LC_ALL, "ru");
	int numOfValues;
	double eps1, eps2, eps;
	int M, N;
	cout << "Введите количество переменных в функции" << endl;
	cin >> numOfValues;
	double* x = new double[numOfValues];
	double* x1 = new double[numOfValues];
	double* deltaX = new double[numOfValues];
	double t = 0;
	double c = 0;
	for (int i = 0; i < numOfValues; i++) {
		cout << "Введите x" << i << endl;
		cin >> x[i];
	}
	cout << "Введите eps" << endl;
	cin >> eps;
	cout << "Введите максимальное число итераций для метода наискорейшего спуска" << endl;
	cin >> M;
	cout << "Введите максимальное число итераций для метода штрафов" << endl;
	cin >> N;
	double r = 1, rk = 0;
	cout << "Введите r" << endl;
	cin >> r;
	cout << "Введите число множитель штрафа" << endl;
	cin >> c;
	int k = 0;
	while (k < N){
		cout << endl;
		cout << "шаг " << k << endl;
		cout << "r: " <<endl<< r << endl;
		x = findMinXFastestDescend(x, r, numOfValues, eps, M);
		rk = r;
		r = r * c;
		k++;
		cout<<"P = " << P(x, rk) << endl;
		cout << "В точке:" << endl;
		cout << "(";
		for (int i = 0; i < numOfValues; i++) {
			cout << x[i] << "; ";
		}
		cout << ")";
		cout << endl;
		if (abs(P(x, rk)) < eps)
			break;
		cout << "Текущие значение f(x):" << endl;
		cout << f(x) << endl;
	}
	cout << "Текущие значение f(x):" << endl;
	cout << f(x) << endl;

}