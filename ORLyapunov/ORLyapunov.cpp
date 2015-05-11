// ORLyapunov.cpp : コンソール アプリケーションのエントリ ポイントを定義します。
//

#include "stdafx.h"


#include <iostream>
#define _USE_MATH_DEFINES
#include <math.h>
#include <fstream>

const int N = 100;
const double Mx = 32.0;
const double Mn = 0.0009765625;

using namespace std;

//ベクトルクラス
class vector3
{
public:
	double v[3];
	static int dim;
	//friend class 宣言
	friend class matrix3;

	//コンストラクタ
	vector3()
	{
		v[0] = v[1] = v[2] = 0.0;
	}
	vector3(double v0, double v1, double v2)
	{
		v[0] = v0;
		v[1] = v1;
		v[2] = v2;

	}

	//メンバ関数(引数が1つだからクラス内定義できる)
	vector3 operator+=(const vector3& a)
	{
		for (int i = 0; i < dim; i++)
			v[i] += a.v[i];
		return *this;
	}
	vector3 operator-=(const vector3& a)
	{
		for (int i = 0; i < dim; i++)
			v[i] -= a.v[i];
		return *this;
	}
};
int vector3::dim = 3;

//オペレータ関数(引数2つ)、その他関数
vector3 operator+(const vector3& a, const vector3& b)
{
	vector3 ans;
	for (int i = 0; i < a.dim; i++)
		ans.v[i] = a.v[i] + b.v[i];
	return ans;
}
vector3 operator-(const vector3& a, const vector3& b)
{
	vector3 ans;
	for (int i = 0; i < a.dim; i++)
		ans.v[i] = a.v[i] - b.v[i];
	return ans;
}
vector3 operator*(double a, const vector3& b)
{
	vector3 ans;
	for (int i = 0; i < b.dim; i++)
		ans.v[i] = a*b.v[i];
	return ans;
}
vector3 operator/(const vector3& a, double b)
{
	vector3 ans;
	for (int i = 0; i < a.dim; i++)
		ans.v[i] = a.v[i] / b;
	return ans;
}
double abs(const vector3& a) //大きさ（ノルム）
{
	double ans = 0.0;
	for (int i = 0; i < a.dim; i++)
		ans += a.v[i] * a.v[i];
	ans = sqrt(ans);
	return ans;
}
double operator*(const vector3& a, const vector3 b) //内積
{
	double ans = 0.0;
	for (int i = 0; i < a.dim; i++)
		ans += a.v[i] * b.v[i];
	return ans;
}
void gram_schmidt(vector3 e[], vector3 ed[]) //グラムシュミットの直行化法
{
	vector3 d;
	ed[0] = e[0];
	for (int n = 1; n < e[0].dim; n++)
	{
		d = vector3(0.0, 0.0, 0.0);
		for (int i = 0; i < n; i++)
			d += (ed[i] * e[n]) / (ed[i] * ed[i])*ed[i];
		ed[n] = e[n] - d;
	}
}

double O(double x)
{
	double Ox;
	if (fabs(fmod(x, 2 * Mx)) == Mx)
		Ox = -Mx;
	else if (fabs(fmod(x, 2 * Mx)) < Mx)
		Ox = fmod(x, 2 * Mx);
	else if (fmod(x, 2 * Mx) > Mx)
		Ox = fmod(x, 2 * Mx) - 2 * Mx;
	else
		Ox = 2 * Mx - fabs(fmod(x, 2 * Mx));

	return Ox;
}

double Od(double x)
{
	double Ox = 1.0;
	if (fabs(fmod(x, 64)) == Mx)
		Ox = -1.0;
	return Ox;
}

double R(double x)
{
	double Rx;
	if (fmod(x, Mn) == 0)
		Rx = x;
	else if (x > 0)
		Rx = x - fmod(x, Mn);
	else
		Rx = x - fmod(x, Mn) - Mn;
	return Rx;
}

double Rd(double x)
{
	double Rx = 0.0;
	if (fabs(fmod(x, Mn)) == 0.0)
		Rx = 1.0;
	return Rx;
}

double g(double x)
{
	double gx;
	const double K = 1.0;
	const double SIG = 10.0;
	const double e = 1.0;

	if (x < -e)
		gx = O(O(R(K*x)) + SIG);
	else if (x >= -e && x <= e)
		gx = 0.0;
	else
		gx = O(O(R(K*x)) - SIG);

	return gx;
}

double gd(double x)
{
	double gx;
	double K = 1.0;
	const double SIG = 10.0;
	double e = 1.0;
	if (x < -e)
		gx = Od(O(R(K*x)) + SIG)*Od(R(K*x))*Rd(K*x)*K;
	else if (x >= -e && x <= e)
		gx = 0.0;
	else
		gx = Od(O(R(K*x)) - SIG)*Od(R(K*x))*Rd(K*x)*K;

	return gx;
}


class matrix3
{
public:
	double m[3][3];
	static int dim;
	//friend class 宣言
	friend class vector3;

	//コンストラクタ
	matrix3()
	{
		for (int i = 0; i < matrix3::dim; i++)
		for (int j = 0; j < matrix3::dim; i++)
			m[i][j] = 0.0;
	}
	matrix3(double v00, double v01, double v02, double v10, double v11, double v12, double v20, double v21, double v22)
	{
		m[0][0] = v00; m[0][1] = v01; m[0][2] = v02;
		m[1][0] = v10; m[1][1] = v11; m[1][2] = v12;
		m[2][0] = v20; m[2][1] = v21; m[2][2] = v22;
	}
};
int matrix3::dim = 3;

//行列とベクトルのかけ算
vector3 operator*(const matrix3& a, const vector3& b)
{
	vector3 ans;
	for (int i = 0; i < ans.dim; i++)
	for (int j = 0; j < b.dim; j++)
		ans.v[i] += a.m[i][j] * b.v[j];
	return ans;
}

int main()
{
	int n;

	short s[N];
	double x1, x2, x3, x11, x22, x33;
	double a = 1.0;
	double h1[3];
	double h2[3][3];
	double h123;
	double theta1, theta2;

	theta1 = theta2 = 8.81738;
	h123 = h1[1] = h1[2] = 30.1504;

	for (int i = 0; i < 3; i++)
	for (int j = 0; j < 3; j++)
		h2[i][j] = 30.1504;

	//メッセージ入力（任意の定数）
	for (int i = 0; i < N; i++)
		s[i] = 0.0;//(i % 64) - 32;


	vector3 x[3], xd[3], l, ramda, u1(1.0, 0.0, 0.0), u2(0.0, 1.0, 0.0), u3(0.0, 0.0, 1.0);

	ofstream ofs("ORLyapunov.txt");

	for (h1[0] = -32; h1[0] < 39.9; h1[0] += 0.1)
	{
		x1 = 0.0;
		x2 = 0.0;
		x3 = 0.0;
		u1 = vector3(1.0, 0.0, 0.0);
		u1 = vector3(0.0, 1.0, 0.0);
		u1 = vector3(0.0, 0.0, 1.0);
		l = vector3(0.0, 0.0, 0.0);
		for (n = 0; n < N; n++)
		{

			x11 = O(s[n] - g(x1) + O(R(a*x3)) + theta1);

			double sigma1, sigma2;
			sigma1 = O(R(h1[0] * x1)) + O(R(h1[1] * x2)) + O(R(h1[2] * x3));
			sigma2 = O(R(h2[0][0] * O(R(x1 * x1)))) +
				O(R(h2[0][1] * O(R(x1 * x2)))) +
				O(R(h2[0][2] * O(R(x1 * x3)))) +
				O(R(h2[1][0] * O(R(x2 * x1)))) +
				O(R(h2[1][1] * O(R(x2 * x2)))) +
				O(R(h2[1][2] * O(R(x2 * x3)))) +
				O(R(h2[2][0] * O(R(x3 * x1)))) +
				O(R(h2[2][1] * O(R(x3 * x2)))) +
				O(R(h2[2][2] * O(R(x3 * x3))));

			x22 = O(theta2 + sigma1 + O(R(O(R(O(R(h123*x11))*x2))*x3)) + sigma2);

			x33 = x2;

			matrix3 j(Od(s[n] - g(x1) + O(R(a*x3)) + theta1)*(-gd(x1)),
				0.0,
				Od(s[n] - g(x1) + O(R(a*x3)) + theta1)*Od(R(a*x3))*Rd(a*x3)*a,
				Od(theta2 + sigma1 + O(R(O(R(O(R(h123*x1))*x2))*x3)) + sigma2)
				*(Od(R(h1[0] * x1))*Rd(h1[0] * x1)*h1[0]
				+ Od(R(O(R(O(R(h123*x1))*x2))*x3))*Rd(O(R(O(R(h123*x1))*x2))*x3)*Od(R(O(R(h123*x1))*x2))*Rd(O(R(h123*x1))*x2)*Od(R(h123*x1))*Rd(h123*x1)*h123*x2*x3
				+ Od(R(h2[0][0] * O(R(x1 * x1))))*Rd(h2[0][0] * O(R(x1 * x1)))*Od(R(x1 * x1))*Rd(x1 * x1)*h2[0][0] * 2 * x1
				+ Od(R(h2[0][1] * O(R(x1 * x2))))*Rd(h2[0][1] * O(R(x1 * x2)))*Od(R(x1 * x2))*Rd(x1 * x2)*h2[0][1] * x2
				+ Od(R(h2[0][2] * O(R(x1 * x3))))*Rd(h2[0][2] * O(R(x1 * x3)))*Od(R(x1 * x3))*Rd(x1 * x3)*h2[0][2] * x3
				+ Od(R(h2[1][0] * O(R(x2 * x1))))*Rd(h2[1][0] * O(R(x2 * x1)))*Od(R(x2 * x1))*Rd(x2 * x1)*h2[1][0] * x2
				+ Od(R(h2[2][0] * O(R(x3 * x1))))*Rd(h2[2][0] * O(R(x3 * x1)))*Od(R(x3 * x1))*Rd(x3 * x1)*h2[2][0] * x3),
				Od(theta2 + sigma1 + O(R(O(R(O(R(h123*x1))*x2))*x3)) + sigma2)
				*(Od(R(h1[1] * x2))*Rd(h1[1] * x2)*h1[1]
				+ Od(R(O(R(O(R(h123*x1))*x2))*x3))*Rd(O(R(O(R(h123*x1))*x2))*x3)*Od(R(O(R(h123*x1))*x2))*Rd(O(R(h123*x1))*x2)*O(R(h123*x1))*x3
				+ Od(R(h2[0][1] * O(R(x1 * x2))))*Rd(h2[0][1] * O(R(x1 * x2)))*Od(R(x1 * x2))*Rd(x1 * x2)*h2[0][1] * x1
				+ Od(R(h2[1][0] * O(R(x2 * x1))))*Rd(h2[1][0] * O(R(x2 * x1)))*Od(R(x2 * x1))*Rd(x2 * x1)*h2[1][0] * x1
				+ Od(R(h2[1][1] * O(R(x2 * x2))))*Rd(h2[1][1] * O(R(x2 * x2)))*Od(R(x2 * x2))*Rd(x2 * x2)*h2[1][1] * 2 * x2
				+ Od(R(h2[1][2] * O(R(x2 * x3))))*Rd(h2[1][2] * O(R(x2 * x3)))*Od(R(x2 * x3))*Rd(x2 * x3)*h2[1][2] * x3
				+ Od(R(h2[2][1] * O(R(x3 * x2))))*Rd(h2[2][1] * O(R(x3 * x2)))*Od(R(x3 * x2))*Rd(x3 * x2)*h2[2][1] * x3),
				Od(theta2 + sigma1 + O(R(O(R(O(R(h123*x1))*x2))*x3)) + sigma2)
				*(Od(R(h1[2] * x3))*Rd(h1[2] * x3)*h1[2]
				+ Od(R(O(R(O(R(h123*x1))*x2))*x3))*Rd(O(R(O(R(h123*x1))*x2))*x3)*O(R(O(R(h123*x1))*x2))
				+ Od(R(h2[0][2] * O(R(x1 * x3))))*Rd(h2[0][2] * O(R(x1 * x3)))*Od(R(x1 * x3))*Rd(x1 * x3)*h2[0][2] * x1
				+ Od(R(h2[1][2] * O(R(x2 * x3))))*Rd(h2[1][2] * O(R(x2 * x3)))*Od(R(x2 * x3))*Rd(x2 * x3)*h2[1][2] * x2
				+ Od(R(h2[2][0] * O(R(x3 * x1))))*Rd(h2[2][0] * O(R(x3 * x1)))*Od(R(x3 * x1))*Rd(x3 * x1)*h2[2][0] * x1
				+ Od(R(h2[2][1] * O(R(x3 * x2))))*Rd(h2[2][1] * O(R(x3 * x2)))*Od(R(x3 * x2))*Rd(x3 * x2)*h2[2][1] * x2
				+ Od(R(h2[2][2] * O(R(x3 * x3))))*Rd(h2[2][2] * O(R(x3 * x3)))*Od(R(x3 * x3))*Rd(x3 * x3)*h2[2][2] * 2 * x3),
				0.0, 1.0, 0.0);

			//            matrix3 j(Od(x1)*(-gd(x1)),0.0,Od(x3)*Od(x3)*Rd(x3)*a,
			//                        Od(x1)*(Od(x1)*Rd(x1)*h1[0]
			//                                +Od(x1)*Rd(x1)*Od(x1)*Rd(x1)*Od(x1)*Rd(x1)*h123*x2*x3
			//                                +Od(x1)*Rd(x1)*Od(x1)*Rd(x1)*h2[0][0]*2*x1
			//                                +Od(x1)*Rd(x1)*Od(x1)*Rd(x1)*h2[0][1]*x2
			//                                +Od(x1)*Rd(x1)*Od(x1)*Rd(x1)*h2[0][2]*x3
			//                                +Od(x1)*Rd(x1)*Od(x1)*Rd(x1)*h2[1][0]*x2
			//                                +Od(x1)*Rd(x1)*Od(x1)*Rd(x1)*h2[2][0]*x3),
			//                        Od(x2)*(Od(x2)*Rd(x2)*h1[1]
			//                                +Od(x2)*Rd(x2)*Od(x2)*Rd(x2)*Od(x2)*Rd(x2)*h123*x1*x3
			//                                +Od(x2)*Rd(x2)*Od(x2)*Rd(x2)*h2[0][1]*x1
			//                                +Od(x2)*Rd(x2)*Od(x2)*Rd(x2)*h2[1][0]*x1
			//                                +Od(x2)*Rd(x2)*Od(x2)*Rd(x2)*h2[1][1]*2*x2
			//                                +Od(x2)*Rd(x2)*Od(x2)*Rd(x2)*h2[1][2]*x3
			//                                +Od(x2)*Rd(x2)*Od(x2)*Rd(x2)*h2[2][1]*x3),
			//                        Od(x3)*(Od(x3)*Rd(x3)*h1[2]
			//                                +Od(x3)*Rd(x3)*Od(x3)*Rd(x3)*Od(x3)*Rd(x3)*h123*x1*x2
			//                                +Od(x3)*Rd(x3)*Od(x3)*Rd(x3)*h2[0][2]*x1
			//                                +Od(x3)*Rd(x3)*Od(x3)*Rd(x3)*h2[1][2]*x2
			//                                +Od(x3)*Rd(x3)*Od(x3)*Rd(x3)*h2[2][0]*x1
			//                                +Od(x3)*Rd(x3)*Od(x3)*Rd(x3)*h2[2][1]*x2
			//                                +Od(x3)*Rd(x3)*Od(x3)*Rd(x3)*h2[2][2]*2*x3),
			//                        0.0,1.0,0.0);


			x1 = x11;
			x2 = x22;
			x3 = x33;

			x[0] = j * u1;
			x[1] = j * u2;
			x[2] = j * u3;

			gram_schmidt(x, xd);

			cout << j.m[0][0] << ' ' << j.m[0][1] << ' ' << j.m[0][2] << endl;
			cout << j.m[1][0] << ' ' << j.m[1][1] << ' ' << j.m[1][2] << endl;
			cout << j.m[2][0] << ' ' << j.m[2][1] << ' ' << j.m[2][2] << endl << endl;


			cout << x[0].v[0] << ' ' << x[0].v[1] << ' ' << x[0].v[2] << endl;
			cout << x[1].v[0] << ' ' << x[1].v[1] << ' ' << x[1].v[2] << endl;
			cout << x[2].v[0] << ' ' << x[2].v[1] << ' ' << x[2].v[2] << endl << endl;

			cout << xd[0].v[0] << ' ' << xd[0].v[1] << ' ' << xd[0].v[2] << endl;
			cout << xd[1].v[0] << ' ' << xd[1].v[1] << ' ' << xd[1].v[2] << endl;
			cout << xd[2].v[0] << ' ' << xd[2].v[1] << ' ' << xd[2].v[2] << endl << endl;

			for (int i = 0; i < 3; i++)
			if (abs(xd[i]) == 0)
			{
				cout << h1[0] << endl << n << endl << "error" << i << endl;
				return -1;
			}

			l.v[0] += log(abs(xd[0]));
			l.v[1] += log(abs(xd[1]));
			l.v[2] += log(abs(xd[2]));
			u1 = xd[0] / abs(xd[0]);
			u2 = xd[1] / abs(xd[1]);
			u3 = xd[2] / abs(xd[2]);



		}

		ramda = l / N;
		ofs << h1[0] << ' ' << ramda.v[0] << ' ' << ramda.v[1] << ' ' << ramda.v[2] << '\n';

	}

	ofs.close();

	return 0;

}

