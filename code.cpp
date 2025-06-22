#define _USE_MATH_DEFINES

#include <iostream>
#include <math.h>
#include <vector>
#include <complex>
#include <algorithm>
#include <Eigen/Dense>
#include <string>
#include "matplotlibcpp.h"
#include <random>

using namespace std;
using namespace Eigen;

namespace plt = matplotlibcpp;

string in_bits = "00100000110100101000100110101010111101101000"; // 44 bits

int Tc; // time (seconds)
int fs; // sample rate (hz)
int N; // buffer size
int W;
int M; // bits cnt
double Tb;
int Tbp; // samples num per bit
double fn;
double fn1;
double fn2;

vector<double> ASK(double A1, double A2, vector<bool> bn) {
	vector<double> signal;
	int offset = 0;
	double A = A1;
	for (int n = 0; n < N; n++)
	{
		double t = (double)n / fs;
		if (n % Tbp == 0) {
			int bit = bn[offset++];
			if (bit == 0) {
				A = A1;
			}
			else {
				A = A2;
			}
		}
		signal.push_back(A * sin(2 * M_PI * fn * t));
	}

	return signal;
}

vector<double> PSK(double A1, double A2, vector<bool> bn) {
	vector<double> signal;
	int offset = 0;
	double phase = 0; // init phase

	for (int n = 0; n < N; n++) {
		double t = (double)n / fs;
		if (n % Tbp == 0) {
			int bit = bn[offset++];
			if (bit == 1) {
				phase = M_PI;
			}
			else {
				phase = 0;
			}
		}
		signal.push_back(sin(2 * M_PI * fn * t + phase));
	}

	return signal;
}

vector<double> FSK(double A1, double A2, vector<bool> bn) {
	vector<double> signal;
	int offset = 0;
	double A = A1;
	double fnk = fn1;
	for (int n = 0; n < N; n++)
	{
		double t = (double)n / fs;
		if (n % Tbp == 0) {
			int bit = bn[offset++];
			if (bit == 0) {
				fnk = fn1;
			}
			else {
				fnk = fn2;
			}
		}
		signal.push_back(sin(2 * M_PI * fnk * t));
	}

	return signal;
}

vector<double> calcIntegral(vector<double> x) {
	vector<double> p;
	// version with height limit
	for (int ofs = 0; ofs < N; ofs++)
	{
		double s = 0;
		for (int n = 0; n < Tbp; n++) {
			if (ofs + n < N) // h
				s += x[ofs + n];

			p.push_back(s);
		}
		ofs += Tbp;
	}

	return p;
}

vector<int> bitConverter(vector<double> c) {
	vector<int> b;
	double m = 0;

	double d = 0;
	for (int ofs = 0; ofs < N; ofs += Tbp)
	{
		// remove extra bit 
		if (ofs + Tbp > N) {
			break;
		}
		for (int n = ofs; n < ofs + Tbp; n++)
		{
			if (n < N)
				d += c[n];
		}

		m = d / Tbp;

		if (m > 0.4) {
			b.push_back(0);
		}
		else {
			b.push_back(1);
		}

		d = 0;
	}

	return b;
}

vector<double> askComparator(vector<double> p, double h) {
	vector<double> c;
	for (int i = 0; i < N; i++) {
		if (p[i] < h) {
			c.push_back(0.0);
		}
		else {
			c.push_back(1.0);
		}
	}

	return c;
}

vector<int> askDemodulator(vector<double> za) {
	vector<double> x;

	double A = 1;
	for (int n = 0; n < N; n++)
	{
		double t = (double)n / fs;
		x.push_back(A * sin(2 * M_PI * fn * t) * za[n]);
	}

	vector<double> p = calcIntegral(x);

	// get max value from p
	double h = ceil(*max_element(p.begin(), p.end()) / 2);

	vector<double> c = askComparator(p, h);

	vector<int> bits;
	bits = bitConverter(c);

	return bits;
}

vector<double> pskComparator(vector<double> p) {
	vector<double> c;
	for (int i = 0; i < N; i++) {
		if (p[i] > 0) {
			c.push_back(1.0);
		}
		else {
			c.push_back(0.0);
		}
	}

	return c;
}

vector<int> pskDemodulator(vector<double> zp) {
	vector<double> x;

	double A = 1;
	for (int n = 0; n < N; n++)
	{
		double t = (double)n / fs;
		x.push_back(A * sin(2 * M_PI * fn * t) * zp[n]);
	}

	vector<double> p = calcIntegral(x);

	vector<double> c = pskComparator(p);

	vector<int> bits;
	bits = bitConverter(c);

	return bits;
}

vector<double> fskSum(vector<double> p1, vector<double> p2) {
	vector<double> p;

	for (int i = 0; i < N; i++)
	{
		p.push_back(p2[i] - p1[i]);
	}

	return p;
}

vector<double> fskComparator(vector<double> p) {
	vector<double> c;

	for (int i = 0; i < N; i++) {
		if (p[i] < 0) {
			c.push_back(1.0);
		}
		else {
			c.push_back(0.0);
		}
	}

	return c;
}

vector<int> fskDemodulator(vector<double> zp) {
	vector<double> x1;
	double A = 1;
	for (int n = 0; n < N; n++)
	{
		double t = (double)n / fs;
		x1.push_back(A * sin(2 * M_PI * fn1 * t) * zp[n]);
	}

	vector<double> x2;
	for (int n = 0; n < N; n++)
	{
		double t = (double)n / fs;
		x2.push_back(A * sin(2 * M_PI * fn2 * t) * zp[n]);
	}

	vector<double> p1 = calcIntegral(x1);
	vector<double> p2 = calcIntegral(x2);

	vector<double> p = fskSum(p1, p2);

	vector<double> c = fskComparator(p);

	vector<int> bits;
	bits = bitConverter(c);

	return bits;
}

vector<bool> hammingCoder74(vector<bool> bn) {
	vector<bool> codedBits(7);

	codedBits[2] = bn[0];
	codedBits[4] = bn[1];
	codedBits[5] = bn[2];
	codedBits[6] = bn[3];

	codedBits[0] = codedBits[2] ^ codedBits[4] ^ codedBits[6];
	codedBits[1] = codedBits[2] ^ codedBits[5] ^ codedBits[6];
	codedBits[3] = codedBits[4] ^ codedBits[5] ^ codedBits[6];

	return codedBits;
}

vector<bool> hammingDecoder74(vector<bool> codedBits) {
	vector<bool> decodedBits(4);

	int S = 0;
	S = (codedBits[0] ^ (codedBits[2] ^ codedBits[4] ^ codedBits[6])) * pow(2, 0) +
		(codedBits[1] ^ (codedBits[2] ^ codedBits[5] ^ codedBits[6])) * pow(2, 1) +
		(codedBits[3] ^ (codedBits[4] ^ codedBits[5] ^ codedBits[6])) * pow(2, 2);

	cout << "\nS: " << S << endl;

	if (S)
		codedBits[S - 1] = !codedBits[S - 1];

	decodedBits[0] = codedBits[2];
	decodedBits[1] = codedBits[4];
	decodedBits[2] = codedBits[5];
	decodedBits[3] = codedBits[6];

	return decodedBits;
}

vector<bool> strToVec(string& bitString) {
	vector<bool> bits;
	for (char c : bitString)
		bits.push_back(c == '1');

	return bits;
}

string vecToStr(vector<bool> bits) {
	string bitString;
	for (bool bit : bits)
		bitString += bit ? '1' : '0';

	return bitString;
}

double getBER(vector<bool>& bn, vector<bool>& obtBits) {
	double errBits = 0;

	for (int i = 0; i < bn.size(); i++) {
		if (bn[i] != obtBits[i])
			errBits++;
	}

	return (double)errBits / bn.size();
}

void generateWhiteNoise(double a, vector<double>& signal) {
	/// white noise ///
	random_device randDev;
	mt19937 gen(randDev());
	normal_distribution<> d(0, 1);
	/// white noise ///

	vector<double> modifiedSignal;
	for (int i = 0; i < N; i++)
		signal[i] += a * d(gen);
}

void generateWhiteNoiseMod(double b, vector<double>& signal) {
	/// white noise ///
	random_device randDev;
	mt19937 gen(randDev());
	normal_distribution<> d(0, 1);
	/// white noise ///

	vector<double> modifiedSignal;
	for (int i = 0; i < N; i++)
		signal[i] *= d(gen) * exp(-b * i);
}

string transSystem(string in_bits, double alfa, double beta, int config) {
	vector<bool> bits = strToVec(in_bits);
	vector<bool> codedBits;
	vector<bool> decodedBits;
	string result;

	/// HAMMING CODER
	for (int i = 0; i < bits.size(); i += 4) {
		vector<bool> group(bits.begin() + i, bits.begin() + i + 4);
		vector<bool> codedGroup = hammingCoder74(group);

		cout << "beforeCodedGroup: ";
		for (auto bit : group) {
			cout << bit;
		}
		cout << endl;

		for (int j = 0; j < codedGroup.size(); j++) {
			codedBits.push_back(codedGroup[j]);
		}

		cout << "codedGroup: ";
		for (int j = 0; j < codedGroup.size(); j++) {
			cout << codedGroup[j];
		}
		cout << " Size: " << codedGroup.size() << endl;
	}
	cout << "codedBits: ";
	for (auto bit : codedBits) {
		cout << bit;
	}
	cout << " Size: " << codedBits.size() << endl;
	/// HAMMING CODER

	/// VARIABLES
	M = codedBits.size();
	Tc = 1;
	fs = 44100;
	N = Tc * fs;
	W = 2;
	Tb = Tc / (double)M;
	Tbp = (int)(Tb * fs);
	fn = W / Tb;
	fn1 = (W + 1) / Tb;
	fn2 = (W + 2) / Tb;
	/// VARIABLES

	/// MODULATOR
	vector<double> modulatorOut = ASK(1, 0.5, codedBits);
	/// MODULATOR

	/// POINT A
	if (config == 1) {
		generateWhiteNoise(alfa, modulatorOut);
		generateWhiteNoiseMod(beta, modulatorOut);
	}
	else {
		generateWhiteNoiseMod(beta, modulatorOut);
		generateWhiteNoise(alfa, modulatorOut);
	}
	/// POINT B

	/// DEMODULATOR
	vector<int> demodulatorOut = askDemodulator(modulatorOut);
	cout << "demodulatorOut: ";
	for (int i = 0; i < demodulatorOut.size(); i++) {
		cout << demodulatorOut[i];
	}
	cout << " Size:" << demodulatorOut.size() << endl;
	/// DEMODULATOR

	/// HAMMING DECODER
	for (int i = 0; i < demodulatorOut.size(); i += 7) {
		vector<int> group(demodulatorOut.begin() + i, demodulatorOut.begin() + i + 7);
		vector<bool> groupBool;
		for (int j = 0; j < group.size(); j++) {
			groupBool.push_back(group[j]);
		}
		vector<bool> decodedGroup = hammingDecoder74(groupBool);

		cout << "decodedGroup: ";
		for (int j = 0; j < decodedGroup.size(); j++) {
			cout << decodedGroup[j];
		}
		cout << " Size: " << decodedGroup.size() << endl;

		for (int j = 0; j < decodedGroup.size(); j++) {
			decodedBits.push_back(decodedGroup[j]);
		}
	}
	/// HAMMING DECODER

	/// BER
	double BER = getBER(bits, decodedBits);
	cout << "BER: " << BER << endl;
	/// BER

	result = vecToStr(decodedBits);
	return result;
}

int main()
{
	// bits[>42], alfa[0,1], beta[0,10], config[1/2]
	string result;
	vector<double> alfa = { 0, 0.1, 0.5, 1.0 };
	vector<double> beta = { 0, 1, 5, 10 };

	for (auto a : alfa) {
		for (auto b : beta) {
			result = transSystem(in_bits, a, b, 1);
			cout << "Result: " << result << endl;
		}
	}

}