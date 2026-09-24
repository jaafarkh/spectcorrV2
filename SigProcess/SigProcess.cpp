// SigProcess.cpp : Defines the entry point for the console application.
// Written by Jaafar Alsalaet
//Last edited on 30-8-2021
//c++ implementation of spectral correlation/coherence
//using ACP, FastACP, FAM, FSC and fast Dirichlet kernel SC
//
//FFTW library is required to compile and run the project. Depending on your buil configuration CPU type, 
//rename either libfftw3-3-x86.dll or libfftw3-3-x64.dll to libfftw3-3.dll and put it in your application .exe folder
//for more information, see http://www.fftw.org/install/windows.html
//I use x86 build setting


#include "stdafx.h"
#include "sigpro.h"
#include <fftw3.h>
#include <windows.h>
#include <sstream>
#define USEFFTW 1
#define DbgMsg( s )            \
{                             \
std::wostringstream os_;    \
   os_ << s;                   \
   OutputDebugStringW( os_.str().c_str() );  \
}

const double pi = 3.1415926535897932384626433832795;

polar operator +(polar const& c1, polar const& c2) {
	double real, imag;
	real = c1.amp * cos(c1.angle) + c2.amp * cos(c2.angle);
	imag = c1.amp * sin(c1.angle) + c2.amp * sin(c2.angle);	
	return polar(sqrt(imag*imag + real*real), atan2(imag, real));;
}

polar operator -(polar const& c1, polar const& c2) {
	complex c3;
	c3.real = c1.amp * cos(c1.angle) - c2.amp * cos(c2.angle);
	c3.imag = c1.amp * sin(c1.angle) - c2.amp * sin(c2.angle);
	return c3.toPolar();
}

polar operator *(polar const& c1, polar const& c2) {
	polar c3;
	c3.amp = c1.amp * c2.amp;
	c3.angle = c1.angle + c2.angle;
	return c3;
}

polar operator /(polar const& c1, polar const& c2) {
	polar c3;
	c3.amp = c1.amp / c2.amp;
	c3.angle = c1.angle - c2.angle;
	return c3;
}

complex operator +(complex const& c1, complex const& c2) {
	complex c3;
	c3.real = c1.real + c2.real;
	c3.imag = c1.imag + c2.imag;
	return c3;
}

complex operator -(complex const& c1, complex const& c2) {
	complex c3;
	c3.real = c1.real - c2.real;
	c3.imag = c1.imag - c2.imag;
	return c3;
}

complex operator *(complex const& c1, complex const& c2) {
	complex c3;
	c3.real = c1.real*c2.real - c1.imag*c2.imag;
	c3.imag = c1.real*c2.imag + c2.real*c1.imag;
	return c3;
	//return complex(c1.real*c2.real - c1.imag*c2.imag, c1.real*c2.imag + c2.real*c1.imag);
}
complex operator *(float const& c1, complex const& c2) {
	complex c3;
	c3.real = c1 * c2.real;
	c3.imag = c1 * c2.imag;
	return c3;	
}
complex operator *(double const& c1, complex const& c2) {
	complex c3;
	c3.real = c1 * c2.real;
	c3.imag = c1 * c2.imag;
	return c3;
}
complex operator /(complex const& c1, complex const& c2) {
	complex c3;
	double r2;
	r2 = c2.real*c2.real + c2.imag*c2.imag;
	c3.real = (c1.real*c2.real + c1.imag*c2.imag) / r2;
	c3.imag = (c2.real*c1.imag - c1.real*c2.imag) / r2;
	return c3;
}


//define sqrt of a complex variable using relations that do not need arctan function
complex csqrt(complex z) {
	// 1. Handle the edge case of a pure zero input
	if (z.real == 0.0 && z.imag == 0.0) {
		return complex(0.0, 0.0);
	}
	// 2. Calculate the magnitude (absolute value) |z| = sqrt(x^2 + y^2)
	// Using std::hypot prevents overflow during intermediate squaring
	double r = std::hypot(z.real, z.imag);
	double u, v;
	// 3. Branching to maintain precision and correct signs
	if (z.real >= 0.0) {
		u = std::sqrt((r + z.real) * 0.5);
		v = z.imag / (2.0 * u);
	}
	else {
		v = (z.imag >= 0.0 ? 1.0 : -1.0) * std::sqrt((r - z.real) * 0.5);
		u = z.imag / (2.0 * v);
	}
	return complex(u, v);
}

void fft(complex *A, int Nb, int stat)
{
	//stat=1 forward transform, -1 backward transform
	// starting index for the input and output data stream is 0
	int lpk, l, M, me1, k, j, nbd2, nbm1, N;
	complex u1, w1;
	complex to1;
	//float uu2 ;

	N = (int)(log2((double)Nb));

	//		BitReversal(A, N); //~same speed as the following code

	nbd2 = Nb / 2;
	nbm1 = Nb - 1;
	j = 0;
	for (l = 0; l < nbm1; l++)
	{
		if (l < j)
		{
			to1 = A[j];
			A[j] = A[l];
			A[l] = to1;
		}

		k = nbd2;
	label3:
		if (k > j)
		{
			j = j + k;
		}
		else
		{
			j = j - k;
			k = k / 2;
			goto label3;
		}
	}

	for (M = 1; M <= N; M++)
	{
		u1 = complex(1.0, 0.0);
		me1 = 1 << M;
		k = me1 >> 1;
		w1 = complex(cos(pi / k), -1.0 * stat * sin(pi / k));
		for (j = 0; j < k; j++)
		{
			for (l = j; l< Nb; l += me1)
			{
				lpk = l + k;
				to1 = A[lpk] * u1;
				A[lpk] = A[l] - to1;
				A[l] = A[l] + to1;

			}//end of L			
			u1 = u1 * w1;
		} //end of j
	} //end {of m}


}
static int nextPow2(int n) {
	//int p = 1;
	//while (p < n) p <<= 1;	
	//return p;
	//fastest method using builtin c++ function
	unsigned long idx;
	_BitScanReverse(&idx, n - 1);
	return 1u << (idx + 1);
}

// FFT of 'in' (length n <= m), zero-padded to length m (power of two).
// Result written into 'out' (length m). 'in' and 'out' may NOT alias.
static void fftPad(complex* in, int n, complex* out, int m, int st) {
	for (int i = 0; i < n; i++) out[i] = in[i];
	for (int i = n; i < m; i++) out[i] = complex(0.0, 0.0);
	fft(out, m, st);
}
void ffthelper(complex* A, complex* B, int N)
{//helper function to keep the convention of fftw
	int i;
	for (i = 0; i < N; i++)
		B[i] = A[i];
	fft(B, N, 1);
}
int FSCoh(float* x, float &da, float amax, int L, int Noverlap, int Nw, int opt, float* SCo)
{


	/*
	calculate spectral correlation/coherance using Fast ACP

	da: delta alpha relative to Fs
	x() input data
	amax: maximum cyclic frequency (normalized by Fs)
	L: input data length
	to get best speed and accuracy, make sure L gives power of two for KN , where KN = (L - Nw + R) / R;=> L = R*KN + Nw - R
	Nw: window length        '
	opt:0 :SC no negative frequencies correlation, 1 : Scoh, 2: SC with negative frequencies, 3: Scoh with negative freq

	SCo() output data: of size Nw/2 * a2
	*/
	int nfft, hnfft, KN, i, k, j, R, Nb, ri, ci, mri, inc;
	bool Cohr;
	double mk, cp1;
	R = Nw - Noverlap;
	nfft = Nw;
	hnfft = nfft / 2;
	auto windat = new double[Nw];
	fftw_complex *in, *out;
	fftw_plan p;
	complex *sptr;
	//auto a = new complex[Nw];
	//auto b = new complex[Nw];
	complex* a = (complex*)fftw_malloc(sizeof(complex) * Nw); //new complex[Nw];
	complex* b = (complex*)fftw_malloc(sizeof(complex) * Nw); //new complex[Nw];
	in = (fftw_complex*)&a[0]; //(fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
	out = (fftw_complex*)&b[0];
	p = fftw_plan_dft_1d(Nw, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	//Prep WinData;
	for (i = 0; i < Nw; i++)
		windat[i] = 2.0 * (0.5 - 0.5 * cos(2 * i * pi / Nw));

	KN = (int)((L - Noverlap) / (double)(Nw - Noverlap));
	Cohr = false;
	if ((opt == 1) || (opt == 3)) Cohr = true;
	//std::cout << "KN: " << KN << "\n";	

	auto dshift = new complex[Nw];
	complex sfa, sfm, kfk;
	auto CPY2 = new double[hnfft](); //putting () will initialize all elements to zero	
	double KMU = 0.0;
	for (i = 0; i < Nw; i++)
		KMU = KMU + windat[i] * windat[i];
	KMU = KN * KMU;
	Nb = (int)round((1.0f / da / Nw));
	DbgMsg("FastACP Nb: " << Nb << "\n");
	if (Nb < 1) {
		//("Window length must be smaller than signal length");
		return -1;
	}
	da = 1.0F / (Nw * Nb);
	int a2 = (int)(amax / da);
	auto SCoCmplx = new complex[hnfft * a2];
	auto SCoX2 = new double[hnfft * a2]();
	auto Kft = new complex[nfft * Nb];

	for (i = 0; i < Nw; i++)
		dshift[i] = complex(cos(2 * pi * da * i), sin(2 * pi * da * i));
	sfm = complex(cos(2 * pi * R * da), sin(2 * pi * R * da));
	for (k = 0; k < KN; k++)
	{
		for (i = 0; i < nfft; i++)
			a[i] = complex(x[i + k * R] * windat[i], 0.0);
#ifdef USEFFTW
		fftw_execute(p);
#else
		ffthelper(a, b, Nw);
#endif // USEFFTW	

		for (i = 0; i < hnfft; i++)
		{
			Kft[i] = b[i];
			if (i > 0) Kft[nfft - i] = Kft[i].conj();
			mk = b[i].real; cp1 = b[i].imag;
			CPY2[i] = CPY2[i] + mk * mk + cp1 * cp1;
		}

		for (j = 1; j < Nb; j++) // calculate complete mesh across one deltaF
		{

			for (i = 0; i < Nw; i++)
				a[i] = a[i] * dshift[i]; // apply shifting
										 //fftw_execute(p); 
#ifdef USEFFTW
			fftw_execute(p);
#else
			ffthelper(a, b, Nw);
#endif // USEFFTW	
			inc = j * nfft;
			for (i = 0; i < nfft; i++)
				Kft[i + inc] = b[i];
		}

		for (i = 0; i < Nw; i++)
			dshift[i] = dshift[i] * sfm;//account for time of the next block

		for (i = 0; i < hnfft; i++)
		{
			sfa = Kft[i];
			mri = i;
			ci = 0;
			mk = -1.0;
			inc = 0;
			sptr = &SCoCmplx[i*a2];
			for (j = 0; j < a2; j++)
			{
				if (ci == Nb) {
					ci = 0;
					mri = mri - 1;
					sfa = sfa * complex(cos(2 * pi * (k * R) * Nb * da), sin(2 * pi * (k * R) * Nb * (mk * da)));
					if (mri == 0) {
						mk = 1.0F;
						sfa = sfa.conj();
					}
					if (mri == -1) inc = nfft;
					if ((mri < 0) && (opt == 0)) break; //prevent corr with negative spectrum
				}
				ri = inc + mri;
				kfk = complex(Kft[ri + ci * nfft].real, mk * Kft[ri + ci * nfft].imag);

				*sptr = *sptr + sfa * kfk;
				sptr++;
				ci++;
				if (Cohr)
					SCoX2[i*a2 + j] = SCoX2[i*a2 + j] + kfk.real * kfk.real + kfk.imag * kfk.imag;
				//SCoCmplx[i + j * hnfft] = SCoCmplx[i + j * hnfft] + kfk * sfa;

			}
		}
	}

	if (Cohr)
	{
		for (i = 0; i < hnfft; i++)
		{
			for (j = 0; j < a2; j++)
			{
				if (j == 0) {
					SCo[i + j * hnfft] = 0.0f; //nullify alpha=0 component to protect SCho map
				}
				else {
					SCo[i + j * hnfft] = (float)(SCoCmplx[i*a2 + j].amp() / sqrt(CPY2[i] * SCoX2[i*a2 + j])); // pick the nearest bin
				}
			}
		}
	}
	else
	{
		for (i = 0; i < hnfft; i++)
		{
			for (j = 0; j < a2; j++)
				SCo[i + j * hnfft] = (float)(SCoCmplx[i*a2 + j].amp() / KMU);
		}
	}

	fftw_destroy_plan(p);
	fftw_free(a);
	fftw_free(b);
	delete[] windat; delete[] SCoCmplx; delete[] SCoX2;
	delete[] dshift; delete[] CPY2; delete[] Kft;
	return 0;
}
int SCoh(float* x, float da, float amax, int L, int Noverlap, int Nw, int opt, float* SCo)
{
	/*
	calculate spectral correlation/coherance using ACP

	da: delta alpha relative to Fs
	x() input data
	amax: maximum cyclic frequency (normalized by Fs)
	L: input data length
	to get best speed and accuracy, make sure L gives power of two for KN , where KN = (L - Nw + R) / R;=> L = R*KN + Nw - R
	Nw: window length        '
	opt:0 :SC no negative frequencies correlation, 1 : Scoh, 2: SC with negative frequencies, 3: Scoh with negative freq

	SCo() output data: of size Nw/2 * a2
	*/
	int nfft, hnfft, KN, i, k, j, R;
	R = Nw - Noverlap;
	nfft = Nw;
	hnfft = nfft / 2;
	auto windat = new double[Nw];
	fftw_complex *in, *out;
	fftw_plan p;

	complex* a = (complex*)fftw_malloc(sizeof(complex) * Nw); //new complex[Nw];
	complex* b = (complex*)fftw_malloc(sizeof(complex) * Nw); //new complex[Nw];
	in = (fftw_complex*)&a[0]; //(fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
	out = (fftw_complex*)&b[0];
	p = fftw_plan_dft_1d(Nw, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	//Prep WinData;
	for (i = 0; i < Nw; i++)
		windat[i] = 2.0 * (0.5 - 0.5 * cos(2 * i * pi / Nw));

	KN = (int)((L - Noverlap) / (double)(Nw - Noverlap));
	//std::cout << "KN: " << KN << "\n";
	auto yws = new complex[hnfft * KN];
	auto dshift = new complex[L];
	//complex sfa, sfm;
	auto CPY2 = new double[hnfft];
	auto CPX2 = new double[hnfft];
	auto CPS = new complex[hnfft];
	auto xc = new complex[L];
	double KMU = 0.0;
	int a2 = (int)(amax / da);
	for (i = 0; i < L; i++) {
		xc[i] = complex(x[i], 0.0);
		dshift[i] = complex(cos(2 * pi * da * i), sin(2 * pi * da * i));
	}
	for (i = 0; i < Nw; i++)
		KMU = KMU + windat[i] * windat[i];
	KMU = KN * KMU;

	for (i = 0; i < hnfft; i++)
	{
		CPY2[i] = 0.0;
	}
	for (k = 0; k < KN; k++) // calculate initial spectral data
	{

		for (i = 0; i < nfft; i++)
			a[i] = complex(x[i + k * R] * windat[i], 0.0);
#ifdef USEFFTW
		fftw_execute(p);
#else
		ffthelper(a, b, Nw);
#endif // USEFFTW	

		for (i = 0; i < hnfft; i++)
		{
			yws[i + k * hnfft] = b[i];
			CPY2[i] = CPY2[i] + b[i].real * b[i].real + b[i].imag * b[i].imag;
		}
	}


	for (j = 0; j < a2; j++)
	{
		for (i = 0; i < hnfft; i++)
		{
			CPX2[i] = 0.0;
			CPS[i] = complex(0.0, 0.0);
		}
		for (k = 0; k < KN; k++) {
			for (i = 0; i < nfft; i++)
				a[i] = windat[i] * xc[i + k * R];
#ifdef USEFFTW
			fftw_execute(p);
#else
			ffthelper(a, b, Nw);
#endif // USEFFTW	
			for (i = 0; i < hnfft; i++)
			{
				if ((j * da <= ((float)i / Nw)) || (opt > 0)) {
					CPS[i] = CPS[i] + yws[i + k * hnfft] * (b[i].conj());
				}
				CPX2[i] = CPX2[i] + b[i].real * b[i].real + b[i].imag * b[i].imag;
			}
		}
		for (i = 0; i < hnfft; i++) {
			if ((opt == 0) || (opt == 2)) {
				SCo[i + j * hnfft] = (float)(CPS[i].amp() / KMU);
			}
			else {
				if (j == 0) {
					SCo[i + j * hnfft] = 0.0f; //nullify alpha=0 component to protect SCho map
				}
				else {
					SCo[i + j * hnfft] = (float)(CPS[i].amp() / sqrt((CPX2[i] * CPY2[i])));
				}
			}

		}
		for (i = 0; i < L; i++)
			xc[i] = xc[i] * dshift[i];

	}



	fftw_destroy_plan(p);
	fftw_free(a);
	fftw_free(b);
	delete[] windat;
	delete[] yws; delete[] dshift; delete[] xc; delete[] CPY2; delete[] CPX2; delete[] CPS;
	return 0;
}

//========================================================================================================================

int AntoniFSC2(float* x, float &da, float amax, int L, int Nw, int opt, float* SCo)
{
	/*
	calculate spectral correlation/coherance using Antoni et al. Fast SC
	This is two sided, p = -P ... P (according to matlab code)
	da: delta alpha relative to Fs
	x() input data
	amax: maximum cyclic frequency (normalized by Fs)
	L: input data length
	to get best speed and accuracy, make sure L gives power of two for KN , where KN = (L - Nw + R) / R;=> L = R*KN + Nw - R
	if KN is not power-of-two, the function will adjust it by zero padding and update da
	Nw: window length        '
	opt:0 :SC no negative frequencies correlation, 1 : Scoh, 2: SC with negative frequencies, 3: Scoh with negative freq
	output:
	adjusted da
	SCo() output data
	*/
	int nfft, hnfft, KN, KNP2, i, k, j, R, P, M, ir, ir2;
	complex sfi, sfj, sfj1;
	R = (int)floor(1.0f / (2.0f * amax));
	if (R > (Nw / 4))
		R = Nw / 4;
	//std::cout << "R: " << R << "\n";

	DbgMsg("AntoniFSC R: " << R << "\n");
	nfft = Nw;
	hnfft = nfft / 2;
	auto windat = new double[Nw];
	fftw_complex *in, *out, *ink, *outk;
	fftw_plan fp, fp2;

	complex* a = (complex*)fftw_malloc(sizeof(complex) * Nw); //new complex[Nw];
	complex* b = (complex*)fftw_malloc(sizeof(complex) * Nw); //new complex[Nw];
	in = (fftw_complex*)&a[0]; //(fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
	out = (fftw_complex*)&b[0];
	fp = fftw_plan_dft_1d(Nw, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	//Prep WinData;
	for (i = 0; i < Nw; i++)
		windat[i] = (0.5 - 0.5 * cos(2 * i * pi / Nw));

	//std::cout << "P: " << P << "\n";
	KN = (int)((L - Nw + R) / R);
	//KNP2 = (int)pow(2.0, ceil(log2((float)KN)));//next power of two
	KNP2 = nextPow2(KN);
	//std::cout << "KN, KNP2: " << KN<<", "<< KNP2 << "\n";
	R = (int)ceil((float)(L - Nw) / (float)(KNP2 - 1)); //updated R
	if (R < 1) R = 1;
	KN = (int)((L - Nw + R) / R); //KN will be less than KNP2 but larger than previous KN
	P = (int)floor((float)Nw / (2.0f * R));
	DbgMsg("AntoniFSC KN, KNP2: " << KN << ", " << KNP2 << "\n");
	da = 1.0F / (KNP2 * R); //updated normalized da
							//std::cout << "updated da: " << da << "\n";
	int a2 = (int)(amax / da);
	complex* Ak = (complex*)fftw_malloc(sizeof(complex) * KNP2); //new complex[Nw];
	complex* Bk = (complex*)fftw_malloc(sizeof(complex) * KNP2); //new complex[Nw];
	ink = (fftw_complex*)&Ak[0]; //create alias
	outk = (fftw_complex*)&Bk[0];
	fp2 = fftw_plan_dft_1d(KNP2, ink, outk, FFTW_FORWARD, FFTW_ESTIMATE);
	//std::cout << "KN: " << KN << "\n";
	auto yws = new complex[hnfft * KN];
	auto SCoCmplx = new complex[hnfft * a2];
	auto CPY2 = new double[hnfft];
	auto Rw = new double[KNP2]();
	auto SRw = new double[KN / 2]();
	double KMU = 0.0, RW0, w;

	for (i = 0; i < Nw; i++)
		KMU = KMU + windat[i] * windat[i];
	RW0 = KMU;
	KMU = KN * KMU;
	M = Nw / 2;
	for (i = 0; i < KNP2; i++) //generate Rw(alpha) from 0 to KN-1
	{
		for (j = 0; j < Nw; j++)
			Rw[i] = Rw[i] + windat[j] * windat[j] * cos(2.0f * pi * (j - M) * i * da);
	}

	for (i = 0; i < hnfft; i++)
	{
		CPY2[i] = 0.0;
	}
	for (k = 0; k < KN; k++) // calculate initial spectral data
	{

		for (i = 0; i < nfft; i++)
			a[i] = complex(x[i + k * R] * windat[i], 0.0);
#ifdef USEFFTW
		fftw_execute(fp);
#else
		ffthelper(a, b, Nw);
#endif // USEFFTW	

		for (i = 0; i < hnfft; i++)
		{
			yws[i + k * hnfft] = b[i];
			CPY2[i] = CPY2[i] + b[i].real * b[i].real + b[i].imag * b[i].imag;
		}
	}
	if ((opt == 1) || (opt == 3)) {
		for (k = 0; k < KN; k++) {
			for (i = 0; i < hnfft; i++)
				yws[i + k * hnfft] = yws[i + k * hnfft] / (sqrt(CPY2[i] / KMU));
		}
	}
	for (k = 0; k < KNP2; k++) {
		Ak[k] = complex(0.0f, 0.0f);
	}
	M = KN / 2;
	if (a2 < M) M = a2;
	w = 2.0f * pi * (Nw / 2 - 1);
	sfi = complex(cos(w * da), -sin(w * da));
	for (j = 0; j <= P; j++) {
		sfj1 = complex(cos(w * j / Nw), sin(w * j / Nw));
		for (k = 0; k < hnfft; k++) {
			if (k >= j) {
				for (i = 0; i < KN; i++) {
					Ak[i] = yws[k + i * hnfft] * yws[(k - j) + i * hnfft].conj();
				}
#ifdef USEFFTW
				fftw_execute(fp2);
#else
				ffthelper(Ak, Bk, KNP2);
#endif // USEFFTW	

				sfj = sfj1;
				for (i = 0; i < M; i++) {
					//Bk[i] = Bk[i] * complex(cos(w * (i * da - (float)j / Nw)), -sin(w * (i * da - (float)j / Nw)));
					SCoCmplx[k * a2 + i] = SCoCmplx[k * a2 + i] + Bk[i] * sfj;
					sfj = sfj * sfi;
				}
				//now will process p < 0
				if (j > 0) {
					for (i = 0; i < KN; i++) {
						Ak[i] = yws[k + i * hnfft].conj() * yws[(k - j) + i * hnfft];
					}
#ifdef USEFFTW
					fftw_execute(fp2);
#else
					ffthelper(Ak, Bk, KNP2);
#endif // USEFFTW	
					sfj = sfj1.conj();
					for (i = 0; i < M; i++) {
						//Bk[i] = Bk[i] * complex(cos(w * (i * da + (float)j / Nw)), -sin(w * (i * da + (float)j / Nw)));
						SCoCmplx[k * a2 + i] = SCoCmplx[k * a2 + i] + Bk[i] * sfj;
						sfj = sfj * sfi;
					}
				}

			}
			else {//aliased (for coherence)
				if (opt > 1) {
					for (i = 0; i < KN; i++) {
						Ak[i] = yws[k + i * hnfft] * yws[(j - k) + i * hnfft];
					}
#ifdef USEFFTW
					fftw_execute(fp2);
#else
					ffthelper(Ak, Bk, KNP2);
#endif // USEFFTW	
					sfj = sfj1.conj();
					for (i = 0; i < M; i++) {
						//Bk[i] = Bk[i] * complex(cos(w * (i * da + (float)j / Nw)), -sin(w * (i * da + (float)j / Nw)));
						Bk[i] = Bk[i] * sfj;
						sfj = sfj * sfi;
						SCoCmplx[k * a2 + i] = SCoCmplx[k * a2 + i] + Bk[i];
					}
				}
			}
		}
	}
	for (i = 0; i < M; i++)
		SRw[i] = Rw[i]; //Rw(alpha - 0)
	for (j = 1; j <= P; j++) {
		for (i = 0; i < M; i++) {
			ir = i - (int)(j / (da * Nw)); //calculate shifted Rw(alpha - p * deltaF)
			ir2 = i + (int)(j / (da * Nw)); //calculate shifted Rw(alpha - p * deltaF) when p < 0
			if (ir < 0)  ir = -ir; //Rw(alpha) is even function
			SRw[i] = Rw[ir] + Rw[ir2] + SRw[i];
		}
	}
	for (i = 0; i < hnfft; i++)
	{
		for (j = 0; j < M; j++) {
			if (j == 0) {
				SCo[i + j * hnfft] = 0.0f; //nullify alpha=0 component to protect SCho map
			}
			else {
				SCo[i + j * hnfft] = RW0 * SCoCmplx[i * a2 + j].amp() / KMU / SRw[j];
			}
		}

	}
	fftw_destroy_plan(fp);
	fftw_free(a);
	fftw_free(b);
	fftw_destroy_plan(fp2);
	fftw_free(Ak);
	fftw_free(Bk);
	delete[] windat;
	delete[] yws; delete[] Rw; delete[] SRw; delete[] CPY2;
	delete[] SCoCmplx;
	return 0;
}


int FAM(float* x, float &da, float amax, int L, int Nw, int opt, float* SCo)
{
	/*
	calculate spectral correlation/coherance using using FFT Accumulation Method, Roberts et al.

	da: delta alpha relative to Fs
	x[]: input data
	amax: maximum cyclic frequency (normalized by Fs)
	L: input data length
	to get best speed and accuracy, make sure L gives power of two for KN , where KN = (L - Nw + R) / R;=> L = R*KN + Nw - R
	if KN is not power-of-two, the function will adjust it by zero padding and update da
	Nw: window length        '
	Nw: window length        '
	opt:0 :SC no negative frequencies correlation, 1 : Scoh, 2: SC with negative frequencies, 3: Scoh with negative freq
	output:
	updated da,
	SCo[] output data: size Nw*a2
	*/
	int nfft, hnfft, KN, KNP2, i, k, j, R, ap, i1, j1, q, qr, R0, k1, aindex, findex, mf;
	float alpha0;

	fftw_complex *in, *out, *ink2, *outk2;
	complex *a, *b, *Ak, *Bk;
	fftw_plan fp, fp2;
	in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * Nw); //new complex[Nw];
	out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * Nw); //new complex[Nw];
	a = (complex*)&in[0];
	b = (complex*)&out[0];

	fp = fftw_plan_dft_1d(Nw, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	auto windat = new double[Nw];


	nfft = Nw;
	hnfft = nfft / 2;
	//Prep WinData;
	for (i = 0; i < Nw; i++)
		windat[i] = (0.54 - 0.46 * cos(2 * i * pi / Nw)); //hamming window
	R = Nw / 4; //initial
	R0 = R;
	KN = (int)((L - Nw + R) / R);
	//KNP2 = (int)pow(2.0, ceil(log2(KN)));//next power of two	
	KNP2 = nextPow2(KN);
	qr = KNP2 *  R / (2 * Nw);

	ink2 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * KNP2); //new complex[KNP2];
	outk2 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * KNP2); //new complex[KNP2];
	Ak = (complex*)&ink2[0]; //create complex alias
	Bk = (complex*)&outk2[0];
	fp2 = fftw_plan_dft_1d(KNP2, ink2, outk2, FFTW_FORWARD, FFTW_ESTIMATE);

	da = 1.0F / (KNP2 * R); //updated normalized da
	int a2 = (int)(amax / da);
	//std::cout << "KN, KNP2: " << KN << ", " << KNP2 << "\n";
	//std::cout << "updated da: " << da << "\n";
	auto yws = new complex[nfft * KN];
	auto SCoCmplx = new complex[(nfft + 1) * a2];
	auto CPY2 = new double[nfft];
	double KMU = 0.0, w;

	for (i = 0; i < (nfft + 1) * a2; i++)
		SCoCmplx[i] = complex(0.0f, 0.0f);
	for (i = 0; i < Nw; i++)
		KMU = KMU + windat[i] * windat[i];
	KMU = KN * KMU / 2.0f;


	for (i = 0; i < hnfft; i++)
	{
		CPY2[i] = 0.0;
	}
	w = 2.0f * pi * R / Nw;
	for (k = 0; k < KN; k++) // calculate initial spectral data
	{

		for (i = 0; i < nfft; i++)
			if ((i + k * R) < L) {
				a[i] = complex(x[i + k * R] * windat[i], 0.0);
			}
			else {
				a[i] = complex(0.0, 0.0);
			}

#ifdef USEFFTW
			fftw_execute(fp);
#else
			ffthelper(a, b, Nw);
#endif // USEFFTW	

			for (i = 0; i < nfft; i++)
			{
				yws[i + k * nfft] = b[i] * complex(cos(w * i * k), -sin(w * i * k)); //phase compensation;
				CPY2[i] = CPY2[i] + b[i].real * b[i].real + b[i].imag * b[i].imag;
			}
	}
	if ((opt == 1) || (opt == 3)) {
		for (k = 0; k < KN; k++) {
			for (i = 0; i < nfft; i++)
				yws[i + k * nfft] = yws[i + k * nfft] / (sqrt(CPY2[i] / KMU));//to calculate coherence
		}
	}
	for (k = 0; k < KNP2; k++) {
		Ak[k] = complex(0.0f, 0.0f);
	}

	mf = nfft;
	if (opt < 2) mf = hnfft;
	for (j = 0; j < mf; j++) {
		for (i = 0; i < mf; i++) {

			for (k = 0; k < KN; k++) {
				Ak[k] = yws[j + k * nfft] * yws[i + k * nfft].conj();
			}
#ifdef USEFFTW
			fftw_execute(fp2);
#else
			ffthelper(Ak, Bk, KNP2);
#endif // USEFFTW	
			if (i >= hnfft) {
				i1 = i - nfft;
			}
			else {
				i1 = i;
			}
			if (j >= hnfft) {
				j1 = j - nfft;
			}
			else {
				j1 = j;
			}
			//f0 = (i1 + j1) / 2;
			findex = i1 + j1; //FAM provides double resolution in spectral frequency, so index is multipled by 2
			alpha0 = (float)(j1 - i1) / (Nw * da);
			//
			ap = (int)floor(alpha0);
			//DbgMsg(" ap: " << ap);
			for (q = -qr; q < qr; q++) {
				aindex = ap + q;
				if ((aindex >= 0) && (findex >= 0) && (aindex < a2)) {
					if (q < 0) {
						k1 = KNP2 + q;
					}
					else {
						k1 = q;
					}

					SCoCmplx[findex * a2 + aindex] = Bk[k1] + SCoCmplx[findex * a2 + aindex];
				}

			}

		}
	}

	for (i = 0; i < nfft; i++)
	{
		for (j = 0; j < a2; j++) {
			if (j == 0) {
				SCo[i + j * hnfft] = 0.0f; //nullify alpha=0 component to protect SCho map
			}
			else {
				SCo[i + j * nfft] = (float)(SCoCmplx[i * a2 + j].amp() / KMU);
			}
		}

	}

	fftw_destroy_plan(fp);
	fftw_free(in);
	fftw_free(out);
	fftw_destroy_plan(fp2);
	fftw_free(ink2);
	fftw_free(outk2);

	delete[] windat;
	delete[] yws; delete[] CPY2;
	delete[] SCoCmplx;
	return 0;
}
//=============================================================================================
int AFAM(float* x, float &da, float amax, int L, int Nw, int opt, float* SCo)
{
	/*
	calculate spectral correlation/coherance using using Aggregated FFT Accumulation Method, J et al.

	da: delta alpha relative to Fs
	x(): input data
	amax: maximum cyclic frequency (normalized by Fs)
	L: input data length
	to get best speed, make sure L gives power of two for KN , where KN = (L - Nw + R) / R;=> L = R*KN + Nw - R
	If KN is not power-of-two, the function will adjust it and update da
	Nw: window length        '
	opt:0 :SC no negative frequencies correlation, 1 : Scoh, 2: SC with negative frequencies, 3: Scoh with negative freq
	output:
	updated da
	SCo() output data: size Nw*a2
	*/
	int nfft, hnfft, KN, KNP2, i, k, j, R, ap, i1, j1, q, qr, R0, k1, aindex, findex, mf, M;
	float alpha0;

	fftw_complex *in, *out, *ink2, *outk2;
	complex *a, *b, *Ak, *Bk;
	fftw_plan fp, fp2;
	in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * Nw); //new complex[Nw];
	out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * Nw); //new complex[Nw];
	a = (complex*)&in[0];
	b = (complex*)&out[0];

	fp = fftw_plan_dft_1d(Nw, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	auto windat = new double[Nw];


	nfft = Nw;
	hnfft = nfft / 2;
	//Prep WinData;
	for (i = 0; i < Nw; i++)
		windat[i] = (0.5 - 0.5 * cos(2 * i * pi / Nw)); //hamming window
	R = Nw / 4; //initial

	KN = (int)((L - Nw + R) / R);
	//KNP2 = (int)pow(2, ceil(log2(KN)));//next power of two	
	KNP2 = nextPow2(KN);
	DbgMsg("AFAM R: " << R << ", KN: " << KN << ", KNP2: " << KNP2 << "\n";);

	qr = KNP2 / 2; // 'in original FAM is is limited by KNP2 * R / (2 * Nw)	

	ink2 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * KNP2); //new complex[KNP2];
	outk2 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * KNP2); //new complex[KNP2];
	Ak = (complex*)&ink2[0]; //create complex alias
	Bk = (complex*)&outk2[0];
	fp2 = fftw_plan_dft_1d(KNP2, ink2, outk2, FFTW_FORWARD, FFTW_ESTIMATE);

	da = 1.0F / (KNP2 * R); //updated normalized da
							//std::cout << "KN, KNP2: " << KN << ", " << KNP2 << "\n";
							//std::cout << "updated da: " << da << "\n";
	int a2 = (int)(amax / da);
	auto yws = new complex[nfft * KN];
	auto SCoCmplx = new complex[(nfft + 1) * a2];
	auto SRw = new double[(nfft + 1) * KNP2];
	auto Rw = new double[qr + 1];
	auto CPY2 = new double[nfft];
	double w, KMU;
	complex sfi, sfa;
	for (i = 0; i < (nfft + 1) * a2; i++)
		SCoCmplx[i] = complex(0.0f, 0.0f);
	KMU = 0.0f;
	for (i = 0; i < Nw; i++)
		KMU = KMU + windat[i] * windat[i];
	KMU = (KN * KMU / 2.0f);
	M = Nw / 2;
	for (q = 0; q <= qr; q++) {
		Rw[q] = 0;
		for (i = 0; i < Nw; i++) {
			Rw[q] = Rw[q] + windat[i] * windat[i] * cos(2 * pi * (i - M) * q * da); //this is even function
		}
	}
	for (i = 0; i < nfft; i++) {
		for (j = 0; j < KNP2; j++) {
			SRw[i * KNP2 + j] = 0.0;
		}
	}

	for (i = 0; i < hnfft; i++)
	{
		CPY2[i] = 0.0;
	}
	w = 2.0f * pi * R / Nw;
	for (k = 0; k < KN; k++) // calculate initial spectral data
	{

		for (i = 0; i < nfft; i++)
			if ((i + k * R) < L) {
				a[i] = complex(x[i + k * R] * windat[i], 0.0);
			}
			else {
				a[i] = complex(0.0, 0.0);
			}

#ifdef USEFFTW
			fftw_execute(fp);
#else
			ffthelper(a, b, Nw);
#endif // USEFFTW	

			for (i = 0; i < nfft; i++)
			{
				yws[i + k * nfft] = b[i] * complex(cos(w * i * k), -sin(w * i * k)); //phase compensation;
				CPY2[i] = CPY2[i] + b[i].real * b[i].real + b[i].imag * b[i].imag;
			}
	}
	if ((opt == 1) || (opt == 3)) {
		for (k = 0; k < KN; k++) {
			for (i = 0; i < nfft; i++)
				yws[i + k * nfft] = yws[i + k * nfft] / (sqrt(CPY2[i] / KMU));//to calculate coherence
		}
	}
	for (k = 0; k < KNP2; k++) {
		Ak[k] = complex(0.0f, 0.0f);
	}

	mf = nfft;
	w = 2 * pi * hnfft * da;
	if (opt < 2) { mf = hnfft; }
	//DbgMsg("AFAM R: " << R << ", mf: " << mf << ", hnfft: " << hnfft << ", opt: " << opt << "\n";);
	sfi = complex(cos(w), -sin(w));
	for (j = 0; j < mf; j++) {// fj
		for (i = 0; i < mf; i++) { //fm
			if (i >= hnfft) {
				i1 = i - nfft;
			}
			else {
				i1 = i;
			}
			if (j >= hnfft) {
				j1 = j - nfft;
			}
			else {
				j1 = j;
			}
			//f0 = (i1 + j1) / 2;
			findex = j;
			alpha0 = (float)(j1 - i1) / (Nw * da);
			//
			ap = (int)lround(alpha0);
			if ((ap + qr - 1) < 0 || (ap - qr) > a2) //check boundaries
			{
				continue;
			}
			for (k = 0; k < KN; k++) {
				Ak[k] = yws[j + k * nfft] * yws[i + k * nfft].conj();
			}
#ifdef USEFFTW
			fftw_execute(fp2);
#else
			ffthelper(Ak, Bk, KNP2);
#endif // USEFFTW	
			//DbgMsg(" ap: " << ap);
			sfa = complex(cos(w * qr), sin(w * qr));
			for (q = -qr; q < qr; q++) {
				aindex = ap + q;
				if ((aindex >= 0) && (aindex < a2)) {
					if (q < 0) {
						k1 = KNP2 + q;
					}
					else {
						k1 = q;
					}
					if (aindex < KNP2) { //SRw is periodic with period KNP2
						SRw[findex * KNP2 + aindex] += Rw[abs(q)];
					}

					SCoCmplx[findex * a2 + aindex] = SCoCmplx[findex * a2 + aindex] + Bk[k1] * sfa;
				}
				sfa = sfa * sfi; //update the phase
			}

		}
	}

	for (i = 0; i < hnfft; i++) {
		for (j = 0; j < a2; j++) {
			if (j == 0) {
				SCo[i + j * hnfft] = 0.0f; //nullify alpha=0 component to protect SCho map
			}
			else {
				j1 = j % KNP2; //SRw is periodic with period KNP2
				w = abs(SRw[i * KNP2 + j1])*KN / KNP2; //KN/KNP2 account for zero-padding
				if (w > 1e-15) {
					SCo[i + j * hnfft] = (float)(SCoCmplx[i * a2 + j].amp() / KN / w);
				}
				else {
					SCo[i + j * hnfft] = 0.0f;
				}
			}
		}

	}

	fftw_destroy_plan(fp);
	fftw_free(in);
	fftw_free(out);
	fftw_destroy_plan(fp2);
	fftw_free(ink2);
	fftw_free(outk2);

	delete[] windat;
	delete[] yws; delete[] CPY2;
	delete[] SCoCmplx;
	delete[] SRw;
	delete[] Rw;
	return 0;
}

// ---------------------------------------------------------------------
// Choose the optimal OLS block length M for FastSCP (Eq. 23):
//   M is a power of 2, M >= Nb, M <= 2^ceil(log2(inputLen))
//   minimizing C = (2*Ns+1)*M*log2(M) + Ns*M
// ---------------------------------------------------------------------
static int chooseOptimalM(int inputLen, int Nb) {
	int Mmin = nextPow2(Nb);
	int Mmax = nextPow2(inputLen);
	if (Mmax < Mmin) Mmax = Mmin;
	int bestM = Mmin;
	double bestCost = -1.0;
	for (int M = Mmin; M <= Mmax; M <<= 1) {
		int R = M - Nb + 1;
		if (R <= 0) continue;
		int Ns = (inputLen + R - 1) / R; // ceil(inputLen / R)
		double cost = (2.0 * Ns + 1.0) * M * std::log2((double)M) + (double)Ns * M;
		if (bestCost < 0.0 || cost < bestCost) {
			bestCost = cost;
			bestM = M;
		}
	}
	return bestM;
}

// Gamma must be caller-allocated with size Na*K (same layout as Sc).
// ------------------------------------------------------------------
void computeSpectralCoherence(const complex* Sc, int Na, int K, double dAlpha, double dF, complex* Gamma)
{

	auto interpS0 = [&](double f) -> complex {
		//if (f < 0.0) return complex(0.0, 0.0);
		if (f < 0.0) f = -f;
		if (f < dAlpha) return complex(0.0, 0.0);
		double pos = f / dF;
		int k0 = (int)std::floor(pos);
		if (k0 >= K - 1) return (k0 < K) ? Sc[k0] : complex(0.0, 0.0);
		double frac = pos - k0;
		return Sc[k0] * (1.0 - frac) + Sc[k0 + 1] * frac;

		//int k0 = (int)std::llround(pos);
		//if (k0 > K - 1) return complex(0.0, 0.0);
		//return Sc[k0];
	};

	for (int k = 0; k < K; ++k) {
		double f = k * dF;
		complex s0f = Sc[k];
		for (int a = 0; a < Na; ++a) {
			double alpha = a * dAlpha;
			complex s0fma = interpS0(f - alpha);
			complex denom = csqrt(s0f * s0fma);
			complex num = Sc[k + a*K];
			Gamma[k + a*K] = (denom.amp() > 0.0) ? num / denom : complex(0.0, 0.0);
		}
	}

}
// ---------------------------------------------------------------------

// ---------------------------------------------------------------------
int FastSCP(const float* x, float &da, int Lsig, int Nw, int Nb, float amax, int opt, float* SC) {
	//opt=0: SC, 1: SCoh
	// FastSCP, by:
	// Y. Chen, J. Wang, L. Qiu, G. Liang, Y. Li, 
	// Fast computation of the spectral correlation via frequency-averaging, Mech. Syst. Signal Process. 223 (2025) 111851.
	// Inputs
	//   x[]      : real input signal, length Lsig
	//   Lsig     : signal length L
	//   Nw       : used only to set spectral resolution, the number of spectral frequencies K = Nw/2
	//   Nb       : smoothing point count (frequency-smoothing window length)
	//   amax     : max. alpha (normalized by Fs)
	//
	// Output
	//   SC       : float real output buffer of size K*Na (K = Nw/2), pre-allocated by the caller.
	//               data interleaved as:
	//                  SC[k + K * aindex] = |Ŝ_x(alpha_aindex, f_k)|
	//              i.e. for fixed spectral frequency f_k, the Na cyclic-
	//              frequency bins are contiguous.
	//              completly written from scratch by Jaf

	const int K = Nw / 2; // number of (positive) spectral frequencies

	const int L = Lsig;
	fftw_complex *XLf;
	complex* XL;
	fftw_plan fp;
	XLf = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * L); //fftw data type	
	XL = (complex*)&XLf[0];
	fp = fftw_plan_dft_1d(L, XLf, XLf, FFTW_FORWARD, FFTW_ESTIMATE);
	for (int n = 0; n < L; n++) XL[n] = complex(x[n], 0.0);
	//fft(XL, L, 1);
	fftw_execute(fp); //we must use FFTW since L could be not power-of-two to keep da and prevent smearing due to zero-padding
	for (int n = L / 2 + 1; n < L; n++) XL[n] = complex(0.0, 0.0); //negelect negative spectrum
	da = 1.0f / (float)L; //adjusted dAlpha


						  // ---- Choose OLS parameters (Eq. 23) ----
	int Na = (int)(amax / da);
	const int inputLen = Na + Nb - 1;      // length of U(fk)
	int M = chooseOptimalM(inputLen, Nb);
	//M = 512;
	const int R = M - Nb + 1;             // block shift
	const int Ns = (inputLen + R - 1) / R; // number of segments (ceil)
										   //std::cout << "M: " << M << ", R: " << R << ", Ns: " << Ns << " \n";

										   // ---- Scratch buffers reused across the k loop ----
	complex* U = new complex[inputLen]; // U(fk), Eq. (22) top row
	complex* V = new complex[Nb];       // V(fk), Eq. (22) bottom row	
	complex* S2x = new complex[K*Na];
	fftw_complex *Hf, *Ublockf, *Yf;
	complex *H, *Ublock, *Y;
	fftw_plan fp2, fp3, fi;
	Hf = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * M);        // FFT{V(fk), M}
	Ublockf = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * M);        // i-th zero-padded block of U(fk)
	Yf = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * M);        // IFFT{ FFT{Ublock,M} .* H, M }
	H = (complex*)&Hf[0];
	Ublock = (complex*)&Ublockf[0];
	Y = (complex*)&Yf[0];
	fp2 = fftw_plan_dft_1d(M, Hf, Hf, FFTW_FORWARD, FFTW_ESTIMATE); //estimate H
	fp3 = fftw_plan_dft_1d(M, Ublockf, Ublockf, FFTW_FORWARD, FFTW_ESTIMATE); //estimate fft(U)
	fi = fftw_plan_dft_1d(M, Yf, Yf, FFTW_BACKWARD, FFTW_ESTIMATE); //inverse transform ifft(Y)
	for (int k = 0; k < K; k++) {
		// spectral-frequency bin f_k = k * Fs/Nw, expressed as a bin
		// index of the L-point spectrum: fk_bin = k * L / Nw
		const int fk_bin = (int)std::llround((double)k * L / (double)Nw);

		// ---- Build U(fk), Eq. (22): U[i] = X_L(fk - i*Δα) ----
		for (int i = 0; i < inputLen; i++) {
			int m = (fk_bin - i) % L; // Circular fetch of XL at an arbitrary (possibly negative) bin index.
			if (m < 0) m += L;
			U[i] = XL[m];
		}

		// ---- Build V(fk), Eq. (22): V[j] = conj(X_L(fk - (Nb-1-j)*Δα)) ----
		for (int j = 0; j < Nb; j++) {
			int m = (fk_bin - (Nb - 1 - j)) % L;
			if (m < 0) m += L;
			V[j] = XL[m].conj();
		}

		// ---- Step 2: H = FFT{V(fk), M} ----		
		for (int i = 0; i < Nb; i++) H[i] = V[i];
		for (int i = Nb; i < M; i++) H[i] = complex(0.0, 0.0);

#ifdef USEFFTW // USEFFTW	
		fftw_execute(fp2);
#else
		fft(H, M, 1);
#endif 

		// ---- Steps 1,3,4: partition U(fk) into Ns blocks and OLS-convolve ----
		for (int seg = 0; seg < Ns; seg++) {
			const int start = seg * R;

			// Block i of U(fk): U[start .. start+M-1], zero-padded past the end.
			for (int j = 0; j < M; j++) {
				const int srcIdx = start + j;
				Ublock[j] = (srcIdx < inputLen) ? U[srcIdx] : complex(0.0, 0.0);
			}

			// Y = IFFT{ FFT{Ui(fk), M} .* H, M }

#ifdef USEFFTW // USEFFTW	
			fftw_execute(fp3);
#else
			fft(Ublock, M, 1);
#endif 
			for (int j = 0; j < M; j++) Y[j] = Ublock[j] * H[j];
#ifdef USEFFTW // USEFFTW	
			fftw_execute(fi);
#else
			fft(Y, M, -1);//inverse transform
#endif 			
			for (int j = 0; j < M; ++j) Y[j] = Y[j] / (double)M;

			// Discard first Nb-1 points; keep Y[Nb-1 .. M-1] (length R)
			// and map them to cyclic-frequency indices seg*R .. seg*R+R-1.
			for (int j = 0; j < R; j++) {
				const int alphaIdx = start + j;
				if (alphaIdx >= Na) break; // last block may run past Na
				S2x[k + K * alphaIdx] = Y[Nb - 1 + j];
			}
		}
	}
	if ((opt == 0) || (opt == 2)) { //SC
		for (int k = 0; k < K; k++) {
			for (int a = 0; a < Na; a++) {
				SC[k + K * a] = S2x[k + K * a].amp();
			}
		}
		// --- Step 8: magnitude adjustment,  Sc <- Sc / (L * Delta_f) ---> proper adjustment without Fs, use L * Nb
		const double norm = 1.0 / ((double)L * (double)Nb);
		for (int i = 0; i < Na * K; ++i) SC[i] *= norm;
	}
	else { //Spectral coherence
		complex* SCoh = new complex[K*Na];
		computeSpectralCoherence(S2x, Na, K, (double)da, (double)(1.0 / Nw), SCoh);
		for (int k = 0; k < K; k++) {
			for (int a = 0; a < Na; a++) {
				if (a == 0) {
					SC[k + K * a] = 0.0f; //nullify alpha=0 SCho data to protect the map
				}
				else {
					SC[k + K * a] = SCoh[k + K * a].amp();
				}
			}
		}
		//no need for magnitude adjustment for Coherence
		delete[] SCoh;
	}

	delete[] U;
	delete[] V;
	delete[] S2x;
	fftw_destroy_plan(fp);
	fftw_free(XLf);
	fftw_destroy_plan(fp2);
	fftw_destroy_plan(fp3);
	fftw_destroy_plan(fi);
	fftw_free(Yf);
	fftw_free(Ublockf);
	fftw_free(Hf);
	return 0;
}
int FDirSC(float* x, float &da, float amax, int L, int Nw, int opt, float* SCo)
{
	/*
	calculate spectral correlation/coherance using  Faster SC, Borghesani and Antoni 2018

	da: delta alpha relative to Fs
	x() input data
	amax: maximum cyclic frequency (normalized by Fs)
	L: input data length
	to get best speed and accuracy, make sure L gives power of two for KN , where KN = (L - Nw + R) / R;=> L = R*KN + Nw - R
	if KN is not power-of-two, the function will adjust it by zero padding and update da
	Nw: window length        '
	opt:0 :SC no negative frequencies correlation, 1 : Scoh, 2: SC with negative frequencies, 3: Scoh with negative freq
	output:
	adjusted da
	SCo() output data
	*/
	int nfft, hnfft, KN, KNP2, i, k, j, R, P, M, ri, KNR;
	float Nb;
	R = (int)floor(1.0f / (2.0f * amax));
	//std::cout << "FDirSC, R: " << R << "\n";
	if (R > (Nw / 4)) R = Nw / 4;
	DbgMsg("FDirSC, R: " << R << "\n");

	nfft = Nw;
	hnfft = nfft / 2;
	auto windat = new double[Nw];
	fftw_complex *in, *out, *ink, *outk, *inr;
	fftw_plan fp, fp2, fp3;

	complex* a = (complex*)fftw_malloc(sizeof(complex) * Nw); //new complex[Nw];
	complex* b = (complex*)fftw_malloc(sizeof(complex) * Nw); //new complex[Nw];
	in = (fftw_complex*)&a[0]; //(fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
	out = (fftw_complex*)&b[0];
	fp = fftw_plan_dft_1d(Nw, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	//Prep WinData;
	for (i = 0; i < Nw; i++)
		windat[i] = (0.5 - 0.5 * cos(2 * i * pi / Nw));

	//std::cout << "P: " << P << "\n";
	KN = (int)((L - Nw + R) / R);
	//KNP2 = (int) pow(2.0, ceil(log2(KN)));//next power of two
	KNP2 = nextPow2(KN);
	//std::cout << "KN, KNP2: " << KN << ", " << KNP2 << "\n";
	R = (int)ceil((float)(L - Nw) / (float)(KNP2 - 1)); //updated R
	if (R < 1) R = 1;
	KN = (int)((L - Nw + R) / R); //KN will be less than KNP2 but larger than previous KN
	P = (int)floor((float)Nw / (2.0f * R));
	DbgMsg("FDirSC KN, KNP2: " << KN << ", " << KNP2 << "\n");
	da = 1.0F / ((float)KNP2 * (float)R); //updated normalized da
										  //std::cout << "updated da: " << da << "\n";
	int a2 = (int)(amax / da);
	Nb = (1.0f / da) / Nw;
	complex* Ak = (complex*)fftw_malloc(sizeof(complex) * KNP2); //new complex[Nw];
	complex* Bk = (complex*)fftw_malloc(sizeof(complex) * KNP2); //new complex[Nw];
	ink = (fftw_complex*)&Ak[0]; //(fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
	outk = (fftw_complex*)&Bk[0];
	fp2 = fftw_plan_dft_1d(KNP2, ink, outk, FFTW_FORWARD, FFTW_ESTIMATE);
	//std::cout << "KN: " << KN << "\n";
	auto yws = new complex[hnfft * KN];
	auto ywsd = new complex[hnfft * KN];
	auto SCoCmplx = new complex[hnfft * a2];
	auto CPY2 = new double[hnfft];
	auto Dp = new complex[Nw];
	double KMU = 0.0f, w;
	complex s1;
	KNR = KNP2*R;
	complex* Ar = (complex*)fftw_malloc(sizeof(complex) * KNR); //new complex[Nw];	
	inr = (fftw_complex*)&Ar[0]; //(fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);

	fp3 = fftw_plan_dft_1d(KNR, inr, inr, FFTW_FORWARD, FFTW_ESTIMATE);//in-place transform
	for (i = 0; i < Nw; i++)
		KMU = KMU + windat[i] * windat[i];

	M = Nw / 2 - 1;
	w = 2.0f * pi / Nw;
	for (i = 0; i < Nw; i++) //generate Dirichlet kernel
	{
		s1 = complex(0, 0);
		for (j = -P; j <= P; j++)
			s1 = s1 + complex(cos(w * j * (i - hnfft)), sin(w * j * (i - hnfft)));
		Dp[i] = s1;
	}
	for (i = 0; i < KNR; i++) //generate calibration data
	{
		if (i < Nw) {
			Ar[i] = windat[i] * windat[i] * Dp[i];
		}
		else {
			Ar[i] = complex(0.0f, 0.0f);
		}
	}
	fftw_execute(fp3); //calibration data now in Ar[]

	for (i = 0; i < hnfft; i++)
	{
		CPY2[i] = 0.0;
	}
	for (k = 0; k < KN; k++) // calculate initial spectral data
	{

		for (i = 0; i < nfft; i++)
			a[i] = complex(x[i + k * R] * windat[i], 0.0);
#ifdef USEFFTW
		fftw_execute(fp);
#else
		ffthelper(a, b, Nw);
#endif // USEFFTW	
		for (i = 0; i < hnfft; i++)
		{
			yws[i + k * hnfft] = b[i];
			CPY2[i] = CPY2[i] + b[i].real * b[i].real + b[i].imag * b[i].imag;
		}
		for (i = 0; i < nfft; i++)
			a[i] = a[i] * Dp[i];
		//fftw_execute(fp);
		ffthelper(a, b, Nw);
		for (i = 0; i < hnfft; i++)
		{
			ywsd[i + k * hnfft] = b[i];
		}
	}
	/*if ((opt == 1) || (opt == 3)) {  //not working as expected, use standard coherence function
	for (k = 0; k < KN; k++) {
	for (i = 0; i < hnfft; i++)
	yws[i + k * hnfft] = yws[i + k * hnfft] / (sqrt(CPY2[i] / KMU));
	}
	}*/

	for (k = 0; k < KNP2; k++) {
		Ak[k] = complex(0.0f, 0.0f);
	}
	M = KN / 2;
	if (a2 < M) M = a2;
	for (k = 0; k < hnfft; k++) {
		for (i = 0; i < KN; i++) {
			Ak[i] = yws[k + i * hnfft] * ywsd[k + i * hnfft].conj();
		}
#ifdef USEFFTW
		fftw_execute(fp2);
#else
		ffthelper(Ak, Bk, KNP2);
#endif // USEFFTW	
		for (i = 0; i < M; i++)
			SCoCmplx[k + i*hnfft] = SCoCmplx[k + i*hnfft] + Bk[i] / Ar[i];

	}

	if ((opt == 0) || (opt == 2))
	{

		for (i = 0; i < hnfft; i++)
		{
			for (j = 0; j < a2; j++)
				SCo[i + j * hnfft] = (float)(SCoCmplx[i + j*hnfft].amp() / KN);
		}
	}
	else
	{
		for (i = 0; i < hnfft; i++)
		{
			for (j = 0; j < a2; j++)
			{
				ri = (int)round(i - (j / Nb));
				if (ri < 0) ri = -ri;
				if (j == 0) {
					SCo[i + j * hnfft] = 0.0;
				}
				else {
					SCo[i + j * hnfft] = (float)(KMU * SCoCmplx[i + j*hnfft].amp() / sqrt(CPY2[i] * CPY2[ri])); // pick the nearest bin
				}
			}
		}
		//complex* SCoh = new complex[hnfft*a2];
		//computeSpectralCoherence(SCoCmplx, a2, hnfft, (double)da, (double)(1.0/Nw), SCoh);
		//for (i = 0; i < hnfft; i++)
		//{
		//	for (j = 0; j < a2; j++)
		//	{					
		//		SCo[i + j * hnfft] = SCoh[i + j*hnfft].amp(); // pick the nearest bin
		//	}
		//}
		////no need for magnitude adjustment for Coherence
		//delete[] SCoh;

	}
	fftw_destroy_plan(fp);
	fftw_free(a);
	fftw_free(b);
	fftw_destroy_plan(fp2);
	fftw_free(Ak);
	fftw_free(Bk);
	fftw_destroy_plan(fp3);
	fftw_free(Ar);

	delete[] windat;
	delete[] yws; delete[] ywsd; delete[] Dp; delete[] CPY2;
	delete[] SCoCmplx;
	return 0;
}


void printStat(double T[][22], char* label, int alpharange, int runrange)
{
	int i, j;//i is the alpha index, j is the run index
	double minT, maxT,avg;
	std::cout << "============== Method: "<<label<<" ===================\n";
	std::printf("%s min Times:\n",label);
	for (i = 0; i < alpharange; ++i) {
		minT = 1e200;
		for (j = 0; j < runrange; ++j) {
			if (T[i][j] < minT) minT = T[i][j];			
		}
		std::printf(" %.4f,",  minT);
	}
	std::printf("\n %s max Times:\n", label);
	for (i = 0; i < alpharange; ++i) {
		maxT = 0;
		for (j = 0; j < runrange; ++j) {
			if (T[i][j] > maxT) maxT = T[i][j];
		}
		std::printf(" %.4f,", maxT);
	}
	std::printf("\n %s Average Times:\n", label);
	for (i = 0; i < alpharange; ++i) {
		avg = 0;
		for (j = 0; j < runrange; ++j) {
			avg +=T[i][j] ;
		}
		std::printf(" %.4f,", avg/(double)runrange);
	}
	std::printf("\n");
	
}

int main()
{
	float da, wdf,SamRate,amax,Fs;
	int a1, a2,L, Nw, Noverlap, i,ia,Rfam,Lfam,Nb,nf,ainc,runid;
	float fc, fm;
	double f2, fr, ffault, zeta, Tp, b, ph2, pph2, dt, fring;
	int npos;
	double T1[5][22], T2[5][22], T3[5][22], T4[5][22], T5[5][22], T6[5][22];

	int checkid = 2; //checking type: 0: nothing, 1: check estimation values,
					 // 2: check run time for alpha 0.0325 to 0.5 L=2^17, 3: check run time alpha = 0.05, L =2^18
	L = 1 << 17; //131072;
	if (checkid == 3) {
		L = 1 << 18; //262144;
	}
	
	SamRate = 32768;
	Fs = SamRate;
	Nw = 512;
	Nb = 256; //FastSCP smoothing length
	srand(time(NULL));   // Initialization, should only be called once.
	
	
	zeta = 0.16;
	dt = 1 / Fs;
	fring = (Fs/Nw)*25; //ringing freq
	fr = 20; //rotational speed of shaft
	ffault = 5; //bearing fault frequency in X
	f2 = ffault * fr; //fault rate in Hz
	
	int* iPos = new int[L];
	
	
	nf = Nw / 2;
	Rfam = Nw / 4;
	Lfam = L + Nw - Rfam; //used to get ready power-of-two KN since KN = (int)((L - Nw + R) / R);
	Noverlap = (int)(Nw * 0.75);//or R = 0.333Nw for ACP and Fast ACP
	
	wdf = SamRate / Nw;//delta freq
	da = 1.0f / (float)L; //relative to Sample rate
	a1 = 0;
	amax = 0.05;
	a2 = (int)ceil(amax / da);
	auto SCo = new float[Nw  * L];//(L/2) is maximum alpha count, but we perserve double memory for safe operation
	auto EES = new float[L];//Envelope spectrum
	auto x = new float[Lfam+32]();//putting () will initilize to zero all elements

	fc = (Fs / Nw) * 25; //carrier freq (spectral) for simple case
	fm = 10;              //modulating frequency (alpha) for simple case
	Tp = (1 / f2) / dt;
	npos = 0;
	b = 0;
	ph2 = 0;
	for (int i = 0; i < L; ++i) {
		//double t = i / Fs;
		//double am = 1;
		//if (std::cos(2 * pi * fm * t) < 0) am = -1;
		//x[n] = (1.0 + 0.8 * std::cos(2 * pi * fm * t)+ 0.4 * std::cos(4 * pi * fm * t)) * std::cos(2 * pi * fc * t);
		//x[n] = (1.0 + 0.5 * am) * std::cos(2 * pi * fc * t);

		//gen simulated bearing signal
		pph2 = ph2;
		ph2 = ph2 + 2 * pi * f2 * dt;
		if (b == 0) {
			npos += 1; // number of pulses
			iPos[npos] = i;
		}
		if (npos < 6) {
			for (int j = 1; j <= npos; ++j) {
				x[i] += std::cos(2 * pi * fring * (i - iPos[j]) * dt) * std::exp(-(i - iPos[j]) * dt * 2 * pi * fring * zeta); // smooth fault
			}
		}
		else {
			for (int j = npos - 5; j <= npos; ++j) {
				x[i] += std::cos(2 * pi * fring * (i - iPos[j]) * dt) * std::exp(-(i - iPos[j]) * dt * 2 * pi * fring * zeta); // smooth fault
			}
		}

		b += 1;
		if (std::sin(ph2) >= 0 && std::sin(pph2) < 0) {
			if (b > Tp / 2) { // to suppress early switching	
				b = 0;
			}
		}
	}
	if (checkid == 1) { //check estimaed values for bearing signal 
		int kIdx = (int)std::llround(fring / (Fs / Nw));
		SCoh(x, da, amax, L, Noverlap, Nw, 0, SCo); //standard ACP		
		ainc = (int)std::llround(f2 / (da*Fs));
		std::cout << "Spectral Freq index=" << kIdx << ", a inc index:" << ainc << "\n";
		for (int a = 0; a < a2; a += ainc) {
			//std::cout << "alpha=" << (a*da*Fs) << " Hz  |S|=" << (SCo[kIdx + nf * a]) << "\n";
		}
		for (int i = 0; i < a2; ++i) {
			EES[i] = 0;
			float s = 0;
			for (int k = 0; k < nf; ++k) {
				s = s + SCo[k + nf * i];
			}
			EES[i] = s / (float)nf;
		}
		for (int a = 0; a < a2; a += ainc) {
			std::cout << "alpha=" << (a*da*Fs) << " Hz  ES (Standard ACP) =" << (EES[a]) << "\n";
		}

		//FSCoh(x, da, amax, L, Noverlap, Nw, 0, SCo);//FastACP
		//AFAM(x, da, amax, Lfam, Nw, 0, SCo);
		FastSCP(x, da, L, Nw, Nb, amax, 0, SCo);
		ainc = (int)std::llround(f2 / (da*Fs));
		std::cout << "Spectral Freq index=" << kIdx << ", a inc index:" << ainc << "\n";
		for (int a = 0; a < a2; a += ainc) {
			//std::cout << "alpha=" << (a*da*Fs) << " Hz  |S|=" << (SCo[kIdx + nf * a]) << "\n";
		}
		for (int i = 0; i < a2; ++i) {
			EES[i] = 0;
			float s = 0;
			for (int k = 0; k < nf; ++k) {
				s = s + SCo[k + nf * i];
			}
			EES[i] = s / (float)nf;
		}
		for (int a = 0; a < a2; a += ainc) {
			std::cout << "alpha=" << (a*da*Fs) << " Hz  ES (other method) =" << (EES[a]) << "\n";
		}
	}

	
	if (checkid == 2) { //2: check run time for alpha 0.0325 to 0.5 L=2^17
		for (runid = 0; runid < 20; ++runid) {
			std::cout << "Run id " << runid <<":\n";
			amax = 0.03125;
			for (ia = 0; ia < 5; ia++) { //to 5
				a2 = (int)(amax*L);
				std::cout << "amax, alpha count: " << amax << ", " << a2 << "\n";
				
				for (i = 0; i < Lfam; i++) {
					x[i] = (1.0f + sin(2.0 * pi * i * 10 / Nw)) * (float)rand() / 32767.0f;
				}

				//auto start = std::chrono::high_resolution_clock::now();		
				//SCoh(x, da, a1, a2, L, Noverlap, Nw, 0, SCo);
				//auto finish = std::chrono::high_resolution_clock::now();
				//std::chrono::duration<double> elapsed1 = finish - start;
				//std::cout << "Elapsed time normal ACP: " << elapsed1.count() << " s\n";

				auto start2 = std::chrono::high_resolution_clock::now();
				FAM(x, da, amax, Lfam, Nw, 0, SCo);
				auto finish2 = std::chrono::high_resolution_clock::now();
				std::chrono::duration<double> elapsed2 = finish2 - start2;
				std::cout << "Elapsed time FAM: " << elapsed2.count() << " s\n";
				T1[ia][runid] = elapsed2.count();

				auto start3 = std::chrono::high_resolution_clock::now();
				AntoniFSC2(x, da, amax, L, Nw, 0, SCo);
				auto finish3 = std::chrono::high_resolution_clock::now();
				std::chrono::duration<double> elapsed3 = finish3 - start3;
				std::cout << "Elapsed time Antoni FSC2: " << elapsed3.count() << " s\n";
				T2[ia][runid] = elapsed3.count();

				auto start4 = std::chrono::high_resolution_clock::now();
				FDirSC(x, da, amax, L, Nw, 0, SCo);
				auto finish4 = std::chrono::high_resolution_clock::now();
				std::chrono::duration<double> elapsed4 = finish4 - start4;
				std::cout << "Elapsed time Fast Dir SC: " << elapsed4.count() << " s\n";
				T3[ia][runid] = elapsed4.count();

				auto start5 = std::chrono::high_resolution_clock::now();
				FSCoh(x, da, amax, L, Noverlap, Nw, 0, SCo);
				auto finish5 = std::chrono::high_resolution_clock::now();
				std::chrono::duration<double> elapsed5 = finish5 - start5;
				std::cout << "Elapsed time FastACP: " << elapsed5.count() << " s\n";
				T4[ia][runid] = elapsed5.count();

				auto start6 = std::chrono::high_resolution_clock::now();
				AFAM(x, da, amax, Lfam, Nw, 0, SCo);
				auto finish6 = std::chrono::high_resolution_clock::now();
				std::chrono::duration<double> elapsed6 = finish6 - start6;
				std::cout << "Elapsed time AFAM: " << elapsed6.count() << "s \n";
				T5[ia][runid] = elapsed6.count();

				auto start7 = std::chrono::high_resolution_clock::now();
				FastSCP(x, da, L, Nw, Nb, amax, 0, SCo);
				auto finish7 = std::chrono::high_resolution_clock::now();
				std::chrono::duration<double> elapsed7 = finish7 - start6;
				std::cout << "Elapsed time FastSCP: " << elapsed7.count() << "s \n";
				T6[ia][runid] = elapsed7.count();


				amax = 2.0f*amax;
				if (amax > 0.5) break;
				std::cout << "============================================================ \n";
			}
		}
		printStat(T1, "FAM",5,20);
		printStat(T2, "Antoni FSC",5,20);
		printStat(T3, "Fast Dir",5,20);
		printStat(T4, "FastACP",5,20);
		printStat(T5, "AFAM",5,20);
		printStat(T6, "FstSCP",5,20);
	}
	if (checkid == 3) { //3: check run time for alpha 0.0325,  L=2^17 to 2^21
		amax = 0.0625;		
		for (ia = 0; ia < 5; ia++) { //to 5
			L = 1 << (17 + ia);
			std::cout << "Signal Length L: " << L << "\n";
			a2 = (int)(amax*L);
			std::cout << "amax, alpha count: " << amax << ", " << a2 << "\n";
			Lfam = L + Nw - Rfam; //used to get ready power-of-two KN since KN = (int)((L - Nw + R) / R);
			delete[] x;
			//x.resize(Lfam, 0.0);
			x = new float[Lfam + 32]();//putting () will initilize to zero all elements
		    for (runid = 0; runid < 20; ++runid) {
			 std::cout << "Run id " << runid << ":\n";									

				for (i = 0; i < Lfam; i++) {
					x[i] = (1.0f + sin(2.0 * pi * i * 10 / Nw)) * (float)rand() / 32767.0f;
				}
				//auto start = std::chrono::high_resolution_clock::now();		
				//SCoh(x, da, a1, a2, L, Noverlap, Nw, 0, SCo);
				//auto finish = std::chrono::high_resolution_clock::now();
				//std::chrono::duration<double> elapsed1 = finish - start;
				//std::cout << "Elapsed time normal ACP: " << elapsed1.count() << " s\n";

				auto start2 = std::chrono::high_resolution_clock::now();
				FAM(x, da, amax, Lfam, Nw, 0, SCo);
				auto finish2 = std::chrono::high_resolution_clock::now();
				std::chrono::duration<double> elapsed2 = finish2 - start2;
				std::cout << "Elapsed time FAM: " << elapsed2.count() << " s\n";
				T1[ia][runid] = elapsed2.count();

				auto start3 = std::chrono::high_resolution_clock::now();
				AntoniFSC2(x, da, amax, L, Nw, 0, SCo);
				auto finish3 = std::chrono::high_resolution_clock::now();
				std::chrono::duration<double> elapsed3 = finish3 - start3;
				std::cout << "Elapsed time Antoni FSC2: " << elapsed3.count() << " s\n";
				T2[ia][runid] = elapsed3.count();

				auto start4 = std::chrono::high_resolution_clock::now();
				FDirSC(x, da, amax, L, Nw, 0, SCo);
				auto finish4 = std::chrono::high_resolution_clock::now();
				std::chrono::duration<double> elapsed4 = finish4 - start4;
				std::cout << "Elapsed time Fast Dir SC: " << elapsed4.count() << " s\n";
				T3[ia][runid] = elapsed4.count();

				auto start5 = std::chrono::high_resolution_clock::now();
				FSCoh(x, da, amax, L, Noverlap, Nw, 0, SCo);
				auto finish5 = std::chrono::high_resolution_clock::now();
				std::chrono::duration<double> elapsed5 = finish5 - start5;
				std::cout << "Elapsed time FastACP: " << elapsed5.count() << " s\n";
				T4[ia][runid] = elapsed5.count();

				auto start6 = std::chrono::high_resolution_clock::now();
				AFAM(x, da, amax, Lfam, Nw, 0, SCo);
				auto finish6 = std::chrono::high_resolution_clock::now();
				std::chrono::duration<double> elapsed6 = finish6 - start6;
				std::cout << "Elapsed time AFAM: " << elapsed6.count() << "s \n";
				T5[ia][runid] = elapsed6.count();

				auto start7 = std::chrono::high_resolution_clock::now();
				FastSCP(x, da, L, Nw, Nb, amax, 0, SCo);
				auto finish7 = std::chrono::high_resolution_clock::now();
				std::chrono::duration<double> elapsed7 = finish7 - start6;
				std::cout << "Elapsed time FastSCP: " << elapsed7.count() << "s \n";
				T6[ia][runid] = elapsed7.count();								
				std::cout << "============================================================ \n";
			}
		}
		printStat(T1, "FAM", 5,20);
		printStat(T2, "Antoni FSC", 5,20);
		printStat(T3, "Fast Dir", 5,20);
		printStat(T4, "FastACP", 5,20);
		printStat(T5, "AFAM", 5,20);
		printStat(T6, "FstSCP", 5,20);
	}
	delete[] iPos;
	delete[] x;
	delete[] SCo;
	delete[] EES;
	
	getchar();
    return 0;
}

