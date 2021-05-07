//64 2 13 -1 104 "hex_64_2_13.txt" "pat_64_2_13.txt" "per_64_2_13.txt"

#include <fstream>
#include <immintrin.h>
#include <chrono>

#include "ipp.h"

using namespace std::chrono;

#define n128 2

#define randBlock 1000000

#define nExtraAttempts 200 

struct mpattern {
	__m128i mu[n128];
};

int shiftOneBit(mpattern m1, mpattern* m2);

class BComb {
public:
	BComb(int inOne, int inGaps, int inPeriod, int inMax, int inReadLen, char* ioutputHexFile, char* ioutputPatternFile, char* ioutputPeriodFile);
	~BComb();

	int posLevel[257];

	int find32();
	//int find32Given();

	bool isMultiple(unsigned int nn);

	int setParameters();
	int startProcessing();
	int findTestPatterns2();
	bool checkPattern(int nS);

	int printPattern(mpattern x);
	int myRand();

	mpattern* vpat, mCand[128 * n128];

	char outputHexFile[1000], outputPatternFile[1000], outputPeriodFile[1000];
	int nPeriod, rndPos, nReads, maxPatternLength;
	int nCombinations, nGaps, nRLen, nOne, nMax;

	bool readFromFile;

	FILE* fIn, * foH, * foPat, * foPer;

	Ipp16s* rVector;
	IppsRandUniState_16s* pRndObj;

	mpattern* vM;
};

int BComb::myRand() {
	rndPos++;
	if (rndPos >= randBlock) {
		rndPos = 0;
		unsigned int tnow = (unsigned int)(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
		ippsRandUniformInit_16s(pRndObj, 0/*low*/, nRLen - 1/*high*/, tnow ^ (*(unsigned int*)(rVector)));
		ippsRandUniform_16s(rVector, randBlock, pRndObj);
	}
	return rVector[rndPos];
}


int BComb::printPattern(mpattern x) {
	unsigned int a, b[4];
	int ncPer, ncPat;
	ncPer = 0;
	ncPat = 0;
	for (int s = 0; s < n128; s++) {
		b[0] = _mm_extract_epi32(x.mu[s], 0);
		b[1] = _mm_extract_epi32(x.mu[s], 1);
		b[2] = _mm_extract_epi32(x.mu[s], 2);
		b[3] = _mm_extract_epi32(x.mu[s], 3);
		for (int i = 0; i < 4; i++) {
			a = b[i];
			for (int k = 0; k < 32; k++) {
				if ((a & 1) == 1) {
					if (ncPat < nOne)fprintf(foPat, "1");
					ncPat++;
					if (ncPat == nOne)fprintf(foPat, "\n");
					if (ncPer < nPeriod)fprintf(foPer, "1");

				}
				else {
					if (ncPat < nOne)fprintf(foPat, "0");
					if (ncPer < nPeriod)fprintf(foPer, "0");
				}
				ncPer++;
				if (ncPer == nPeriod)fprintf(foPer, "\n");
				a = (a >> 1);
			}
		}
		fprintf(foH, "%08x\t%08x\t%08x\t%08x", b[0], b[1], b[2], b[3]);
		if (s == n128 - 1) {
			fprintf(foH, "\n");
		}
		else {
			fprintf(foH, "\t");
		}
	}
	//printf("\n");
	return 0;
}

int shiftOneBit(mpattern m1, mpattern* m2) {
	__m128i a1, a2, a3, a4;

	a4 = _mm_set_epi32(0, 0, 0, 0);

	for (int s = 0; s < n128; s++) {
		a1 = _mm_slli_epi32(m1.mu[s], 1);
		a2 = _mm_srli_epi32(m1.mu[s], 31);
		a3 = _mm_bslli_si128(a2, 4);
		m2->mu[s] = _mm_or_si128(a4, _mm_or_si128(a1, a3));
		a4 = _mm_bsrli_si128(_mm_srli_epi32(m1.mu[s], 31), 12);
	}
	return 0;
}


BComb::BComb(int inOne, int inGaps, int inPeriod, int inMax, int inReadLen, char* ioutputHexFile, char* ioutputPatternFile, char* ioutputPeriodFile) {
	nOne = inOne;
	nGaps = inGaps;
	nPeriod = inPeriod;
	nMax = inMax;
	nRLen = inReadLen;
	sprintf(outputHexFile, ioutputHexFile);
	sprintf(outputPatternFile, ioutputPatternFile);
	sprintf(outputPeriodFile, ioutputPeriodFile);

	vM = nullptr;
	fIn = nullptr;
	foH = nullptr;
	foPat = nullptr;
	foPer = nullptr;

	readFromFile = false;
	rVector = ippsMalloc_16s(randBlock);
	int sizeRndObj;
	ippsRandUniformGetSize_16s(&sizeRndObj);
	pRndObj = (IppsRandUniState_16s*)ippsMalloc_8u(sizeRndObj);

	unsigned int tnow = (unsigned int)(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now().time_since_epoch()).count());

	ippsRandUniformInit_16s(pRndObj, 0/*low*/, nRLen - 1/*high*/, tnow);
	ippsRandUniform_16s(rVector, randBlock, pRndObj);
	rndPos = 0;
	vpat = nullptr;
}

BComb::~BComb() {
	if (vM != nullptr) {
		free(vM); vM = nullptr;
	}
	if (fIn != nullptr) {
		fclose(fIn); fIn = nullptr;
	}
	if (foH != nullptr) {
		fclose(foH); foH = nullptr;
	}
	if (foPat != nullptr) {
		fclose(foPat); foPat = nullptr;
	}
	if (foPer != nullptr) {
		fclose(foPer); foPer = nullptr;
	}
	if (rVector != nullptr) {
		ippsFree(rVector); rVector = nullptr;
	}
	if (pRndObj != nullptr) {
		ippsFree(pRndObj); pRndObj = nullptr;
	}
	if (vpat != nullptr) {
		free(vpat); vpat = nullptr;
	}
}

int BComb::setParameters() {
	unsigned int a, b, c, d;
	maxPatternLength = __min(128 * n128, nRLen);

	nCombinations = 100000;

	vpat = (mpattern*)malloc(sizeof(mpattern) * nCombinations);
	if (vpat == nullptr) {
		printf("Error: cannot allocate memory for combinations.\n");
		return -1;
	}

	foH = fopen(outputHexFile, "w");
	foPat = fopen(outputPatternFile, "w");
	foPer = fopen(outputPeriodFile, "w");
	if (readFromFile) {
		fIn = fopen("in.txt", "r");

		nReads = 0;
		while (fscanf(fIn, "%08x\t%08x\t%08x\t%08x", &a, &b, &c, &d) == 4) {
			nReads++;
		}

		nReads /= n128;

		printf("Number of input candidates: %i\n", nReads);

		vM = (mpattern*)malloc(sizeof(mpattern) * nReads);
		if (vM == nullptr) {
			return -1;
		}

		rewind(fIn);

		for (int i = 0; i < nReads; i++) {
			for (int j = 0; j < n128; j++) {
				fscanf(fIn, "%08x\t%08x\t%08x\t%08x", &a, &b, &c, &d);
				vM[i].mu[j] = _mm_set_epi32(d, c, b, a);
			}
		}

		fclose(fIn); fIn = nullptr;
	}
	return 0;
}

int BComb::findTestPatterns2() {
	int ii, kk;
	int* ipos;
	unsigned int uOne[4];
	__m128i mOne;
	mpattern* vOne, mMask;
	bool t;

	vOne = (mpattern*)malloc(sizeof(mpattern) * 128 * n128);

	mOne = _mm_set_epi32(0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff);

	for (int k = 0; k < n128; k++) {
		for (int kkk = 0; kkk < 128; kkk++) {
			for (int j = 0; j < n128; j++) {
				vOne[k * 128 + kkk].mu[j] = mOne;
			}
		}
	}

	for (int k = 0; k < n128; k++) {
		for (int i = 0; i < 128; i++) {
			uOne[0] = 0xffffffff;
			uOne[1] = 0xffffffff;
			uOne[2] = 0xffffffff;
			uOne[3] = 0xffffffff;
			ii = i / 32;
			kk = i % 32;
			uOne[ii] ^= (1 << kk);
			vOne[k * 128 + i].mu[k] = _mm_set_epi32(uOne[3], uOne[2], uOne[1], uOne[0]);
		}
	}


	for (int i = 0; i < n128; i++) {
		mMask.mu[i] = mOne;
	}

	for (int k = 0; k < n128; k++) {
		for (int i = 0; i < 128; i++) {
			if (k * 128 + i < nRLen) continue;
			for (int s = 0; s < n128; s++) {
				mMask.mu[s] = _mm_and_si128(vOne[k * 128 + i].mu[s], mMask.mu[s]);
			}
		}
	}

	ipos = (int*)malloc(sizeof(int) * nGaps);

	for (int i = 0; i < nCombinations; i++) {
		vpat[i] = mMask;

		ipos[0] = myRand();
		for (int j = 1; j < nGaps; j++) {
			do {
				t = false;
				ipos[j] = myRand();
				for (int s = 0; s < j; s++) {
					if (ipos[s] == ipos[j]) {
						t = true;
						break;
					}
				}
			} while (t);
		}
		ippsSortAscend_32s_I(ipos, nGaps);

		for (int p = 0; p < nGaps; p++) {
			for (int s = 0; s < n128; s++) {
				vpat[i].mu[s] = _mm_and_si128(vpat[i].mu[s], vOne[ipos[p]].mu[s]);
			}
		}

	}

	free(ipos); ipos = nullptr;
	free(vOne); vOne = nullptr;
	return 0;
}

/*/
int BComb::find32Given() {
	unsigned int u[4 * n128], ua;
	int mt, iLevel, nS;
	int cTrue, nCand;
	bool isFound;
	mpattern mc;

	cTrue = 0;
	nCand = 0;

	mt = 0;
	for (int i = 0; i < nReads; i++) {

		mc = vM[i];

		for (int s = 0; s < n128; s++) {
			u[4 * s + 0] = _mm_extract_epi32(mc.mu[s], 0);
			u[4 * s + 1] = _mm_extract_epi32(mc.mu[s], 1);
			u[4 * s + 2] = _mm_extract_epi32(mc.mu[s], 2);
			u[4 * s + 3] = _mm_extract_epi32(mc.mu[s], 3);
		}

		for (int s = 4 * n128 - 1; s >= 0; s--) {
			iLevel = s * 32;
			if (u[s] == 0) continue;
			ua = u[s] >> 1;
			while (ua != 0) {
				iLevel++;
				ua >>= 1;
			}
			break;
		}

		nS = nRLen - iLevel;

		mCand[0] = mc;

		for (int j = 1; j < nS; j++) {
			shiftOneBit(mCand[j - 1], mCand + j);
		}

		isFound = checkPattern(nS);

		if (isFound) {
			nCand++;
			bool t;
			for (int m = 0; m < nAttempts; m++) {
				findTestPatterns2();
				t = checkPattern(nS);
				if (!t)break;
			}
			if (t) {
				cTrue++;
				printf("Mt: %i\n", mt);
				printPattern(mc);
				printf("%i\n", iLevel + 1);
			}
		}
		mt++;
		if (cTrue == nMax) return cTrue;
	}
	printf("Total candidates: %i\n", mt);
	printf("Found patterns: %i\n", cTrue);
	return cTrue;
}
*/

bool BComb::isMultiple(unsigned int nn) {
	unsigned int mm, i2, ut, uq, ucheck;
	mm = (unsigned int)(ceil(sqrt(double(nPeriod))));

	ut = 1 | (nn << 1);

	for (unsigned int i1 = 2; i1 <= mm; i1++) {
		i2 = nPeriod / i1;
		if (i2 * i1 != nPeriod) continue;

		uq = (ut << (32 - i1)) >> (32 - i1);
		ucheck = 0;
		for (unsigned int i = 0; i < i2; i++) {
			ucheck |= (uq << (i1 * i));
		}
		if (ucheck == ut) return true;

		uq = (ut << (32 - i2)) >> (32 - i2);
		ucheck = 0;
		for (unsigned int i = 0; i < i1; i++) {
			ucheck |= (uq << (i2 * i));
		}
		if (ucheck == ut) return true;
	}
	return false;
}

int BComb::find32() {
	unsigned int u[4 * n128];
	int ind2, mt, ik, ik1, ik2, iLevel, nS, nT, nR, nr;
	int cTrue, nRT;
	unsigned int nC, n1, ua, ia;
	bool isFound, t;

	cTrue = 0;

	ind2 = nPeriod - 2;
	nC = 1 << ind2;
	mt = 0;
	for (unsigned int i = 0; i < nC; i++) {
		if (isMultiple(i)) continue;
		n1 = _mm_popcnt_u32(i);
		if (n1 > 30) continue;
		nR = nOne % (n1 + 1);
		nT = nOne / (n1 + 1);

		if (nR == 0) {
			nT--;
			nR = n1 + 1;
		}

		nRT = nPeriod * nT + 1;

		if (nRT > maxPatternLength) continue;

		if (nR > 1) {
			ia = nR - 1;
			while (ia > 0) {
				nRT++;
				ia >>= 1;
			}
			if (nRT > maxPatternLength) continue;
		}

		for (int s = 0; s < 4 * n128; s++) {
			u[s] = 0;
		}

		for (int k = 0; k < nT; k++) {
			ik = k * nPeriod;
			ik1 = ik / 32;
			ik2 = ik % 32;
			u[ik1] |= (1 << ik2);
		}

		for (int j = 0; j < ind2; j++) {
			if (((i >> j) & 1) == 0) continue;

			for (int k = 0; k < nT; k++) {
				ik = k * nPeriod + j + 1;
				ik1 = ik / 32;
				ik2 = ik % 32;
				u[ik1] |= (1 << ik2);
			}
		}

		ik = nT * nPeriod;
		ik1 = ik / 32;
		ik2 = ik % 32;
		u[ik1] |= (1 << ik2);
		nr = 1;

		if (nr < nR) {
			for (int j = 0; j < ind2; j++) {
				if (((i >> j) & 1) == 0) continue;
				ik = nT * nPeriod + j + 1;
				ik1 = ik / 32;
				ik2 = ik % 32;
				u[ik1] |= (1 << ik2);
				nr++;
				if (nr == nR) break;
			}
		}

		for (int s = 0; s < n128; s++) {
			mCand[0].mu[s] = _mm_set_epi32(u[4 * s + 3], u[4 * s + 2], u[4 * s + 1], u[4 * s + 0]);
		}

		for (int s = 4 * n128 - 1; s >= 0; s--) {
			iLevel = s * 32;
			if (u[s] == 0) continue;
			ua = u[s] >> 1;
			while (ua != 0) {
				iLevel++;
				ua >>= 1;
			}
			break;
		}

		if (iLevel > maxPatternLength) continue;

		nS = nRLen - iLevel;

		for (int j = 1; j < nS; j++) {
			shiftOneBit(mCand[j - 1], mCand + j);
		}

		isFound = checkPattern(nS);

		if (isFound) {
			for (int m = 0; m < nExtraAttempts; m++) {
				findTestPatterns2();
				t = checkPattern(nS);
				if (!t)break;
			}
			if (t) {
				cTrue++;
				printPattern(mCand[0]);
				printf("Mt: %i, %i\n", mt, iLevel + 1);
			}
		}
		mt++;
		if (cTrue == nMax) break;
	}
	printf("Total candidates: %i\n", mt);
	printf("Found patterns: %i\n", cTrue);
	return cTrue;
}

bool checkPair(mpattern mA, mpattern mB) {
	__m128i mC, mD, mR, m1, m2, m3, m4;
	for (int i = 0; i < n128; i++) {
		mC = _mm_or_si128(mA.mu[i], mB.mu[i]);
		mR = _mm_cmpeq_epi32(mA.mu[i], mC);
		if (i == 0) {
			mD = mR;
		}
		else {
			mD = _mm_and_si128(mD, mR);
		}
	}

	m1 = _mm_bsrli_si128(mD, 8);
	m2 = _mm_and_si128(mD, m1);
	m3 = _mm_bsrli_si128(m2, 4);
	m4 = _mm_and_si128(m2, m3);

	return !(_mm_extract_epi32(m4, 0) == 0);
}

bool BComb::checkPattern(int nS) {
	bool isOk;

	for (int i = 0; i < nCombinations; i++) {
		isOk = false;

		for (int j = 0; j < nS; j++) {
			isOk = checkPair(vpat[i], mCand[j]);
			if (isOk) break;
		}

		if (!isOk) {
			return false;
		}
	}

	return true;
}

int BComb::startProcessing() {
	int ires;
	if (setParameters() == -1) return -2;
	if (findTestPatterns2() != 0) return -1;

	//if (readFromFile) {
//		ires = find32Given();
	//}
	//else {
		ires = find32();
//	}

	return ires;
}

int main(int argc, char* argv[]) {
	int ires;
	int nOne, nGaps, nPeriod, nMax, nReadLen;
	char outputHexFile[1000], outputPatternFile[1000], outputPeriodFile[1000];

	nOne = atoi(argv[1]);
	nGaps = atoi(argv[2]);
	nPeriod = atoi(argv[3]);
	nMax = atoi(argv[4]);
	nReadLen = atoi(argv[5]);
	sprintf(outputHexFile, argv[6]);
	sprintf(outputPatternFile, argv[7]);
	sprintf(outputPeriodFile, argv[8]);

	BComb* bc;
	bc = new BComb(nOne, nGaps, nPeriod, nMax, nReadLen, outputHexFile, outputPatternFile, outputPeriodFile);
	ires = bc->startProcessing();

	delete bc; bc = nullptr;

	if (ires < 0) return -1;

	return ires;
}
