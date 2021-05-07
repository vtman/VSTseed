//C:\Temp2\Genome\Ref37\human.bin D:\Temp2\Genome\ref64

#include <fstream>

#include <immintrin.h>

//#include <iostream>
//#include <chrono>

//to be modified >>
#define ni32 3
#define no8 4
#define no32 1
//<<



#define ntot (no8*8)

#define no8m (2*no8-1)
#define nt8 (4 + no8m)

#define n32 1000000
#define nIblock (n32 + ni32 +1)
#define nOblock 1000000 
#define nFiles 256

int findPosition(char* text, int n, char sym);

class R2S {
public:
	R2S();
	~R2S();

	char inputRefFile[1000], outputFolder[1000], outputTemplate[1000], outputFile[1000], filePrefix[100];
	
	FILE* fi, *pfo[nFiles];
	long long fileSize;

	char* pvo[nFiles], *vi;
	unsigned char sigPos[no32*8 + 4];
	int ncount[nFiles];

	__m128i m1[ni32], m2[ni32], m3[ni32], mres[no32], mACCG[no32];

	int findACCG(__m128i* v1, __m128i* v2);
	int countOne(__m128i* vM);
	int allocateMemory();
	int openFiles();
	int closeFiles();
	int startProcessing();
	int setSettings(int n, char * va[]);
	int printUsage();
	int printInputParameters();
	int freeMemory();
};

int spaced2contig(__m128i* m, __m128i* res) {

	//32 ni=3
	__m128i c, t, s;
c = _mm_set1_epi32(0xb8083a77);
res[0] = _mm_and_si128(m[0], c);
c = _mm_set1_epi32(0x94800000);
t = _mm_and_si128(m[1], c);
s = _mm_srli_epi32(t, 5);
res[0] = _mm_or_si128(res[0], s);
c = _mm_set1_epi32(0x01400100);
t = _mm_and_si128(m[1], c);
res[0] = _mm_or_si128(res[0], t);
c = _mm_set1_epi32(0x08004091);
t = _mm_and_si128(m[1], c);
s = _mm_slli_epi32(t, 3);
res[0] = _mm_or_si128(res[0], s);
c = _mm_set1_epi32(0x00000042);
t = _mm_and_si128(m[1], c);
s = _mm_slli_epi32(t, 19);
res[0] = _mm_or_si128(res[0], s);
c = _mm_set1_epi32(0x0000000e);
t = _mm_and_si128(m[2], c);
s = _mm_slli_epi32(t, 13);
res[0] = _mm_or_si128(res[0], s);


/*
	
	// 64 3 ni=4
	__m128i c, t, s;
	c = _mm_set1_epi32(0xf712e5f7);
	res[0] = _mm_and_si128(m[0], c);
	c = _mm_set1_epi32(0x00808000);
	t = _mm_and_si128(m[2], c);
	s = _mm_srli_epi32(t, 12);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x00441200);
	t = _mm_and_si128(m[2], c);
	res[0] = _mm_or_si128(res[0], t);
	c = _mm_set1_epi32(0x00200400);
	t = _mm_and_si128(m[2], c);
	s = _mm_slli_epi32(t, 6);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x00000054);
	t = _mm_and_si128(m[3], c);
	s = _mm_slli_epi32(t, 17);
	res[0] = _mm_or_si128(res[0], s);

	c = _mm_set1_epi32(0xe5f712e5);
	res[1] = _mm_and_si128(m[1], c);
	c = _mm_set1_epi32(0x12006002);
	t = _mm_and_si128(m[2], c);
	res[1] = _mm_or_si128(res[1], t);
	c = _mm_set1_epi32(0x00010110);
	t = _mm_and_si128(m[2], c);
	s = _mm_slli_epi32(t, 11);
	res[1] = _mm_or_si128(res[1], s);
	c = _mm_set1_epi32(0x000001a3);
	t = _mm_and_si128(m[3], c);
	s = _mm_slli_epi32(t, 3);
	res[1] = _mm_or_si128(res[1], s);
	
	*/

	return 0;
}

int R2S::findACCG(__m128i* v1, __m128i* v2) {
	for (int i = 0; i < no32; i++) {
		v2[i] = _mm_or_si128(v1[i], _mm_bsrli_si128(v1[i], 4));
	}
	return 0;
}

int R2S::countOne(__m128i* vM) {
	__m128i ma, mb, mc, md;
	int n;
	n = 0;
	for (int i = 0; i < no32; i++) {
		ma = _mm_bsrli_si128(vM[i], 8);
		mb = _mm_or_si128(vM[i], ma);
		ma = _mm_bsrli_si128(mb, 4);
		mc = _mm_or_si128(mb, ma);
		n += _mm_popcnt_u32(_mm_extract_epi32(mc, 0));
	}

	return n;
}

R2S::R2S() {
	fi = nullptr;
	for (int i = 0; i < nFiles; i++) {
		pfo[i] = nullptr;
	}
	for (int i = 0; i < nFiles; i++) {
		pvo[i] = nullptr;
		ncount[i] = 0;
	}
	vi = nullptr;
}

R2S::~R2S() {
	if (fi != nullptr) {
		fclose(fi); fi = nullptr;
	}
	for (int i = 0; i < nFiles; i++) {
		free(pvo[i]); pvo[i] = nullptr;
	}
	if (vi != nullptr) {
		free(vi); vi = nullptr;
	}
}

int R2S::closeFiles() {
	fclose(fi); fi = nullptr;
	for (int i = 0; i < nFiles; i++) {
		fclose(pfo[i]); pfo[i] = nullptr;
	}
	return 0;
}

int R2S::printUsage() {
	printf("Usage:\n");
	printf("\t1) Input reference file (binary)\n");
	printf("\t2) Output folder\n");
	printf("\n");
	return 0;
}

int R2S::printInputParameters() {
	printf("1) Input reference file: \"%s\"\n", inputRefFile);
	printf("2) Output folder: \"%s\"\n", outputFolder);
	printf("\n");
	return 0;
}

int R2S::openFiles() {
	fi = fopen(inputRefFile, "rb");
	if (fi == nullptr) {
		printf("Error: cannot read from file \"%s\"\n", inputRefFile);
		return -3;
	}

	_fseeki64(fi, 0, SEEK_END);
	fileSize = _ftelli64(fi);

	_fseeki64(fi, 0, SEEK_SET);

	printf("Input file size: %lld\n", fileSize);

	sprintf(outputTemplate, "%s/%s%%03i.bin", outputFolder, filePrefix);

	for (int i = 0; i < nFiles; i++) {
		sprintf(outputFile, outputTemplate, i);
		pfo[i] = fopen(outputFile, "wb");
		if (pfo[i] == nullptr) {
			printf("Error: cannot open output file \"%s\"\n", outputFile);
			return -1;
		}
	}
	
	return 0;
}

int R2S::allocateMemory() {
	for (int i = 0; i < nFiles; i++) {
		pvo[i] = (char*)malloc(sizeof(char) * nOblock * nt8);
		if (pvo[i] == nullptr) {
			printf("Error: memory allocation for output data.\n");
			return -1;
		}
	}

	vi = (char*)malloc(sizeof(char) * 16 * nIblock);
	if (vi == nullptr) {
		printf("Error: cannot allocate memeory for input data.\n");
		return -2;
	}
	
	return 0;
}

int R2S::freeMemory() {
	free(vi); vi = nullptr;
	return 0;
}

int R2S::startProcessing() {
	long long posStart, posEnd, nStep, nb, n1, pos32, posTotal;
	unsigned char ufile;
	unsigned int uPos, u1, uTotal;
	__m128i* vM;
	__m128i ma, mb;

	if (openFiles() != 0) return -1;
	if (allocateMemory() != 0) return -2;

	nStep = n32 * 16;
	posEnd = fileSize - 16 * (ni32 + 1);

	nb = n32;

	uTotal = 0;

	for (posStart = 0; posStart <= posEnd; posStart += nStep) {
		pos32 = posStart / 16;
		printf("Position: %lld (%6.3f)\n", posStart, 100.0 * double(posStart) / double(fileSize));
		n1 = (posEnd - posStart) / 16 - ni32;
		if (n1 < nb) {
			nb = n1;
		}
		_fseeki64(fi, posStart, SEEK_SET);
		fread(vi, sizeof(char), 16 * (nb + ni32), fi);
		
		for (long long i = 0; i < nb; i++) {
			posTotal = pos32 + i;
			uPos = ((unsigned int)posTotal) * 32;
			vM = (__m128i*)(vi + 16 * i);
			for (int k = 0; k < ni32; k++) {
				m1[k] = _mm_load_si128(vM + k);
			}
			for (int k = 0; k < ni32 - 1; k++) {
				m2[k] = m1[k + 1];
			}
			m2[ni32 - 1] = _mm_load_si128(vM + ni32);


			for (int p = 0; p < 32; p++) {
				for (int k = 0; k < ni32; k++) {
					ma = _mm_srli_epi32(m1[k], p);
					mb = _mm_slli_epi32(m2[k], 32 - p);
					m3[k] = _mm_or_si128(ma, mb);
				}
				spaced2contig(m3, mres);
				if (countOne(mres) != ntot) continue;
				findACCG(mres, mACCG);
				memset(sigPos, 0, no32 * 8 + 4);

				*(unsigned int*)(sigPos) = uPos;

				for (int q = 0; q < no32; q++) {
					*(int *)(sigPos + 4 + 4*q) = _mm_extract_epi32(mACCG[q], 0);
				}

				for (int q = 0; q < no32; q++) {
					*(int*)(sigPos + 4 + no8 + 4 * q) = _mm_extract_epi32(mACCG[q], 1);
				}
				

				ufile = sigPos[nt8];

				memcpy(pvo[ufile] + ncount[ufile] * nt8, sigPos, nt8 * sizeof(char));
				ncount[ufile]++;

				if (ncount[ufile] == nOblock) {
					fwrite(pvo[ufile], sizeof(char), nt8 * nOblock, pfo[ufile]);
					uTotal += nOblock;
					ncount[ufile] = 0;
				}

				uPos++;
			}
			
		}
	}

	for (int i = 0; i < nFiles; i++) {
		fwrite(pvo[i], sizeof(char), nt8 * ncount[i], pfo[i]);
		uTotal += ncount[i];
	}

	printf("Total: %llu\n", uTotal);
	
	closeFiles();
	freeMemory();
	
	return 0;
}

int R2S::setSettings(int n, char* va[]) {
	if (n != 3) {
		printUsage();
		return -1;
	}

	sprintf(inputRefFile, va[1]);
	sprintf(outputFolder, va[2]);
	sprintf(filePrefix, "ref_");

	printInputParameters();

	return 0;
}

int main(int argc, char* argv[]) {
	int ires;
	R2S* fb;
	
	//auto t1 = std::chrono::high_resolution_clock::now();

	fb = new R2S();
	ires = fb->setSettings(argc, argv);
	if (ires == 0) {
		ires = fb->startProcessing();
	}

	delete fb; fb = nullptr;

	/*

	auto t2 = std::chrono::high_resolution_clock::now();
	auto int_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
	
	std::chrono::duration<long, std::micro> int_usec = int_ms;

	std::cout << "f() took " << int_ms.count() << " whole milliseconds "
		<< "(which is " << int_usec.count() << " whole microseconds)" << std::endl;
*/
	return ires;
}
