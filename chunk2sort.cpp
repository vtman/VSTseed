//D:\Temp2\Genome\ref48 D:\Temp2\Genome\ref48sorted D:\Temp2\Genome\48.log 48 1000000 6

#include <fstream>

//#include <iostream>
//#include <chrono>

#include <omp.h>

#define nCut 20

#define nFiles 256

int findPosition(char* text, int n, char sym);

class C2S {
public:
	C2S();
	~C2S();

	char inputFolder[1000], outputFolder[1000], inputTemplate[1000], outputTemplate[1000];
	char logFile[1000];

	long long** nTotal;
	
	int nbits, nMax, nthreads;

	int processFiles(char* inFile, char* outFile, long long *nT, int m);
	int startProcessing();
	int setSettings(int n, char * va[]);
	int printUsage();
	int printInputParameters();

};

C2S::C2S() {
	nTotal = nullptr;
	nthreads = 1;
}

C2S::~C2S() {
	if (nTotal != nullptr) {
		for (int i = 0; i < nthreads; i++) {
			if (nTotal[i] != nullptr) {
				free(nTotal[i]); nTotal[i] = nullptr;
			}
		}
		free(nTotal); nTotal = nullptr;
	}
}

int C2S::printUsage() {
	printf("Usage:\n");
	printf("\t1) Input folder\n");
	printf("\t2) Output folder\n");
	printf("\t3) Log file\n");
	printf("\t4) Weight of the seed\n");
	printf("\t5) Max same\n");
	printf("\t6) Number of threads\n");
	printf("\n");
	return 0;
}

int C2S::printInputParameters() {
	printf("1) Input folder: \"%s\"\n", inputFolder);
	printf("2) Output folder: \"%s\"\n", outputFolder);
	printf("3) Log file: \"%s\"\n", logFile);
	printf("4) Weight of the seed: %i\n", nbits);
	printf("5) Max same: %i\n", nMax);
	printf("6) Number of threads: %i\n", nthreads);
	printf("\n");
	return 0;
}

bool isMore(unsigned char* s1, unsigned char* s2, int nb) {
	for (int i = nb - 1; i >= 0; i--) {
		if (s1[i] == s2[i]) continue;
		if (s1[i] > s2[i]) return true;
		break;
	}
	return false;
}
void swapRecords(unsigned char* s1, unsigned char* s2, unsigned char* vtemp, int nb) {
	memcpy(vtemp, s1, nb * sizeof(char));
	memcpy(s1, s2, nb * sizeof(char));
	memcpy(s2, vtemp, nb * sizeof(char));
}

int sortChunk(unsigned char* v1, unsigned char* v2, unsigned char *vtemp, int indStart, int ntot, int nb, int** pposLocal, int** pncLocal, int iLevel) {
	bool t;
	unsigned char* v_loc, u;
	int* ncLocal, *posLocal;
	v_loc = v1 + indStart * nb;
	if (ntot < nCut || iLevel == 0) {
		
		do {
			t = false;
			
			for (int i = 0; i < ntot - 1; i++) {
				if (isMore(v_loc + i * nb, v_loc + (i + 1) * nb, nb)) {
					t = true;
					swapRecords(v_loc + i * nb, v_loc + (i + 1) * nb, vtemp, nb);
				}
			}
		} while (t);
	}
	else {
		ncLocal = pncLocal[iLevel];
		posLocal = pposLocal[iLevel];
		memcpy(v2, v_loc, ntot * nb * sizeof(char));
		for (int i = 0; i < 256; i++) {
			ncLocal[i] = 0;
		}

		for (int i = 0; i < ntot; i++) {
			ncLocal[v2[i * nb + iLevel]]++;
		}
		posLocal[0] = 0;
		for (int i = 0; i < 256; i++) {
			posLocal[i + 1] = posLocal[i] + ncLocal[i];
		}

		for (int i = 0; i < 256; i++) {
			ncLocal[i] = 0;
		}

		for (int i = 0; i < ntot; i++) {
			u = v2[i * nb + iLevel];
			memcpy(v_loc + (posLocal[u] + ncLocal[u]) * nb, v2 + i * nb, nb * sizeof(char));
			ncLocal[u]++;
		}

		for (int i = 0; i < 256; i++) {
			if (ncLocal[i] <= 1) continue;
			sortChunk(v1, v2, vtemp, indStart + posLocal[i], ncLocal[i], nb, pposLocal, pncLocal, iLevel - 1);
		}

	}
	return 0;
}

int C2S::processFiles(char* inFile, char* outFile, long long *nT, int m) {
	int ncGlobal[65536], **pncLocal, nrec, posGlobal[65537], **pposLocal;
	unsigned short u1;
	FILE* fi, *fo;

	unsigned char* vi, * vo, *vtemp;

	long long fileSize;
	int nb;

	nb = 2 * (nbits / 8) + 3;
	printf("Input file: \"%s\"\n", inFile);
		
	fi = fopen(inFile, "rb");
	if (fi == nullptr) {
		printf("Error: cannot open input file \"%s\"\n", inFile);
		return -1;
	}

	_fseeki64(fi, 0, SEEK_END);
	fileSize = _ftelli64(fi);
	
	_fseeki64(fi, 0, SEEK_SET);

	vi = (unsigned char*)malloc(sizeof(unsigned char) * fileSize);

	if (vi == nullptr) {
		printf("Error: cannot allocate memory (input)\n");
		fclose(fi); fi = nullptr;
		return -2;
	}

	vo = (unsigned char*)malloc(sizeof(unsigned char) * fileSize);

	if (vo == nullptr) {
		printf("Error: cannot allocate memory (output)\n");
		free(vi); vi = nullptr;
		fclose(fi); fi = nullptr;
		return -3;
	}

	fread(vi, sizeof(unsigned char), fileSize, fi);

	fclose(fi); fi = nullptr;


	for (int i = 0; i < 65536; i++) {
		ncGlobal[i] = 0;
	}

	nrec = fileSize / nb;

	for (int i = 0; i < nrec; i++) {
		u1 = *(unsigned short*)(vi + i * nb + nb - 2);
		ncGlobal[u1]++;
	}
		
	posGlobal[0] = 0;

	for (int i = 0; i < 65536; i++) {
		posGlobal[i + 1] = posGlobal[i] + ncGlobal[i];
	}

	for (int i = 0; i < 65536; i++) {
		ncGlobal[i] = 0;
	}

	for (int i = 0; i < nrec; i++) {
		u1 = *(unsigned short*)(vi + i * nb + nb - 2);
		memcpy(vo + (posGlobal[u1] + ncGlobal[u1]) * nb, vi + i * nb, nb * sizeof(char));
		ncGlobal[u1]++;
	}

	vtemp = (unsigned char*)malloc(sizeof(unsigned char) * nb);

	pposLocal = (int**)malloc(sizeof(int*) * (nb - 2));
	pncLocal = (int**)malloc(sizeof(int*) * (nb - 2));

	for (int i = 0; i < nb - 2; i++) {
		pposLocal[i] = (int*)malloc(sizeof(int) * 257);
		pncLocal[i] = (int*)malloc(sizeof(int) * 256);
	}


	for (int i = 0; i < 65536; i++) {
		sortChunk(vo, vi, vtemp, posGlobal[i], ncGlobal[i], nb, pposLocal, pncLocal, nb - 3);
	}


	for (int i = 0; i < nb - 2; i++) {
		free(pposLocal[i]); pposLocal[i] = nullptr;
		free(pncLocal[i]); pncLocal[i] = nullptr;
	}

	free(pncLocal); pncLocal = nullptr;
	free(pposLocal); pposLocal = nullptr;
	free(vtemp); vtemp = nullptr;

	int indS, nq;
	indS = 0;
	nq = 1;

	bool t;

	for (int i = 1; i < nrec; i++) {
		t = true;
		for (int j = 4; j < nb; j++) {
			if (vo[indS * nb + j] != vo[i * nb + j]) {
				t = false;
				break;
			}
		}
		if (t) {
			nq++;
		}
		else {
			if (nq < m) {
				nT[nq]++;
			}
			indS = i;
			nq = 1;
		}
	}

	nT[nq]++;

	fo = fopen(outFile, "wb");
	if (fo == nullptr) {
		printf("Error: cannot open output file \"%s\"\n", outFile);
		return -4;
	}

	fwrite(vo, sizeof(char), fileSize, fo);
	
	fclose(fo); fo = nullptr;

	free(vi); vi = nullptr;
	free(vo); vo = nullptr;

	return 0;
}

int C2S::startProcessing() {

	nTotal = (long long**)malloc(sizeof(long long*) * nthreads);
	for (int i = 0; i < nthreads; i++) {
		nTotal[i] = (long long*)malloc(sizeof(long long*) * nMax);
		memset(nTotal[i], 0, sizeof(long long) * nMax);
	}

	omp_set_num_threads(nthreads);

#pragma omp parallel
	{
		int tid;
		tid = omp_get_thread_num();

		char outputFile[1000], inputFile[1000];

#pragma omp for
		for (int i = 0; i < nFiles; i++) {
			sprintf(inputFile, inputTemplate, i);
			sprintf(outputFile, outputTemplate, i);
			if (processFiles(inputFile, outputFile, nTotal[tid], nMax) != 0) continue;
		}
	}

	for (int i = 1; i < nthreads; i++) {
		for (int j = 0; j < nMax; j++) {
			nTotal[0][j] += nTotal[i][j];
		}
	}

	FILE* flog;
	flog = fopen(logFile, "w");
	if (flog == nullptr) {
		printf("Error: cannot open log file %s\n", logFile);
		return -1;
	}

	for (int i = 0; i < nMax; i++) {
		if (nTotal[0][i] == 0)continue;
		fprintf(flog, "%i\t%lld\n", i, nTotal[0][i]);
	}

	fclose(flog); flog = nullptr;

	return 0;
}

int C2S::setSettings(int n, char* va[]) {
	if (n != 7) {
		printUsage();
		return -1;
	}

	sprintf(inputFolder, va[1]);
	sprintf(outputFolder, va[2]);
	sprintf(logFile, va[3]);
	nbits = atoi(va[4]);
	nMax = atoi(va[5]);
	nthreads = atoi(va[6]);

	printInputParameters();

	sprintf(inputTemplate, "%s/ref_%%03i.bin", inputFolder);
	sprintf(outputTemplate, "%s/ref_%%03i.bin", outputFolder);

	return 0;
}

int main(int argc, char* argv[]) {
	int ires;
	C2S* cs;

	//auto t1 = std::chrono::high_resolution_clock::now();

	cs = new C2S();
	ires = cs->setSettings(argc, argv);
	if (ires == 0) {
		ires = cs->startProcessing();
	}

	/*
	auto t2 = std::chrono::high_resolution_clock::now();
	auto int_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);

	std::chrono::duration<long, std::micro> int_usec = int_ms;

	std::cout << "f() took " << int_ms.count() << " whole milliseconds "
		<< "(which is " << int_usec.count() << " whole microseconds)" << std::endl;

		*/
	delete cs; cs = nullptr;
	return ires;
}
