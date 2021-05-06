// C:\Temp2\Genome\Ref37\GRCh38.p13.genome.fa C:\Temp2\Genome\Ref37\human.bin C:\Temp2\Genome\Ref37\human.txt


#include <fstream>

#define maxN 64
#define n32 10000000
#define nblock (32*n32)

int findPosition(char* text, int n, char sym);

class FA2B {
public:
	FA2B();
	~FA2B();

	char inputFile[1000], outputBinFile[1000], outputIndexFile[1000], mline[500], * bText, * bText2, *vGap;
	int nGap;
	FILE* fi, * foB, * foI;
	long long fileSize;
	unsigned int* uOut;

	int allocateMemory();
	int openFiles();
	int closeFiles();
	int startProcessing();
	int setSettings(int n, char * va[]);
	int printUsage();
	int printInputParameters();
	int freeMemory();
};

int findPosition(char* text, int n, char sym) {
	for (int i = 0; i < n; i++) {
		if (text[i] == sym) return i;
	}
	return -1;
}


FA2B::FA2B() {
	fi = nullptr;
	foB = nullptr;
	foI = nullptr;

	bText = nullptr;
	bText2 = nullptr;
	vGap = nullptr;
	uOut = nullptr;
}

FA2B::~FA2B() {
	if (fi != nullptr) {
		fclose(fi); fi = nullptr;
	}
	if (foB != nullptr) {
		fclose(foB); foB = nullptr;
	}
	if (foI != nullptr) {
		fclose(foI); foI = nullptr;
	}

	if (bText != nullptr) {
		free(bText); bText = nullptr;
	}
	if (bText2 != nullptr) {
		free(bText2); bText2 = nullptr;
	}
	if (uOut != nullptr) {
		free(uOut); uOut = nullptr;
	}
	if (vGap != nullptr) {
		free(vGap); vGap = nullptr;
	}
}

int FA2B::closeFiles() {
	fclose(fi); fi = nullptr;
	fclose(foB); foB = nullptr;
	fclose(foI); foI = nullptr;
	return 0;
}

int FA2B::printUsage() {
	printf("Usage:\n");
	printf("\t1) Input FASTA file\n");
	printf("\t2) Output binary file\n");
	printf("\t3) Output index file\n");
	printf("\n");
	return 0;
}

int FA2B::printInputParameters() {
	printf("1) Input FASTA file: \"%s\"\n", inputFile);
	printf("2) Output binary file: \"%s\"\n", outputBinFile);
	printf("3) Output index file: \"%s\"\n", outputIndexFile);
	printf("\n");
	return 0;
}

int FA2B::openFiles() {
	foB = fopen(outputBinFile, "wb");
	if (foB == nullptr) {
		printf("Error: cannot write to file \"%s\"\n", outputBinFile);
		return -1;
	}

	foI = fopen(outputIndexFile, "w");
	if (foI == nullptr) {
		printf("Error: cannot write to file \"%s\"\n", outputIndexFile);
		return -2;
	}

	fi = fopen(inputFile, "rb");
	if (fi == nullptr) {
		printf("Error: cannot read from file \"%s\"\n", inputFile);
		return -3;
	}

	_fseeki64(fi, 0, SEEK_END);
	fileSize = _ftelli64(fi);

	printf("Input file size: %lld\n", fileSize);

	return 0;
}

int FA2B::allocateMemory() {
	uOut = (unsigned int *)malloc(sizeof(unsigned int) * 4 * n32);
	bText = (char*)malloc(sizeof(char) * nblock);
	bText2 = (char*)malloc(sizeof(char) * nblock);
	vGap = (char*)malloc(sizeof(char) * 16 * nGap);

	if (uOut == nullptr || bText == nullptr || bText2 == nullptr) {
		printf("Error: cannot allocate memory\n");
		return -1;
	}

	memset(vGap, 0, 16 * nGap * sizeof(char));

	return 0;
}

int FA2B::freeMemory() {
	free(bText); bText = nullptr;
	free(bText2); bText2 = nullptr;
	free(uOut); uOut = nullptr;
	free(vGap); vGap = nullptr;
	return 0;
}

int FA2B::startProcessing() {
	int nglen, nC, iki, nt, nt32, pk;
	int index1, index2, iu, index3, bCount;

	long long posStart, posEnd;
	long long groupStart, groupEnd, blockStart, blockEnd, posNN, outPos;
	unsigned int uA, uC, uG, uT, uu[32];


	if (openFiles() != 0) return -1;
	if (allocateMemory() != 0) return -2;

	for (int i = 0; i < 32; i++) {
		uu[i] = 1 << i;
	}
	
	posStart = 0;

	uA = 0;
	uC = 0;
	uG = 0;
	uT = 0;

	groupEnd = -1;

	fwrite(vGap, sizeof(char), 16 * nGap, foB);

	do {
		_fseeki64(fi, posStart, SEEK_SET);
		posEnd = posStart + nblock;
		if (posEnd >= fileSize) {
			posEnd = fileSize - 1;
		}

		bCount = 0;

		fread(bText, sizeof(char), posEnd - posStart, fi);

		index1 = findPosition(bText, posEnd - posStart, '>');
		if (index1 < 0) {
			printf("Error: cannot find \'>\'\n");
			break;
		}

		index2 = findPosition(bText + 1, posEnd - posStart - 1, '>');

		if (index2 < 0 && posEnd < fileSize - 1) {
			printf("Error: buffer is too short\n");
			break;
		}

		

		if (index2 < 0) {
			index2 = posEnd - posStart;
			bText[index2 + 1] = '\n';
			bText[index2 + 2] = '>';
			index2-=2;
		}

		printf("%i\t%i\t", index1, index2);

		index2++;

		for (iu = 0; iu < 200; iu++) {
			if (bText[index1 + iu] == '\n' || bText[index1 + iu] == '\r') break;
		}

		for (iki = 0; iki < iu; iki++) {
			if (bText[index1 + iki] <= 32) break;
		}

		strncpy(mline, bText + index1 + 1, iki - 1);
		mline[iki - 1] = '\0';

		printf("%s\n", mline);

		while (bText[index1 + iu] == '\n' || bText[index1 + iu] == '\r') {
			iu++;
		}

		index1 += iu;

		index3 = index2;

		while (bText[index3] == '\n' || bText[index3] == '\r' || bText[index3] == '>') {
			index3--;
		}

		nglen = 0;

		for (int i = index1; i <= index3; i++) {
			if (bText[i] <= 32) continue;
			bText2[nglen] = bText[i];
			if (bText2[nglen] > 96) bText2[nglen] -= 32;
			nglen++;
		}

		groupStart = groupEnd + 1;
		groupEnd = groupStart + nglen - 1;

		blockStart = 0;

		while (bText2[nglen - 1] == 'N') {
			nglen--;
		}

		for (int i = 0; i < nglen; i++) {
			if (bText2[i] != 'N') break;
			blockStart++;
		}
		blockEnd = blockStart + 1;

		posNN = blockStart;

		for (; blockEnd < nglen; blockEnd++) {
			if (bText2[blockEnd] != 'N') {
				posNN = blockEnd;
				nC = 0;
				continue;
			}
			nC++;
			if (nC == maxN) {
				nt = posNN - blockStart + 1;
				nt32 = nt / 32;
				if (nt % 32 > 0) nt32++;
				for (int p = 0; p < nt32; p++) {
					uA = 0;
					uC = 0;
					uG = 0;
					uT = 0;
					for (int k = 0; k < 32; k++) {
						pk = 32 * p + k;
						if (pk >= nt) break;
						switch (bText2[blockStart + pk]) {
						case 'A': uA |= uu[k]; break;
						case 'C': uC |= uu[k]; break;
						case 'G': uG |= uu[k]; break;
						case 'T': uT |= uu[k]; break;
						}
					}
					uOut[4 * p + 0] = uA;
					uOut[4 * p + 1] = uC;
					uOut[4 * p + 2] = uG;
					uOut[4 * p + 3] = uT;
				}
				outPos = _ftelli64(foB);
				fwrite(uOut, sizeof(unsigned int), 4 * nt32, foB);
				fwrite(vGap, sizeof(char), 16 * nGap, foB);

				fprintf(foI, "%s\t%i\t%lld\t%lld\t%i\t", mline, bCount, outPos, blockStart, nt);
				for (int j = 0; j < 32; j++) {
					fprintf(foI, "%1c", bText2[blockStart + j]);
				}
				fprintf(foI, "\t");
				for (int j = -31; j <= 0; j++) {
					fprintf(foI, "%1c", bText2[posNN + j]);
				}
				fprintf(foI, "\n");

				for (blockStart = blockEnd; blockStart < nglen; blockStart++) {
					if (bText2[blockStart] != 'N') break;
				}
				nC = 0;
				bCount++;
				if (blockStart >= nglen) break;
				blockEnd = blockStart;
			}

		}

		nt = posNN - blockStart + 1;
		nt32 = nt / 32;
		if (nt % 32 > 0) nt32++;
		for (int p = 0; p < nt32; p++) {
			uA = 0;
			uC = 0;
			uG = 0;
			uT = 0;
			for (int k = 0; k < 32; k++) {
				pk = 32 * p + k;
				if (pk >= nt) break;
				switch (bText2[blockStart + pk]) {
				case 'A': uA |= uu[k]; break;
				case 'C': uC |= uu[k]; break;
				case 'G': uG |= uu[k]; break;
				case 'T': uT |= uu[k]; break;
				}
			}
			uOut[4 * p + 0] = uA;
			uOut[4 * p + 1] = uC;
			uOut[4 * p + 2] = uG;
			uOut[4 * p + 3] = uT;
		}
		outPos = _ftelli64(foB);
		fwrite(uOut, sizeof(unsigned int), 4 * nt32, foB);
		fwrite(vGap, sizeof(char), 16 * nGap, foB);

		fprintf(foI, "%s\t%i\t%lld\t%lld\t%i\t", mline, bCount, outPos, blockStart, nt);
		for (int j = 0; j < 32; j++) {
			fprintf(foI, "%1c", bText2[blockStart + j]);
		}
		fprintf(foI, "\t");
		for (int j = -31; j <= 0; j++) {
			fprintf(foI, "%1c", bText2[posNN + j]);
		}
		fprintf(foI, "\n");

		posStart += index2;

		if (index2 <= 0) break;
	} while (posStart < fileSize - 2);
	
	closeFiles();
	freeMemory();
	
	return 0;
}

int FA2B::setSettings(int n, char* va[]) {
	if (n != 4) {
		printUsage();
		return -1;
	}

	sprintf(inputFile, va[1]);
	sprintf(outputBinFile, va[2]);
	sprintf(outputIndexFile, va[3]);

	printInputParameters();

	nGap = 8;

	return 0;
}

int main(int argc, char* argv[]) {
	int ires;
	FA2B* fb;

	fb = new FA2B();
	ires = fb->setSettings(argc, argv);
	if (ires == 0) {
		ires = fb->startProcessing();
	}

	delete fb; fb = nullptr;
	return ires;
}
