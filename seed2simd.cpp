//11011001111010100001101100111101010000110110011110101000011  D:\Genome\simd\testSimd.txt

#include <fstream>

#define maxSeq 10000

class Seed2Simd {
public:
	Seed2Simd(char *iSeed, char *iOutputFile);
	~Seed2Simd();

	int startProcessing();
	int printSeed();
	int formVector();
	int formMatrix();
	int printMatrix(int ny);
	int findSolution();
	int matrixRemoveSame(int* v, int ny, int nx);
	int checkMatrixAllPresent(int* v, int ny, int nx);
	int checkMatrixNoEmptyColumns(int* v, int ny, int nx);
	int initSeq(int n);
	int nextSeq(int n);
	int sortCheckSeq(int n);

	long long iAtt;

	FILE *fo;

	int * iSeq, * iSeqSorted, * iSaveSeq, *indexSorted;
	int indSol;

	int nRows, nCols, nLen, nZero;
	int nGaps, rowFirst, rowLast;
	int nExtraRows, iStep, nSeq, nBin;
	int nLength, nWeight, nLength32, nWeight32;

	int* iCurCol, * iCurLevel, * nRinC, * indexRow;
	int* iMatrix, * iShift, * iGrow, * iErow;
	int* indGap, * indExtra, * iPresent;
	int* iSmallMatrix, * indPosition;

	unsigned int* uBinE, *uBinC;
	unsigned int* uBRef, * uTest;

	bool* bGap, *bVector;
	
	char mSeed[1000], outputFile[1000], mLine[1000], mLine2[1000], *prAllIn, *prAllOut, pr1[32], pr2[32];
};

Seed2Simd::Seed2Simd(char *iSeed, char *iOutputFile) {
	fo = nullptr;
	iAtt = 0;
	prAllIn = nullptr;
	prAllOut = nullptr;
	iCurCol = nullptr;
	iCurLevel = nullptr;
	nRinC = nullptr;
	indexRow = nullptr;
	bVector = nullptr;
	indGap = nullptr;
	indExtra = nullptr;
	iMatrix = nullptr;
	iShift = nullptr;
	iGrow = nullptr;
	iErow = nullptr;
	iSeq = nullptr;
	indexSorted = nullptr;
	iSeqSorted = nullptr;
	iPresent = nullptr;
	indPosition = nullptr;
	iSmallMatrix = nullptr;
	uBinE = nullptr;
	uBinC = nullptr;
	uTest = nullptr;
	bGap = nullptr;
	uBRef = nullptr;
	iSaveSeq = nullptr;
	nSeq = 0;
	indSol = 0;
	sprintf(mSeed, iSeed);
	sprintf(outputFile, iOutputFile);
}


Seed2Simd::~Seed2Simd() {
	if (fo != nullptr) {
		fclose(fo); fo = nullptr;
	}
	if (prAllIn != nullptr) {
		free(prAllIn); prAllIn = nullptr;
	}
	if (prAllOut != nullptr) {
		free(prAllOut); prAllOut = nullptr;
	}
	if (bGap != nullptr) {
		free(bGap); bGap = nullptr;
	}
	if (uTest != nullptr) {
		free(uTest); uTest = nullptr;
	}
	if (uBRef != nullptr) {
		free(uBRef); uBRef = nullptr;
	}
	if (iCurCol != nullptr) {
		free(iCurCol); iCurCol = nullptr;
	}
	if (iCurLevel != nullptr) {
		free(iCurLevel); iCurLevel = nullptr;
	}
	if (nRinC != nullptr) {
		free(nRinC); nRinC = nullptr;
	}
	if (indexRow != nullptr) {
		free(indexRow); indexRow = nullptr;
	}
	if (uBinE != nullptr) {
		free(uBinE); uBinE = nullptr;
	}
	if (uBinC != nullptr) {
		free(uBinC); uBinC = nullptr;
	}
	if (indPosition != nullptr) {
		free(indPosition); indPosition = nullptr;
	}
	if (iSmallMatrix != nullptr) {
		free(iSmallMatrix); iSmallMatrix = nullptr;
	}
	if (bVector != nullptr) {
		free(bVector); bVector = nullptr;
	}
	if (indGap != nullptr) {
		free(indGap); indGap = nullptr;
	}
	if (indExtra != nullptr) {
		free(indExtra); indExtra = nullptr;
	}
	if (iMatrix != nullptr) {
		free(iMatrix); iMatrix = nullptr;
	}
	if (iShift != nullptr) {
		free(iShift); iShift = nullptr;
	}
	if (iGrow != nullptr) {
		free(iGrow); iGrow = nullptr;
	}
	if (iErow != nullptr) {
		free(iErow); iErow = nullptr;
	}
	if (iSeq != nullptr) {
		free(iSeq); iSeq = nullptr;
	}
	if (indexSorted != nullptr) {
		free(indexSorted); indexSorted = nullptr;
	}
	if (iSeqSorted != nullptr) {
		free(iSeqSorted); iSeqSorted = nullptr;
	}
	if (iPresent != nullptr) {
		free(iPresent); iPresent = nullptr;
	}
	if (iSaveSeq != nullptr) {
		free(iSaveSeq); iSaveSeq = nullptr;
	}
}

int Seed2Simd::printMatrix(int ny) {
	int* iM, ii, ntot, iLetter;
	int i1, i2;
	unsigned int uc;
	char mL;
	bool tt;

	indSol++;
	iM = iSmallMatrix + iStep * nCols * ny;

	ntot = 25 + 3 * (nCols - 1) + 4 * nCols;

	for (int i = 0; i < ntot; i++) {
		mLine[i] = '_';
	}
	mLine[ntot] = '\0';

	for (int i = 0; i < ntot; i++) {
		mLine2[i] = '=';
	}
	mLine2[ntot] = '\0';
	fprintf(fo, "%s\n\nSolution #%i\n\n", mLine2, indSol);

	for (int i = 0; i < ny; i++) {
		ii = indexSorted[i];
		fprintf(fo, "%4i | %3i | %3i | %3i | ", ii + 1, iGrow[ii] + 1, iErow[ii] + 1, iShift[ii]);
		for (int j = 0; j < nCols; j++) {
			if (iM[i * nCols + j] == 0) {
				fprintf(fo, "%4i", iM[i * nCols + j]);
			}
			else {
				fprintf(fo, "%4i", indExtra[iM[i * nCols + j]-1]);
			}
			if (j < nCols - 1) fprintf(fo, " & ");
		}
		fprintf(fo, "\n");
	}
	fprintf(fo, "\n");

	for (int i = 0; i < nLength32 * 32; i++) {
		prAllIn[i] = '0';
		prAllOut[i] = '0';
	}

	iLetter = 0;
	
	for (int i = 0; i < nWeight; i++) {
		if (!bVector[i]) continue;
		i1 = i / 32;
		i2 = i % 32;
		prAllIn[i] = 'A' + i1;
		prAllOut[i] = 'A' + i1;
	}
	
	iLetter += nWeight32;

	for (int i = 0; i < ny; i++) {
		ii = indexSorted[i];
		for (int j = 0; j < nCols; j++) {
			i1 = iM[i * nCols + j];
			if (i1 == 0)continue;
			i2 = indExtra[i1 - 1];
			
			prAllIn[i2] = 'A' + iLetter;
			prAllOut[indGap[j]] = 'A' + iLetter;
		}
		iLetter++;
	}

	fprintf(fo, "Input seed:\n\n");

	for (int i = 0; i < nLength32; i++) {
		for (int j = 0; j < 32; j++) {
			fprintf(fo, "%c", prAllIn[32 * i + j]);
			if (j == 7 || j == 15 || j == 23)fprintf(fo, " ");
		}
		fprintf(fo, "\n");
	}

	fprintf(fo, "\n");

	fprintf(fo, "Output seed:\n\n");

	for (int i = 0; i < nWeight32; i++) {
		for (int j = 0; j < 32; j++) {
			fprintf(fo, "%c", prAllOut[32 * i + j]);
			if (j == 7 || j == 15 || j == 23)fprintf(fo, " ");
		}
		fprintf(fo, "\n");
	}

	fprintf(fo, "\n%s\n", mLine);

	for (int i = 0; i < nWeight32; i++) {
		fprintf(fo, "\nRow %2i:", i + 1);
		for (int j = 0; j < 31; j++) {
			fprintf(fo, " ");
		}
		for (int j = 0; j < 32; j++) {
			fprintf(fo, "%c", prAllIn[32 * i + j]);
		}
		fprintf(fo,"\n\n");
		for (int k = 0; k < ny; k++) {
			ii = indexSorted[k];
			if (i != iGrow[ii]) continue;
			i1 = iShift[ii];
			i2 = iErow[ii];
			fprintf(fo, "Ext  %c:", 'A' + nWeight32 + k);
			for (int j = 0; j < 31+i1; j++) {
				fprintf(fo, " ");
			}
			for (int j = 0; j < 32; j++) {
				fprintf(fo, "%c", prAllIn[32 * i2 + j]);
			}
			fprintf(fo, "\n");
		}
	}

	fprintf(fo, "\n%s\n\n", mLine);
	
	for (int i = 0; i < nWeight32; i++) {
		fprintf(fo, "Row %2i:", i + 1);
		for (int j = 0; j < 31; j++) {
			fprintf(fo, " ");
		}
		for (int j = 0; j < 32; j++) {
			fprintf(fo, "%c", prAllIn[32 * i + j]);
		}
		fprintf(fo, "\n\n");
		for (int k = 0; k < ny; k++) {
			ii = indexSorted[k];
			if (i != iGrow[ii]) continue;
			i1 = iShift[ii];
			i2 = iErow[ii];
			mL = 'A' + nWeight32 + k;
			fprintf(fo, "Ext  %c:", mL);
			for (int j = 0; j < 31 + i1; j++) {
				fprintf(fo, " ");
			}
			for (int j = 0; j < 32; j++) {
				if (prAllIn[32 * i2 + j] == mL) {
					fprintf(fo, "%c", mL);
				}
				else {
					fprintf(fo, "0");
				}
			}
			fprintf(fo, "\n");
		}
	}

	fprintf(fo, "\n%s\n\n", mLine);

	fprintf(fo, "__m128i c, t, s;\n");

	for (int i = 0; i < nWeight32; i++) {
		uc = 0;
		for (int j = 0; j < 32; j++) {
			if (prAllIn[32 * i + j] == 'A' + i) {
				uc |= (1 << j);
			}
		}
		fprintf(fo, "c = _mm_set1_epi32(0x%08x);\n", uc);

		fprintf(fo, "res[%i] = _mm_and_si128(m[%i], c);\n", i, i);
		
		for (int k = 0; k < ny; k++) {
			ii = indexSorted[k];
			if (i != iGrow[ii]) continue;

			uc = 0;

			i1 = iShift[ii];
			i2 = iErow[ii];
			mL = 'A' + nWeight32 + k;
			
			for (int j = 0; j < nCols; j++) {
				if (iM[k * nCols + j] == 0) continue;
				uc |= (1 << (indExtra[iM[k * nCols + j] - 1] % 32));
			}
			fprintf(fo, "c = _mm_set1_epi32(0x%08x);\n", uc);
			
			fprintf(fo, "t = _mm_and_si128(m[%i], c);\n", i2);
			if (i1 == 0) {
				fprintf(fo, "res[%i] = _mm_or_si128(res[%i], t);\n", i, i);
			}else if (i1 > 0) {
				fprintf(fo, "s = _mm_slli_epi32(t, %i);\n", i1);
				fprintf(fo, "res[%i] = _mm_or_si128(res[%i], s);\n", i, i);
			}else {
				fprintf(fo, "s = _mm_srli_epi32(t, %i);\n", -i1);
				fprintf(fo, "res[%i] = _mm_or_si128(res[%i], s);\n", i, i);
			}
		}
		fprintf(fo, "\n");
	}

	fprintf(fo, "%s\n\n", mLine2);

	tt = false;

	for (int k = 0; k < nSeq; k++) {
		tt = true;
		for (int i = 0; i < ny; i++) {
			if (iSaveSeq[k*ny + i] == iSeqSorted[i])continue;
			tt = false;
			break;
		}
		if (tt) break;
	}

	if (!tt) {
		for (int i = 0; i < ny; i++) {
			iSaveSeq[nSeq * ny + i] = iSeqSorted[i];
		}
		nSeq++;
	}
		
	return 0;
}


int Seed2Simd::matrixRemoveSame(int *v, int ny, int nx) {
	bool t, tt;
	int nonz, ival, iCount, iCol;

	do {
		do {
			t = false;
			for (int i = 0; i < nx; i++) {
				nonz = 0;
				for (int j = 0; j < ny; j++) {
					if (v[j*nx + i] == 0)continue;
					nonz++;
					ival = v[j*nx + i];
					if (nonz == 2) break;
				}
				if (nonz == 0) return -1;

				if (nonz == 1) {
					for (int jj = 0; jj < ny; jj++) {
						for (int ii = 0; ii < nx; ii++) {
							if (ii == i) continue;
							if (v[jj*nx + ii] == ival) {
								v[jj*nx + ii] = 0;
								t = true;
							}
						}
					}
				}
			}
		} while (t);

		for (int i = 1; i <= nx; i++) {
			iPresent[i] = 0;
		}

		for (int i = 0; i < ny; i++) {
			for (int j = 0; j < nx; j++) {
				if (v[i * nx + j] == 0)continue;
				iPresent[v[i * nx + j]]++;
			}
		}

		for (int i = 1; i <= nx; i++) {
			if (iPresent[i] == 0) return -1;
		}

		tt = false;
		
		for (int i = 1; i <= nx; i++) {
			if (iPresent[i] != 1) continue;

			for (int ii = 0; ii < nx; ii++) {
				for (int jj = 0; jj < ny; jj++) {
					if (v[jj * nx + ii] != i)continue;
					for (int k = 0; k < ny; k++) {
						if (k == jj)continue;
						if (v[k * nx + ii] != 0) {
							iPresent[v[k * nx + ii]]--;
							v[k * nx + ii] = 0;
							tt = true;
						}
					}
				}
			}
		}
		
		
	} while (tt);

	for (int i = 0; i < nx; i++) {
		iPresent[i] = 0;
	}

	for (int i = 0; i < ny; i++) {
		for (int j = 0; j < nx; j++) {
			if (v[i * nx + j] == 0)continue;
			iPresent[j]++;
		}
	}

	iCount = 10000;
	iCol = -1;

	for (int i = 0; i < nx; i++) {
		if (iPresent[i] == 1) continue;
		if (iPresent[i] >= iCount) continue;
		iCount = iPresent[i];
		iCol = i;
	}

	if (iCol == -1) return -2;
	return iCol;
}

int Seed2Simd::checkMatrixAllPresent(int* v, int ny, int nx) {
	for (int i = 1; i <= nx; i++) {
		iPresent[i] = 0;
	}

	for (int i = 0; i < ny; i++) {
		for (int j = 0; j < nx; j++) {
			if (v[i * nx + j] == 0) continue;
			iPresent[v[i * nx + j]]++;
		}
	}
	
	for (int i = 1; i <= nx; i++) {
		if (iPresent[i] > 0)continue;
		return -1;
	}
	return 0;
}

int Seed2Simd::checkMatrixNoEmptyColumns(int* v, int ny, int nx) {
	int s;
	for (int j = 0; j < nx; j++) {
		s = 0;
		for(int i = 0; i < ny; i++){
			s += v[i * nx + j];
		}
		if (s == 0) return -1;
	}
	return 0;
}

int Seed2Simd::initSeq(int n) {
	int m, indd, nO, iL;
	int i1, i2;
	m = 0;

	for (int i = 0; i < nBin*(nCols+1); i++) {
		uTest[i] = 0;
	}

	indd = 0;
	iL = 0;
	do {
		iCurCol[m] = indd;
		iCurLevel[m] = iL;
		bGap[m] = true;
		iSeq[m] = indexRow[iCurCol[m] * nRows + iCurLevel[m]];
		nO = 0;
		for (int i = 0; i < nBin; i++) {
			uTest[(m + 1) * nBin + i] = uTest[m * nBin + i] | uBinC[iSeq[m] * nBin + i];
			nO += __popcnt(uTest[(m + 1) * nBin + i]);
		}
		iL = 0;
		m++;
		if (nO == nCols) break;
		if (m < n) {
			for (int k = 0; k < nCols; k++) {
				i1 = k / 32;
				i2 = k % 32;
				if (((uTest[m * nBin + i1] >> i2) & 1) == 1) continue;
				indd = k;
				break;
			}
		}
		else {
			do {
				m--;
				if (m == -1) return -1;
				indd = iCurCol[m];
				iL = iCurLevel[m] + 1;
			} while (iL == nRinC[iCurCol[m]]);
			
		}
	} while (true);
		

	bool t;

	for (; m < n; m++) {
		bGap[m] = false;
		iCurCol[m] = -1;
		iCurLevel[m] = -1;
		for (int k = 0; k < nRows; k++) {
			t = true;
			for (int s = 0; s < m; s++) {
				if (iSeq[s] == k) {
					t = false;
					break;
				}
			}
			if (t) {
				iSeq[m] = k;
				break;
			}
		}
		if (!t) return -1;
	}


	return 0;
}

int Seed2Simd::nextSeq(int n) {
	int m, indd, nO, i1, i2;
	bool t;
	m = n-1;

	indd = iCurCol[m];
	
	do {
		if (bGap[m]) {
			//printf("Here\n");
			iCurCol[m] = indd;
			iCurLevel[m]++;
			if (iCurLevel[m] == nRinC[iCurCol[m]]) {
				iCurLevel[m] = -1;
				m--;
				if (m == -1) return -1;
				indd = iCurCol[m];
			}
			else {

				iSeq[m] = indexRow[iCurCol[m] * nRows + iCurLevel[m]];
				nO = 0;
				for (int i = 0; i < nBin; i++) {
					uTest[(m + 1) * nBin + i] = uTest[m * nBin + i] | uBinC[iSeq[m] * nBin + i];
					nO += __popcnt(uTest[(m + 1) * nBin + i]);
				}
				m++;
				if (nO == nCols) break;

				if (m < n) {
					for (int k = 0; k < nCols; k++) {
						i1 = k / 32;
						i2 = k % 32;
						if (((uTest[m * nBin + i1] >> i2) & 1) == 1) continue;
						bGap[m] = true;
						indd = k;
						break;
					}
				}
				else {
					m--;
					if (m == -1) return -1;
				}
			}
		}
		else {
			do {
				iSeq[m]++;
				t = true;
				for (int s = 0; s < m; s++) {
					if (iSeq[s] == iSeq[m]) {
						t = false;
						break;
					}
				}
			} while (!t);
			
			if (iSeq[m] == nRows) {
				m--;
				if (m == -1) return -1;
				indd = iCurCol[m];
			}
			else {
				m++;
				break;
			}
		}
	} while (true);

	for (; m < n; m++) {
		bGap[m] = false;
		iCurCol[m] = -1;
		iCurLevel[m] = -1;
		for (int k = 0; k < nRows; k++) {
			t = true;
			for (int s = 0; s < m; s++) {
				if (iSeq[s] == k) {
					t = false;
					break;
				}
			}
			if (t) {
				iSeq[m] = k;
				break;
			}
		}
		if (!t) return -1;
	}
		
	return 0;
}

int Seed2Simd::sortCheckSeq(int n) {
	int nm, nm2;

	nm2 = -1;
	for (int i = 0; i < n; i++) {
		nm = 1000000;
		for (int k = 0; k < n; k++) {
			if (iSeq[k] > nm2 && iSeq[k] < nm) {
				nm = iSeq[k];
				indexSorted[i] = iSeq[k];
			}
		}
		iSeqSorted[i] = nm;
		nm2 = nm;
	}
	

	bool  tt;
	tt = false;

	for (int k = 0; k < nSeq; k++) {
		tt = true;
		for (int i = 0; i < n; i++) {
			if (iSaveSeq[k*n + i] == iSeqSorted[i])continue;
			tt = false;
			break;
		}
		if (tt) break;
	}

	if (tt) return -1;

	return 0;
}

int Seed2Simd::findSolution() {
	int ires, irow, icol, irowNew;
	int i1, i2, iFind;
	int* iM, *iM_old;
	
	unsigned int btest;

	bool isSolFound, tSimple;
	
	isSolFound = false;

	iSaveSeq = (int*)malloc(sizeof(int) * maxSeq * nCols);
	iSmallMatrix = (int*)malloc(sizeof(int) * nCols * nCols * nCols);
	indPosition = (int*)malloc(sizeof(int) * nCols);

	bGap = (bool*)malloc(sizeof(bool) * nCols);
	iSeq = (int*)malloc(sizeof(int) * nCols);
	indexSorted = (int*)malloc(sizeof(int) * nCols);
	iSeqSorted = (int*)malloc(sizeof(int) * nCols);
	iPresent = (int*)malloc(sizeof(int) * (nCols + 1));
	uTest = (unsigned int*)malloc(sizeof(unsigned int) * nBin * (nCols + 1));
	uBRef = (unsigned int*)malloc(sizeof(unsigned int) * nBin);

	iStep = 0;
		
	for (int i = 0; i < nBin; i++) {
		uBRef[i] = 0;
	}
	
	for (int i = 0; i < nCols; i++) {
		i1 = i / 32;
		i2 = i % 32;
		uBRef[i1] |= (1 << i2);
	}
		
	for (int nLevel = 1; nLevel <= nCols; nLevel++) {
		printf("Level: %i\n", nLevel);
		fprintf(fo, "Level: %i\n", nLevel);

		iFind = initSeq(nLevel);

		if (iFind == -1) continue;
		
		do {
			for (int i = 0; i < nBin; i++) {
				btest = 0;
			
				for (int j = 0; j < nLevel; j++) {
					btest |= uBinE[iSeq[j] * nBin + i];
				}

				tSimple = (btest == uBRef[i]);
				if (!tSimple)break;
			}
			if(tSimple){
				ires = sortCheckSeq(nLevel);
				if (ires == 0) {
					for (int j = 0; j < nLevel; j++) {
						memcpy(iSmallMatrix + j * nCols, iMatrix + iSeqSorted[j] * nCols, nCols * sizeof(int));
					}

					iStep = 0;
					iM = iSmallMatrix + iStep * nLevel * nCols;

					do {
						ires = matrixRemoveSame(iM, nLevel, nCols);
						if (ires >= 0) {
							irow = 0;
							for (int i = 0; i < nLevel; i++) {
								if (iM[i * nCols + ires] > 0) {
									irow = i;
									break;
								}
							}
							indPosition[iStep] = irow * nCols + ires;
							iM_old = iM;
							iM = iSmallMatrix + (iStep + 1) * nLevel * nCols;
							memcpy(iM, iM_old, sizeof(int) * nLevel * nCols);
							for (int i = 0; i < nLevel; i++) {
								iM[i * nCols + ires] = 0;
							}
							iM[indPosition[iStep]] = iM_old[indPosition[iStep]];
							iStep++;
						}
						else {
							if (ires == -2) {
								printf("Found\n");
								printMatrix(nLevel);
								isSolFound = true;
							}
							if (iStep == 0) break;
							do {
								iM = iSmallMatrix + iStep * nLevel * nCols;
								iM_old = iSmallMatrix + (iStep - 1) * nLevel * nCols;
								memcpy(iM, iM_old, sizeof(int) * nLevel * nCols);
								irow = indPosition[iStep - 1] / nCols;
								icol = indPosition[iStep - 1] % nCols;
								irowNew = -1;
								for (int i = irow + 1; i < nLevel; i++) {
									if (iM_old[i * nCols + icol] > 0) {
										irowNew = i;
										break;
									}
								}
								if (irowNew == -1) {
									iStep--;
								}
								else {
									indPosition[iStep - 1] = irowNew * nCols + icol;
									for (int i = 0; i < nLevel; i++) {
										iM[i * nCols + icol] = 0;
									}
									iM[indPosition[iStep - 1]] = iM_old[indPosition[iStep - 1]];
									break;
								}
							} while (iStep > 0);
						}

					} while (iStep > 0);
				}
			}

			iFind = nextSeq(nLevel);
		} while (iFind == 0);
		if (isSolFound) {
			printf("Solution found\n");
			break;
		}
	}

	free(uTest); uTest = nullptr;
	
	return 0;
}

int Seed2Simd::formMatrix() {
	int iGs, iGe, iEs, iEe, qGap, qExtra, qtest1;
	int shiftS, shiftE, n;
	int* iMatrix_loc;
	unsigned int* uBinE_loc, * uBinC_loc;
	bool t;
		
	nGaps = nLength - nWeight;

	fprintf(fo, "Number of gaps: %i\n", nGaps);

	n = 0;
	for (int i = 0; i < nWeight; i++) {
		if (bVector[i])continue;
		n++;
	}

	nCols = n;

	indGap = (int*)malloc(sizeof(int) * nCols);
	indExtra = (int*)malloc(sizeof(int) * nCols);

	n = 0;
	for (int i = 0; i < nWeight; i++) {
		if (bVector[i])continue;
		indGap[n] = i;
		n++;
	}

	for (int i = 0; i < n; i++) {
		fprintf(fo, "#%i: %i\n", i + 1, indGap[i]+1);
	}
	fprintf(fo, "\n");

	n = 0;
	for (int i = nWeight; i < nLength; i++) {
		if (!bVector[i])continue;
		indExtra[n] = i;
		n++;
	}

	fprintf(fo, "Extra 1s:\n");
	for (int i = 0; i < n; i++) {
		fprintf(fo, "#%i: %i\n", i + 1, indExtra[i]+1);
	}
	fprintf(fo, "\n");

	for (int ik = 0; ik < 2; ik++) {
		nRows = 0;
		iGs = 0;
		for (;;) {
			qGap = indGap[iGs] / 32;
			iGe = iGs + 1;
			do {
				if (iGe == nCols) break;
				if (indGap[iGe] / 32 > qGap) break;
				iGe++;
			} while (iGe < nCols);
			iGe--;
			if(ik == 0) fprintf(fo, "[%i, %i]\n", iGs+1, iGe+1);

			iEs = 0;
			for (;;) {
				qExtra = indExtra[iEs] / 32;
				iEe = iEs + 1;
				do {
					if (iEe == nCols) break;
					if (indExtra[iEe] / 32 > qExtra) break;
					iEe++;
				} while (iEe < nCols);
				iEe--;
				if (ik == 0) fprintf(fo, "\t[%i, %i]\n", iEs+1, iEe+1);

				shiftS = (indGap[iGs] % 32) - (indExtra[iEe] % 32);
				shiftE = (indGap[iGe] % 32) - (indExtra[iEs] % 32);

				for (int ish = shiftS; ish <= shiftE; ish++) {
					t = false;
					for (int ii = iGs; ii <= iGe; ii++) {
						qtest1 = indGap[ii] % 32 - ish;
						for (int jj = iEs; jj <= iEe; jj++) {
							if (qtest1 == indExtra[jj] % 32) {
								t = true;
								break;
							}
						}
						if (t) break;
					}
					if (!t) continue;
					if (ik == 1) {
						iMatrix_loc = iMatrix + nRows * nCols;
						uBinE_loc = uBinE + nRows * nBin;
						uBinC_loc = uBinC + nRows * nBin;
						for (int m = 0; m < nCols; m++) {
							iMatrix_loc[m] = 0;
						}
						
						for (int ii = iGs; ii <= iGe; ii++) {
							qtest1 = indGap[ii] % 32 - ish;
							for (int jj = iEs; jj <= iEe; jj++) {
								if (qtest1 == indExtra[jj] % 32) {
									iMatrix_loc[ii] = jj + 1;
									uBinE_loc[jj / 32] |= (1 << (jj % 32));
									uBinC_loc[ii / 32] |= (1 << (ii % 32));
								}
							}
						}

						iShift[nRows] = ish;
						iGrow[nRows] = qGap;
						iErow[nRows] = qExtra;
					}
					nRows++;
				}
				iEs = iEe + 1;
				if (iEs == nCols) break;
			}

			iGs = iGe + 1;
			if (iGs == nCols) break;
		}
		
		if (ik == 0) {
			fprintf(fo, "Number of rows in a matrix: %i\n\n", nRows);
			nBin = nCols / 32;
			if (nCols % 32 > 0)nBin++;
			uBinE = (unsigned int*)malloc(sizeof(unsigned int) * nBin * nRows);
			uBinC = (unsigned int*)malloc(sizeof(unsigned int) * nBin * nRows);
			for (int i = 0; i < nBin * nRows; i++) {
				uBinE[i] = 0;
				uBinC[i] = 0;
			}
			iMatrix = (int*)malloc(sizeof(int) * nCols * nRows);
			iShift = (int*)malloc(sizeof(int) * nRows);
			iGrow = (int*)malloc(sizeof(int) * nRows);
			iErow = (int*)malloc(sizeof(int) * nRows);
		}
	}

	for (int i = 0; i < nRows; i++) {
		fprintf(fo, "%4i | %2i %2i %3i |", i+1, iGrow[i]+1, iErow[i]+1, iShift[i]);
		for (int j = 0; j < nCols; j++) {
			if (iMatrix[i * nCols + j] == 0) {
				fprintf(fo, " %3i", iMatrix[i * nCols + j]);
			}
			else {
				fprintf(fo, " %3i", indExtra[iMatrix[i * nCols + j]-1]);
			}
		}
		fprintf(fo, "\n");
	}
	fprintf(fo, "\n");

	iCurCol = (int*)malloc(sizeof(int) * nCols);
	iCurLevel = (int*)malloc(sizeof(int) * nCols);
	nRinC = (int*)malloc(sizeof(int) * nCols);
	indexRow = (int*)malloc(sizeof(int) * nRows * nCols);
	

	for (int i = 0; i < nCols; i++) {
		n = 0;
		for (int j = 0; j < nRows; j++) {
			if (iMatrix[j * nCols + i] == 0)continue;
			indexRow[i * nRows + n] = j;
			n++;
		}
		nRinC[i] = n;
	}
	

	return 0;
}

int Seed2Simd::formVector() {
	int nStart;
	nLength = (int)strlen(mSeed);
	while (mSeed[nLength - 1] != '1') {
		mSeed[nLength - 1] = '\0';
		nLength--;
		if (nLength == 0) {
			printf("Error: no 1s in the seed\n");
			return -1;
		}
	}

	nStart = 0;
	while (mSeed[nStart] != '1') {
		nStart++;
		if (nStart == nLength) {
			printf("Error: no 1s in the seed\n");
			return -2;
		}
	}
	nLength -= nStart;
	for (int i = 0; i < nLength; i++) {
		mSeed[i] = mSeed[i + nStart];
	}
	mSeed[nLength] = '\0';

	fprintf(fo, "Length: %i\n", nLength);

	bVector = (bool*)malloc(sizeof(bool) * nLength);
	if (bVector == nullptr) {
		printf("Error: cannot allocate memory (bVector)\n");
		return -1;
	}
	
	nWeight = 0;

	for (int i = 0; i < nLength; i++) {
		if (mSeed[i] == '1') {
			bVector[i] = true;
			nWeight++;
		}else {
			bVector[i] = false;
		}
	}
	fprintf(fo, "Weight: %i\n\n", nWeight);

	nLength32 = nLength / 32;
	if (nLength % 32) nLength32++;
	nWeight32 = nWeight / 32;
	if (nWeight % 32) nWeight32++;

	prAllIn = (char*)malloc(sizeof(char) * 32 * nLength32);
	prAllOut = (char*)malloc(sizeof(char) * 32 * nLength32);

	printSeed();
	return 0;
}

int Seed2Simd::printSeed() {
	int ijk;
	fprintf(fo, "Corrected seed:\n");
	for (int i = 0; i < nLength; i++) {
		if (bVector[i]) {
			fprintf(fo, "1");
		}
		else {
			fprintf(fo, "0");
		}
	}
	
	fprintf(fo, "\n\n");

	for (int i = 0; i < nLength32; i++) {
		for (int j = 0; j < 4; j++) {
			for (int k = 0; k < 8; k++) {
				ijk = 32 * i + 8 * j + k;
				if (ijk >= nLength) {
					fprintf(fo, "0");
				}
				else {
					if (bVector[ijk]) {
						fprintf(fo, "1");
					}else{
						fprintf(fo, "0");
					}
				}
			}
			if (j < 3)fprintf(fo, " ");
		}
		fprintf(fo, "\n");
	}
	fprintf(fo, "\n");
	return 0;
}

int Seed2Simd::startProcessing() {
	fo = fopen(outputFile, "w");
	if (fo == nullptr) {
		printf("Error: cannot open the output file \"%s\".\n", outputFile);
		return -1;
	}

	fprintf(fo, "Input seed: \"%s\"\n\n", mSeed);
	
	if (formVector() != 0) return -2;
	if (formMatrix() != 0) return -3;
	if (findSolution() != 0) return -4;
	return 0;
}

int main(int argc, char* argv[]) {
	if (argc != 3) {
		printf("Parameters: 1) seed, 2) output file");
		return -1;
	}

	int ires;
	Seed2Simd *s2s;
	s2s = new Seed2Simd(argv[1], argv[2]);
	ires = s2s->startProcessing();
	delete s2s; s2s = nullptr;
	return ires;
}
