#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <Windows.h>
#include <math.h>


//비트맵 헤더를 한묶음으로
typedef struct tagBITMAPHEADER {
	BITMAPFILEHEADER bf;
	BITMAPINFOHEADER bi;
	RGBQUAD hRGB[256]; //이 코드에서는 필요없음 (8bit에만 필요)
}BITMAPHEADER;

//비트맵을 읽어와서 화소정보의 포인터를 리턴
BYTE* loadBitmapFile(int bytesPerPixel, BITMAPHEADER* bitmapHeader, int* imgWidth, int* imgHeight, char* filename);

//비트맵 파일 쓰기
void writeBitmapFile(int bytesPerPixel, BITMAPHEADER outputHeader, BYTE* output, int imgSize, char* filename);

// constructHaarMatrix
double** constructHaarMatrixRecursive(int n);
double** concatenateTwoMatrices(double** hl, double** hr, int n);
double** applyKroneckerProduct(double** A, int n, double a, double b);

void printMatrix(double** A, int m, int n, char name[]);
double** constructIdentity(int k);
double** allocateMemory(int m, int n);
void releaseMemory(double** A, int m);
double** multiplyTwoMatrices(double** A, int m, int n, double** B, int l, int k);
double** transposeMatrix(double** A, int m, int n);
double** NormailzeMatrix(double** A, int m, int n);
double** addTwoMatrices(double** A, int m, int n, double** B, int l, int k);

int CompareDoubleAbsoulte(double x, double y)
{
	double absTolerance = (1.0e-8);
	double diff = x - y;
	if (fabs(diff) <= absTolerance)
		return 1;
	else return 0;
}



int main() {
	/*******************************************************************/
	/*************************** Read image  ***************************/
	/*******************************************************************/
	BITMAPHEADER originalHeader;	//비트맵의 헤더부분을 파일에서 읽어 저장할 구조체
	BITMAPHEADER outputHeader;		//변형을 가한 헤더부분을 저장할 구조체
	int imgSize, imgWidth, imgHeight;					//이미지의 크기를 저장할 변수
	int bytesPerPixel = 3;			//number of bytes per pixel (1 byte for R,G,B respectively)

	BYTE* image = loadBitmapFile(bytesPerPixel, &originalHeader, &imgWidth, &imgHeight, "image_lena_24bit.bmp"); //비트맵파일을 읽어 화소정보를 저장 (불러들이는 이미지는 .c와 같은 폴더에 저장)
	if (image == NULL) return 0;

	imgSize = imgWidth * imgHeight; // total number of pixels
	BYTE* output = (BYTE*)malloc(bytesPerPixel * sizeof(BYTE) * imgSize);				//결과값을 저장할 포인터 선언 및 메모리 할당
	outputHeader = originalHeader;										//헤더정보를 출력헤더정보에 할당





	/*******************************************************************/
	/************************ Perform HWT/IHWT *************************/
	/*******************************************************************/
	//이미지 행렬 A 구성 (RGB값이 있으므로 픽셀당 값 하나씩만 읽어서 imgWidth x imgHeight 행렬 구성)
	double** A; //original image matrix
	double** H; // Haar matrix
	double** Ht;
	double** B;
	double** Bhat;
	double** Ahat;
	int n = imgHeight; //이미지가 정사각형(Height==Width)이라고 가정; n = 2^t,t=0,1,2,...
	A = (double**)malloc(sizeof(double*) * n); // byte
	for (int i = 0; i < n; i++) {
		A[i] = (double*)malloc(sizeof(double) * n);
	}

	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			A[i][j] = image[(i * n + j) * bytesPerPixel];



	//Haar matrix H 구성 (orthonormal column을 갖도록 구성)
	H = (double**)malloc(sizeof(double*) * n);
	H = constructHaarMatrixRecursive(n);
	H = NormailzeMatrix(H, n, n);

	//HWT 수행: 행렬곱 B = H'*A*H
	Ht = transposeMatrix(H, n, n);
	B = multiplyTwoMatrices(multiplyTwoMatrices(Ht, n, n, A, n, n),n, n, H, n, n);
	

	//행렬 B 자르기: B의 upper left corner(subsquare matrix)를 잘라 Bhat에 저장
	//Bhat은 B와 사이즈가 같으며 B에서 잘라서 저장한 부분 외에는 모두 0으로 채워짐
	Bhat = (double**)malloc(sizeof(double*) * n); // byte
	for (int i = 0; i < n; i++) {
		Bhat[i] = (double*)malloc(sizeof(double) * n);
	}

	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			Bhat[i][j] = 0.0;

	int m = 512; //이미지가 정사각형(Height==Width)이라고 가정; m = 2^t,t=0,1,2,...
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < m; j++) {
			Bhat[i][j] = B[i][j];
		}
	}
	

	//IHWT 수행: Ahat = H*Bhat*H'
	Ahat = multiplyTwoMatrices(multiplyTwoMatrices(H, n, n, Bhat, n, n), n, n, Ht, n, n);



	// 3.

	double** Hl = allocateMemory(n / 2, n);
	double** Hh = allocateMemory(n / 2, n);

	double** Hlt = transposeMatrix(Hl, n / 2, n);
	double** Hht = transposeMatrix(Hh, n / 2, n);

	double** ll = multiplyTwoMatrices(multiplyTwoMatrices(Hl, n / 2, n, A, n, n), n / 2, n, Hlt, n, n / 2);
	double** lh = multiplyTwoMatrices(multiplyTwoMatrices(Hl, n / 2, n, A, n, n), n / 2, n, Hht, n, n / 2);
	double** hl = multiplyTwoMatrices(multiplyTwoMatrices(Hh, n / 2, n, A, n, n), n / 2, n, Hlt, n, n / 2);
	double** hh = multiplyTwoMatrices(multiplyTwoMatrices(Hh, n / 2, n, A, n, n), n / 2, n, Hht, n, n / 2);

	double** HtAH = allocateMemory(n, n);
	double** HBHt = allocateMemory(n, n);

	double** llll = allocateMemory(n, n);
	double** llhh = allocateMemory(n, n);
	double** hhll = allocateMemory(n, n);
	double** hhhh = allocateMemory(n, n);

	double** Hll = allocateMemory(n / 4, n);
	double** Hlh = allocateMemory(n / 4, n);
	double** Hllt = allocateMemory(n / 4, n);
	double** Hlht = allocateMemory(n / 4, n);

	double** HltHlAHltHl = allocateMemory(n, n);

	double** llllllll = allocateMemory(n, n);
	double** lllllhlh = allocateMemory(n, n);
	double** lhlhllll = allocateMemory(n, n);
	double** lhlhlhlh = allocateMemory(n, n);

	for (int i = 0; i < n / 2; i++) {
		for (int j = 0; j < n; j++) {
			Hl[i][j] = Ht[i][j];
			Hh[i][j] = Ht[i + n / 2][j];
		}
	}

	for (int i = 0; i < n / 2; i++) {
		for (int j = 0; j < n / 2; j++) {
				HtAH[i][j] = ll[i][j];
		}
	}
	for (int i = 0; i < n / 2; i++) {
		for (int j = n/2; j < n; j++) {
			HtAH[i][j] = lh[i][j - n / 2];
		}
	}
	for (int i = n/2; i < n; i++) {
		for (int j = 0; j < n / 2; j++) {
			HtAH[i][j] = hl[i - n / 2][j];
		}
	}
	for (int i = n / 2; i < n; i++) {
		for (int j = n / 2; j < n; j++) {
			HtAH[i][j] = hh[i - n / 2][j - n / 2];
		}
	}

	// (a)
	char flag = 1;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (B[i][j] != HtAH[i][j]) {
				flag = 0;
				break;
			}
		}
	}
	if (flag) printf("B == HtAH\n");
	else printf("B != HtAH\n");

	A = multiplyTwoMatrices(multiplyTwoMatrices(H, n, n, HtAH, n, n), n, n, Ht, n, n);

	llll = multiplyTwoMatrices(multiplyTwoMatrices(multiplyTwoMatrices(multiplyTwoMatrices(Hlt, n, n / 2, Hl, n / 2, n), n, n, A, n, n), n, n, Hlt, n, n / 2), n, n / 2, Hl, n / 2, n);
	llhh = multiplyTwoMatrices(multiplyTwoMatrices(multiplyTwoMatrices(multiplyTwoMatrices(Hlt, n, n / 2, Hl, n / 2, n), n, n, A, n, n), n, n, Hht, n, n / 2), n, n / 2, Hh, n / 2, n);
	hhll = multiplyTwoMatrices(multiplyTwoMatrices(multiplyTwoMatrices(multiplyTwoMatrices(Hht, n, n / 2, Hh, n / 2, n), n, n, A, n, n), n, n, Hlt, n, n / 2), n, n / 2, Hl, n / 2, n);
	hhhh = multiplyTwoMatrices(multiplyTwoMatrices(multiplyTwoMatrices(multiplyTwoMatrices(Hht, n, n / 2, Hh, n / 2, n), n, n, A, n, n), n, n, Hht, n, n / 2), n, n / 2, Hh, n / 2, n);

	HBHt = addTwoMatrices(addTwoMatrices(addTwoMatrices(llll, n, n, llhh, n, n), n, n, hhll, n, n), n, n, hhhh, n, n);

	// (b)
	flag = 1;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (!CompareDoubleAbsoulte(A[i][j], HBHt[i][j])) {
				flag = 0;
				printf("%lf != %lf\n", A[i][j], HBHt[i][j]);
				break;
			}
		}
	}
	if (flag) printf("A == HBHt\n");
	else printf("A != HBHt\n");

	// (d)
	for (int i = 0; i < n / 4; i++) {
		for (int j = 0; j < n; j++) {
			Hll[i][j] = Ht[i][j];
			Hlh[i][j] = Ht[i + n / 4][j];
		}
	}

	Hllt = transposeMatrix(Hll, n / 4, n);
	Hlht = transposeMatrix(Hlh, n / 4, n);

	llllllll = multiplyTwoMatrices(multiplyTwoMatrices(multiplyTwoMatrices(multiplyTwoMatrices(Hllt, n, n / 4, Hll, n / 4, n), n, n, A, n, n), n, n, Hllt, n, n / 4), n, n / 4, Hll, n / 4, n);
	lllllhlh = multiplyTwoMatrices(multiplyTwoMatrices(multiplyTwoMatrices(multiplyTwoMatrices(Hllt, n, n / 4, Hll, n / 4, n), n, n, A, n, n), n, n, Hlht, n, n / 4), n, n / 4, Hlh, n / 4, n);
	lhlhllll = multiplyTwoMatrices(multiplyTwoMatrices(multiplyTwoMatrices(multiplyTwoMatrices(Hlht, n, n / 4, Hlh, n / 4, n), n, n, A, n, n), n, n, Hllt, n, n / 4), n, n / 4, Hll, n / 4, n);
	lhlhlhlh = multiplyTwoMatrices(multiplyTwoMatrices(multiplyTwoMatrices(multiplyTwoMatrices(Hlht, n, n / 4, Hlh, n / 4, n), n, n, A, n, n), n, n, Hlht, n, n / 4), n, n / 4, Hlh, n / 4, n);

	HltHlAHltHl = addTwoMatrices(addTwoMatrices(addTwoMatrices(llllllll, n, n, lllllhlh, n, n), n, n, lhlhllll, n, n), n, n, lhlhlhlh, n, n);
	flag = 1;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (!CompareDoubleAbsoulte(llll[i][j], HltHlAHltHl[i][j])) {
				flag = 0;
				printf("%lf != %lf\n", llll[i][j], HltHlAHltHl[i][j]);
				break;
			}
		}
	}
	if (flag) printf("llll == HltHlAHltHl\n");
	else printf("llll != HltHlAHltHl\n");

	
	/*******************************************************************/
	/******************* Write reconstructed image  ********************/
	/*******************************************************************/
	//Ahat을 이용해서 위의 image와 같은 형식이 되도록 구성 (즉, Ahat = [a b;c d]면 [a a a b b b c c c d d d]를 만들어야 함)
	BYTE* Are = (BYTE*)malloc(bytesPerPixel * sizeof(BYTE) * imgSize);
	int residx = 0;
	for (int i = 0; i < imgHeight; i++) {
		for (int j = 0; j < imgWidth; j++) {
			for (int k = 0; k < 3; k++) {
				Are[residx] = lhlhlhlh[i][j] ;
				residx++;
			}
		}
	}
	

	writeBitmapFile(bytesPerPixel, outputHeader, Are, imgSize, "lhlhlhlh.bmp");



	// free
	free(image);
	free(output);
	for (int i = 0; i < n; i++) {
		free(A[i]);
		free(H[i]);
		free(Ht[i]);
		free(B[i]);
		free(Bhat[i]);
		free(Ahat[i]);
	}
	free(A);
	free(H);
	free(Ht);
	free(B);
	free(Bhat);
	free(Ahat);
	free(Are);

	for (int i = 0; i < n / 2; i++) {
		free(Hl[i]);
		free(Hh[i]); 
		free(Hlt[i]); 
		free(Hht[i]); 
		free(ll[i]); 
		free(lh[i]); 
		free(hl[i]); 
		free(hh[i]);
	}		
	free(Hl);
	free(Hh);
	free(Hlt);
	free(Hht);
	free(ll);
	free(lh);	
	free(hl);
	free(hh);

	for (int i = 0; i < n; i++) {
		free(HtAH[i]); 
		free(HBHt[i]); 
		free(llll[i]); 
		free(llhh[i]);
		free(hhll[i]); 
		free(hhhh[i]);
	}
	free(HtAH);
	free(HBHt);	
	free(llll);	
	free(llhh);	
	free(hhll);	
	free(hhhh);

	for (int i = 0; i < n / 4; i++) {
		free(Hll[i]); 
		free(Hllt[i]); 
		free(Hlh[i]); 
		free(Hlht[i]);
	}	
	free(Hll);
	free(Hllt);
	free(Hlh);	
	free(Hlht);

	for (int i = 0; i < n; i++) {
		free(HltHlAHltHl[i]);
		free(llllllll[i]);
		free(lllllhlh[i]);
		free(lhlhllll[i]);
		free(lhlhlhlh[i]);
	}
	free(HltHlAHltHl);
	free(llllllll);
	free(lllllhlh);
	free(lhlhllll);
	free(lhlhlhlh);

	return 0;
}

BYTE* loadBitmapFile(int bytesPerPixel, BITMAPHEADER* bitmapHeader, int* imgWidth, int* imgHeight, char* filename)
{
	FILE* fp = fopen(filename, "rb");	//파일을 이진읽기모드로 열기
	if (fp == NULL)
	{
		printf("파일로딩에 실패했습니다.\n");	//fopen에 실패하면 NULL값을 리턴
		return NULL;
	}
	else
	{
		fread(&bitmapHeader->bf, sizeof(BITMAPFILEHEADER), 1, fp);	//비트맵파일헤더 읽기
		fread(&bitmapHeader->bi, sizeof(BITMAPINFOHEADER), 1, fp);	//비트맵인포헤더 읽기
		//fread(&bitmapHeader->hRGB, sizeof(RGBQUAD), 256, fp);	//색상팔렛트 읽기 (24bitmap 에서는 존재하지 않음)

		*imgWidth = bitmapHeader->bi.biWidth;
		*imgHeight = bitmapHeader->bi.biHeight;
		int imgSizeTemp = (*imgWidth) * (*imgHeight);	// 이미지 사이즈를 상위 변수에 할당

		printf("Size of image: Width %d   Height %d\n", bitmapHeader->bi.biWidth, bitmapHeader->bi.biHeight);
		BYTE* image = (BYTE*)malloc(bytesPerPixel * sizeof(BYTE) * imgSizeTemp);	//이미지크기만큼 메모리할당

		fread(image, bytesPerPixel * sizeof(BYTE), imgSizeTemp, fp);//이미지 크기만큼 파일에서 읽어오기

		fclose(fp);
		return image;
	}
}



void writeBitmapFile(int bytesPerPixel, BITMAPHEADER outputHeader, BYTE* output, int imgSize, char* filename)
{
	FILE* fp = fopen(filename, "wb");

	fwrite(&outputHeader.bf, sizeof(BITMAPFILEHEADER), 1, fp);
	fwrite(&outputHeader.bi, sizeof(BITMAPINFOHEADER), 1, fp);
	//fwrite(&outputHeader.hRGB, sizeof(RGBQUAD), 256, fp); //not needed for 24bitmap
	fwrite(output, bytesPerPixel * sizeof(BYTE), imgSize, fp);
	fclose(fp);
}

double** constructHaarMatrixRecursive(int n) {
	double** h;
	if (n > 2)
		h = constructHaarMatrixRecursive(n / 2);
	else {
		//double** h;
		h = allocateMemory(2, 2);
		h[0][0] = 1; h[0][1] = 1; h[1][0] = 1; h[1][1] = -1; //H = [1 1; 1 -1]
		return h;
	}

	// construct the left half (Kronecket product of h and [1;1])
	double** Hl = applyKroneckerProduct(h, n, 1, 1);
	releaseMemory(h, n / 2);

	// construct the right half (Kronecker product of I and [1;-1])
	double** I = constructIdentity(n / 2);
	double** Hr = applyKroneckerProduct(I, n, 1, -1);
	releaseMemory(I, n / 2);


	// merge hl and hr
	double** H = concatenateTwoMatrices(Hl, Hr, n); //H = [Hl Hr]
	releaseMemory(Hl, n);
	releaseMemory(Hr, n);

	return H;
}

double** applyKroneckerProduct(double** A, int n, double a, double b) {
	double** h = allocateMemory(n, n / 2);

	for (int j = 0; j < n / 2; j++) {
		for (int i = 0; i < n / 2; i++) {
			h[2 * i][j] = A[i][j] * a;
			h[2 * i + 1][j] = A[i][j] * b;
		}
	}
	return h;
}

double** concatenateTwoMatrices(double** hl, double** hr, int n) {
	double** H = allocateMemory(n, n);

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (j < n / 2)
				H[i][j] = hl[i][j];
			else
				H[i][j] = hr[i][j - n / 2];
		}
	}
	return H;
}


void printMatrix(double** A, int m, int n, const char name[]) {
	printf("\n%s = \n", name);
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++)
			printf("%lf ", A[i][j]);
		printf("\n");
	}
}

double** constructIdentity(int k) {
	double** I = allocateMemory(k, k);

	for (int i = 0; i < k; i++) {
		for (int j = 0; j < k; j++) {
			if (i != j)
				I[i][j] = 0.0;
			else
				I[i][j] = 1.0;
		}
	}
	return I;
}

double** allocateMemory(int m, int n) {
	double** A;
	A = (double**)malloc(sizeof(double*) * m);
	for (int i = 0; i < m; i++) {
		A[i] = (double*)malloc(sizeof(double) * n);
	}
	return A;
}


void releaseMemory(double** A, int m) {
	for (int i = 0; i < m; i++)
		free(A[i]);
	free(A);
}

double** multiplyTwoMatrices(double** A, int m, int n, double** B, int l, int k) {
	if (n != l) {
		printf("NULL");
		return NULL;
	}

	double** res = allocateMemory(m, k);

	for (int i = 0; i < m; i++) {
		for (int j = 0; j < k; j++) {
			res[i][j] = 0.0;
			for (int x = 0; x < n; x++) {
				res[i][j] += A[i][x] * B[x][j];
			}
		}
	}

	return res;
}

double** transposeMatrix(double** A, int m, int n)
{
	double** B = allocateMemory(n, m);

	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			B[j][i] = A[i][j];

	return B;
}

double** NormailzeMatrix(double** A, int m, int n)
{
	double** B = allocateMemory(m, n);

	for (int i = 0; i < n; i++)
	{
		double len = 0.0;

		for (int j = 0; j < m; j++)
			len += A[j][i] * A[j][i];
		len = sqrt(len);

		for (int j = 0; j < m; j++)
			B[j][i] = A[j][i] / len;
	}

	return B;
}

double** addTwoMatrices(double** A, int m, int n, double** B, int l, int k) {
	if (m != l || n != k) {
		printf("NULL");
		return NULL;
	}

	double** result = allocateMemory(m, n);

	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			result[i][j] = A[i][j] + B[i][j];
		}
	}

	return result;
}
