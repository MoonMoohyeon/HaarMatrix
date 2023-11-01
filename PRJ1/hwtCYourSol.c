#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <Windows.h>
#include <math.h>


//��Ʈ�� ����� �ѹ�������
typedef struct tagBITMAPHEADER {
	BITMAPFILEHEADER bf;
	BITMAPINFOHEADER bi;
	RGBQUAD hRGB[256]; //�� �ڵ忡���� �ʿ���� (8bit���� �ʿ�)
}BITMAPHEADER;

//��Ʈ���� �о�ͼ� ȭ�������� �����͸� ����
BYTE* loadBitmapFile(int bytesPerPixel, BITMAPHEADER* bitmapHeader, int* imgWidth, int* imgHeight, char* filename);

//��Ʈ�� ���� ����
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

int main() {
	/*******************************************************************/
	/*************************** Read image  ***************************/
	/*******************************************************************/
	BITMAPHEADER originalHeader;	//��Ʈ���� ����κ��� ���Ͽ��� �о� ������ ����ü
	BITMAPHEADER outputHeader;		//������ ���� ����κ��� ������ ����ü
	int imgSize, imgWidth, imgHeight;					//�̹����� ũ�⸦ ������ ����
	int bytesPerPixel = 3;			//number of bytes per pixel (1 byte for R,G,B respectively)

	BYTE* image = loadBitmapFile(bytesPerPixel, &originalHeader, &imgWidth, &imgHeight, "image_lena_24bit.bmp"); //��Ʈ�������� �о� ȭ�������� ���� (�ҷ����̴� �̹����� .c�� ���� ������ ����)
	if (image == NULL) return 0;

	imgSize = imgWidth * imgHeight; // total number of pixels
	BYTE* output = (BYTE*)malloc(bytesPerPixel * sizeof(BYTE) * imgSize);				//������� ������ ������ ���� �� �޸� �Ҵ�
	outputHeader = originalHeader;										//��������� ������������ �Ҵ�





	/*******************************************************************/
	/************************ Perform HWT/IHWT *************************/
	/*******************************************************************/
	//�̹��� ��� A ���� (RGB���� �����Ƿ� �ȼ��� �� �ϳ����� �о imgWidth x imgHeight ��� ����)
	double** A; //original image matrix
	double** H; // Haar matrix
	double** Ht;
	A = (double**)malloc(sizeof(double*) * imgHeight); // byte
	for (int i = 0; i < imgHeight; i++) {
		A[i] = (double*)malloc(sizeof(double) * imgWidth);
	}

	for (int i = 0; i < imgHeight; i++)
		for (int j = 0; j < imgWidth; j++)
			A[i][j] = image[(i * imgWidth + j) * bytesPerPixel];



	//Haar matrix H ���� (orthonormal column�� ������ ����)
	int n = imgHeight; //�̹����� ���簢��(Height==Width)�̶�� ����; n = 2^t,t=0,1,2,...
	H = (double**)malloc(sizeof(double*) * n);
	H = constructHaarMatrixRecursive(n);
	/*for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			printf("%lf ", H[i][j]);
		}
		printf("\n");
	}*/
	H = NormailzeMatrix(H, n, n);
	/*for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			printf("%lf ", H[i][j]);
		}
		printf("\n");
	}*/

	//HWT ����: ��İ� B = H'*A*H
	Ht = transposeMatrix(H, n, n);
	double** B = multiplyTwoMatrices(multiplyTwoMatrices(Ht, n, n, A, n, n),n, n, H, n, n);
	/*for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			printf("%lf ", B[i][j]);
		}
		printf("\n");
	}*/
	

	//��� B �ڸ���: B�� upper left corner(subsquare matrix)�� �߶� Bhat�� ����
	//Bhat�� B�� ����� ������ B���� �߶� ������ �κ� �ܿ��� ��� 0���� ä����
	double** Bhat = (double**)malloc(sizeof(double*) * imgHeight); // byte
	for (int i = 0; i < imgHeight; i++) {
		Bhat[i] = (double*)malloc(sizeof(double) * imgWidth);
	}

	for (int i = 0; i < imgHeight; i++)
		for (int j = 0; j < imgWidth; j++)
			Bhat[i][j] = 0.0;

	int m = 64; //�̹����� ���簢��(Height==Width)�̶�� ����; m = 2^t,t=0,1,2,...
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < m; j++) {
			Bhat[i][j] = B[i][j];
		}
	}
	

	//IHWT ����: Ahat = H*Bhat*H'
	double** Ahat = multiplyTwoMatrices(multiplyTwoMatrices(H, n, n, Bhat, n, n), n, n, Ht, n, n);


	/*******************************************************************/
	/******************* Write reconstructed image  ********************/
	/*******************************************************************/
	//Ahat�� �̿��ؼ� ���� image�� ���� ������ �ǵ��� ���� (��, Ahat = [a b;c d]�� [a a a b b b c c c d d d]�� ������ ��)
	BYTE* Are = (BYTE*)malloc(bytesPerPixel * sizeof(BYTE) * imgSize);
	int residx = 0;
	for (int i = 0; i < imgHeight; i++) {
		for (int j = 0; j < imgWidth; j++) {
			for (int k = 0; k < 3; k++) {
				Are[residx] = Ahat[i][j];
				residx++;
			}
		}
	}
		
	

	writeBitmapFile(bytesPerPixel, outputHeader, Are, imgSize, "output1.bmp");


	free(image);
	free(output);
	for (int i = 0; i < imgHeight; i++)
		free(A[i]);
	free(A);
	free(Are);

	return 0;
}

BYTE* loadBitmapFile(int bytesPerPixel, BITMAPHEADER* bitmapHeader, int* imgWidth, int* imgHeight, char* filename)
{
	FILE* fp = fopen(filename, "rb");	//������ �����б���� ����
	if (fp == NULL)
	{
		printf("���Ϸε��� �����߽��ϴ�.\n");	//fopen�� �����ϸ� NULL���� ����
		return NULL;
	}
	else
	{
		fread(&bitmapHeader->bf, sizeof(BITMAPFILEHEADER), 1, fp);	//��Ʈ��������� �б�
		fread(&bitmapHeader->bi, sizeof(BITMAPINFOHEADER), 1, fp);	//��Ʈ��������� �б�
		//fread(&bitmapHeader->hRGB, sizeof(RGBQUAD), 256, fp);	//�����ȷ�Ʈ �б� (24bitmap ������ �������� ����)

		*imgWidth = bitmapHeader->bi.biWidth;
		*imgHeight = bitmapHeader->bi.biHeight;
		int imgSizeTemp = (*imgWidth) * (*imgHeight);	// �̹��� ����� ���� ������ �Ҵ�

		printf("Size of image: Width %d   Height %d\n", bitmapHeader->bi.biWidth, bitmapHeader->bi.biHeight);
		BYTE* image = (BYTE*)malloc(bytesPerPixel * sizeof(BYTE) * imgSizeTemp);	//�̹���ũ�⸸ŭ �޸��Ҵ�

		fread(image, bytesPerPixel * sizeof(BYTE), imgSizeTemp, fp);//�̹��� ũ�⸸ŭ ���Ͽ��� �о����

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