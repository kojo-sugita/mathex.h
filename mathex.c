
/** 
 * @file    mathex.c
 * @brief   math.h�̊g��
 *
 * @author  Kojo Sugita
 * @date    2009-08-05
 */

/* header files */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "mathex.h"

/**
 * ���a���v�Z����
 * @param[in] a ���a���v�Z����z��
 * @param[in] n �z��̗v�f��
 * @return ���a
 */
double Summation(const double *a, size_t n) {
	double sum = 0.0;

	while (n--) 
		sum += a[n];

	return sum;
}

/**
 * �����a���v�Z����
 * @param[in] a �����a���v�Z����z��
 * @param[in] n �z��̗v�f��
 * @return �����a
 */
double SumOfSquares(const double *a, size_t n) {
	double sum = 0.0;

	while (n--)
		sum += (a[n] * a[n]);

	return sum;
}

/**
 * �Ϙa���v�Z����
 * @param[in] a1 �Ϙa���v�Z����z��1
 * @param[in] a2 �Ϙa���v�Z����z��2
 * @param[in] n �z��̗v�f��
 * @return �Ϙa
 */
double SumOfProduct(const double *a1, const double *a2, size_t n) {
	double sum = 0.0;

	while (n--)
		sum += (a1[n] * a2[n]);

	return sum;
}

/**
 * ���ϒl���v�Z����
 * @param[in] a ���͔z��
 * @param[in] n �z��̗v�f��
 * @return ���ϒl
 */
double Average(const double *a, size_t n) {
	size_t i;
	double sum = 0.0;

	if ( n < 1 ) return 0.0;

	for (i = 0; i < n; i++) {
		sum += a[i];
	}
	return sum / (double)n;
}

/**
 * ���U���v�Z����
 * @param[in] a ���͔z��
 * @param[in] n �z��̗v�f��
 * @return ���U
 */
double Variance(const double *a, size_t n) {
	size_t i;
	double sum = 0.0;
	double sqsum = 0.0;
	double avg;

	if ( n < 1 ) return 0.0;

	for (i = 0; i < n; i++) {
		sum += a[i]; /* �a */
		sqsum += a[i] * a[i]; /* �����a */
	}

	avg = sum / (double)n; /* ���� */
	return sqsum / (double)n - avg * avg; /* ���U */
}

/**
 * �s�Ε��U���v�Z����
 * @param[in] a ���͔z��
 * @param[in] n �z��̗v�f��
 * @return �s�Ε��U
 */
double UnbiasedVariance(const double *a, size_t n) {
    size_t i;
    double sum = 0.0;
    double sqsum = 0.0;
    double avg;

    if ( n < 2 ) return 0.0;

    for (i = 0; i < n; i++) {
        sum += a[i]; /* �a */
        sqsum += a[i] * a[i]; /* �����a */
    }

    avg = sum / (double)n; /* ���� */
    return sqsum / (double)(n - 1) - avg * avg; /* �s�Ε��U */
}

/**
 * �����U���v�Z����
 * @param[in] x ���͔z��1
 * @param[in] y ���͔z��2
 * @param[in] n x��y�̗v�f��
 * @return �����U
 */
double Covariance(const double *x, const double *y, int n) {
    int i;
    double sum_x = 0.0, sum_y = 0.0;
    double avg_x, avg_y;
    double sumprod = 0.0;

    if ( n < 1 ) return 0.0;

    for (i = 0; i < n; i++) {
        sum_x += x[i]; /* x�̘a */
        sum_y += y[i]; /* y�̘a */
        sumprod += x[i] * y[i];
    }

    avg_x = sum_x / (double)n; /* x�̕��� */
    avg_y = sum_y / (double)n; /* y�̕��� */
    return sumprod / (double)n - avg_x * avg_y; /* �����U */
}

/**
 * �W���΍����v�Z����
 * @param[in] a ���͔z��
 * @param[in] n �z��̗v�f��
 * @return �W���΍�
 */
double StdDev(const double *a, size_t n) {
	size_t i;
	double sum = 0.0;
	double sqsum = 0.0;
	double avg, var;

	if ( n < 1 ) return 0.0;

	for (i = 0; i < n; i++) {
		sum += a[i]; /* �a */
		sqsum += a[i] * a[i]; /* �����a */
	}

	avg = sum / (double)n; /* ���� */
	var = sqsum / (double)n - avg * avg; /* ���U */

	if ( var < 0.0 ) var = 0.0;
	return sqrt(var); /* �W���΍� */

}

/**
 * �ő�l�����߂�
 * @param[in] a ���͔z��
 * @param[in] n �z��̗v�f��
 * @return �ő�l
 */
double Max(const double *a, size_t n) {
	size_t i;
	double max;

	if ( n < 1 ) return 0.0;

	max = a[0];
	for (i = 1; i < n; i++) {
		if ( max < a[i] ) max = a[i];
	}

	return max;
}

/**
 * �ő�l�����߂�(int�^�z��)
 * @param[in] a ���͔z��
 * @param[in] n �z��̗v�f��
 * @return �ő�l
 */
int MaxInt(const int *a, size_t n) {
	size_t i;
	int max;

	if ( n < 1 ) return 0;

	max = a[0];
	for (i = 1; i < n; i++) {
		if ( max < a[i] ) max = a[i];
	}

	return max;
}


/**
 * �ŏ��l�����߂�
 * @param[in] a ���͔z��
 * @param[in] n �z��̗v�f��
 * @return �ŏ��l
 */
double Min(const double *a, size_t n) {
	size_t i;
	double min;

	if ( n < 1 ) return 0.0;

	min = a[0];
	for (i = 1; i < n; i++) {
		if ( min > a[i] ) min = a[i];
	}

	return min;
}

/**
 * �ŏ��l�����߂�(int�^�z��)
 * @param[in] a ���͔z��
 * @param[in] n �z��̗v�f��
 * @return �ŏ��l
 */
int MinInt(const int *a, size_t n) {
	size_t i;
	int min;

	if ( n < 1 ) return 0;

	min = a[0];
	for (i = 1; i < n; i++) {
		if ( min > a[i] ) min = a[i];
	}

	return min;
}

/**
 * �ő�l�ƍŏ��l�����߂�
 * @param[in] a ���͔z��
 * @param[in] n �z��̗v�f��
 * @param[out] max �ő�l
 * @param[out] min �ŏ��l
 * @retval 1 �ő�l�ƍŏ��l���擾
 * @retval 0 �z��̗v�f��0
 */
int MaxMin(const double *a, size_t n, double *max, double *min) {
	size_t i;

	if ( n < 1 ) return 0;

	*max = a[0];
	*min = a[0];
	for (i = 1; i < n; i++) {
		if ( *max < a[i] ) *max = a[i];
		if ( *min > a[i] ) *min = a[i];
	}

	return 1;
}


/**
 * �K����v�Z����
 * @param[in] n ���� 
 * @return �K��
 */
int Factorial(int n) {
	int fact = n;

	if (n <= 0) return 1;

	while (--n) 
		fact *= n;

	return fact;
}

/**
 * ����g�ݍ��킹�̑������v�Z����
 * @param[in] n nPr��n 
 * @param[in] r nPr��r 
 * @return ����g�ݍ��킹�̑���
 */
int Permutation(int n, int r) {
	if (r <= 0) return 1;
	/* 0���Z�͂��肦�Ȃ� */
	return Factorial(n) / Factorial(n - r);
}

/**
 * �g�ݍ��킹�̑������v�Z����
 * @param[in] n nCr��n 
 * @param[in] r nCr��r 
 * @return �g�ݍ��킹�̑���
 */
int Combination(int n, int r) {
	if (r <= 0) return 1;
	/* 0���Z�͂��肦�Ȃ� */
	return Factorial(n) / (Factorial(r) * Factorial(n - r));
}


/**
 * 2�_�Ԃ̋��������߂�
 * @param[in] p1 �_�̍��W1
 * @param[in] p2 �_�̍��W2
 * @param[in] n ������
 * @return p1��p2�Ԃ̋���
 */
double Distance(const double *p1, const double *p2, size_t n) {
	size_t i;
	double distance = 0.0;

	for (i = 0; i < n; i++) {
		distance += (p2[i] - p1[i]) * (p2[i] - p1[i]);
	}
	return sqrt(distance);
}

/**
 * 2�_p1,p2����Ȃ������m:n�ɓ�������_�����߂� (2����)
 * @param[in] p1 �_�̍��W1
 * @param[in] p2 �_�̍��W2
 * @param[in] m m:n��m
 * @param[in] n m:n��n
 * @param[out] result �����_
 * @retval 1 ����I��
 * @retval 0 �ُ�I��
 */
int DivideInternally2D(const double p1[2], const double p2[2], int m, int n, double result[2]) {
	if (m < 1 && n < 1) return 0;
	
	result[0] = (n * p1[0] + m * p2[0]) / (m + n);
	result[1] = (n * p1[1] + m * p2[1]) / (m + n);

	return 1;
}

/**
 * 2�_p1,p2����Ȃ������m:n�ɊO������_�����߂� (2����)
 * @param[in] p1 �_�̍��W1
 * @param[in] p2 �_�̍��W2
 * @param[in] m m:n��m
 * @param[in] n m:n��n
 * @param[out] result �O���_
 * @retval 1 ����I��
 * @retval 0 �ُ�I��
 */
int DivideExternally2D(const double p1[2], const double p2[2], int m, int n, double result[2]) {
	if (m - n == 0) return 0;

	result[0] = (-n * p1[0] + m * p2[0]) / (m - n);
	result[1] = (-n * p1[1] + m * p2[1]) / (m - n);

	return 1;
}

/**
 * 3�_p1,p2,p3����Ȃ�O�p�`�̏d�S�����߂� (2����)
 * @param[in] p1 �_�̍��W1
 * @param[in] p2 �_�̍��W2
 * @param[in] p3 �_�̍��W3
 * @param[out] result �O�p�`�̏d�S
 * @return result�̐擪�A�h���X
 */
double *Triangle_CenterOfGravity(const double p1[2], const double p2[2], const double p3[2], double result[2]) {
	result[0] = (p1[0] + p2[0] + p3[0]) / 3.0;
	result[1] = (p1[1] + p2[1] + p3[1]) / 3.0;
	return result;
}

/**
 * 2�_�����ԃx�N�g�������߂�
 * @param[in] p1 �_�̍��W1
 * @param[in] p2 �_�̍��W2
 * @param[in] n ������
 * @param[out] vec �x�N�g��
 * @return vec �̐擪�A�h���X
 */
double *ToVector(const double *p1, const double *p2, size_t n, double *vec) {
	size_t i;
	
	for ( i = 0; i < n; i++ ) {
		vec[i] = p2[i] - p1[i];
	}
	return vec;
}

/**
 * �x�N�g���̘a�����߂�
 * @param[in] vec1 �x�N�g��1
 * @param[in] vec2 �x�N�g��2
 * @param[out] sumVec �x�N�g���̘a
 * @param[in] n ������
 * @return sumVec�̐擪�A�h���X
 */
double *VectorSum(const double *vec1, const double *vec2, double *sumVec, size_t n) {
	while (n--) 
		sumVec[n] = vec1[n] + vec2[n];
	return sumVec;
}

/**
 * �x�N�g���̍������߂�
 * @param[in] vec1 �x�N�g��1
 * @param[in] vec2 �x�N�g��2
 * @param[out] sumVec �x�N�g���̘a
 * @param[in] n ������
 * @return diffVec�̐擪�A�h���X
 */
double *VectorDifference(const double *vec1, const double *vec2, double *diffVec, size_t n) {
	while (n--) 
		diffVec[n] = vec1[n] - vec2[n];
	return diffVec;
}

/**
 * �x�N�g���̒���(�m����)���v�Z����
 * @param[in] vec �x�N�g��
 * @param[in] n �x�N�g���̎�����
 * @return vec �̒���
 */
double Norm(const double *vec, size_t n) {
    size_t i;
    double s = 0.0;

    for ( i = 0; i < n; i++ ) {
        s += vec[i] * vec[i];
    }

    return sqrt(s);
}

/**
 * �x�N�g���̓��ς��v�Z����
 * @param[in] vec1 �x�N�g��1
 * @param[in] vec2 �x�N�g��2
 * @param[in] n �x�N�g���̎�����
 * @return vec1 �� vec2 �̓���
 */
double InnerProduct(const double *vec1, const double *vec2, size_t n) {
    size_t i;
    double s = 0.0;

    for ( i = 0; i < n; i++ ) {
        s += vec1[i] * vec2[i];
    }

    return s;
}

/**
 * 2 �̃x�N�g���̂Ȃ��p���v�Z����
 * @param[in] vec1 �x�N�g��1
 * @param[in] vec2 �x�N�g��2
 * @param[in] n �x�N�g���̎�����
 * @return vec1 �� vec2 �̂Ȃ��p
 */
double IncludedAngle(const double *vec1, const double *vec2, size_t n) {
    return acos(InnerProduct(vec1, vec2, n) / 
        (Norm(vec1, n) * Norm(vec2, n) + DBL_MIN));
}

/**
 * �P�ʃx�N�g�������߂�
 * @param[in] vec �x�N�g��
 * @param[out] unitVec �P�ʃx�N�g��
 * @param[in] n �x�N�g���̎�����
 * @retval 1 �P�ʃx�N�g���𐳏�Ɏ擾
 * @retval 0 �x�N�g�����[���x�N�g��
 */
int UnitVector(const double *vec, double *unitVec, size_t n) {
	size_t i;
	double norm = Norm(vec, n);

	/* �x�N�g���̃m������0�łȂ���? */
	if (norm > 0.0) { 
		/* �P�ʃx�N�g�����擾 */
		for (i = 0; i < n; i++) {
			unitVec[i] = vec[i] / Norm(vec, n);
		}
		return 1;

	} else {
		/* ���ʂɂ��ׂ�0���i�[���� */
		for (i = 0; i < n; i++) {
			unitVec[i] = 0;
		}
		return 0;
	}
}

/**
 * 3�����̖@���x�N�g�������߂�
 * @param[in] p1 ���W1
 * @param[in] p2 ���W2
 * @param[in] p3 ���W3
 * @param[in] result �@���x�N�g��
 * @return result�̐擪�A�h���X
 */
double *NormalVector3D(const double *p1, const double *p2, const double *p3, double *result) {
	double vec1[3];
	double vec2[3];

	ToVector(p2, p1, 3, vec1);
	ToVector(p2, p3, 3, vec2);

	CrossProduct3D(vec1, vec2, result);

	return result;
}

/**
 * 3�����x�N�g���̊O�ς��v�Z����
 * @param[in] vec1 �x�N�g��1
 * @param[in] vec2 �x�N�g��2
 * @param[out] result �x�N�g���̊O��
 * @return result�̐擪�A�h���X
 */
double *CrossProduct3D(const double *vec1, const double *vec2, double *result) {
	result[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1];
	result[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2];
	result[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];

	return result;
}

/**
 * 2 �̃x�N�g���̗ގ��x���v�Z����
 * @param[in] vec1 �x�N�g��1
 * @param[in] vec2 �x�N�g��2
 * @param[in] n �x�N�g���̎�����
 * @return vec1 �� vec2 �̗ގ��x (-1.0 <= return <= 1.0) 1.0�ɋ߂��قǗގ��x��
 * @attention �[���x�N�g���������ɓn���ꂽ�ꍇ0��Ԃ�
 */
double VectorSimilarity(const double *vec1, const double *vec2, size_t n) {
    return InnerProduct(vec1, vec2, n) / 
		((Norm(vec1, n) * Norm(vec2, n)) + DBL_MIN);
}

/**
 * ���K���������v�Z����
 * @param[in] vec1 �x�N�g��1
 * @param[in] vec2 �x�N�g��2
 * @param[in] n �x�N�g���̎�����
 * @return vec1 �� vec2 �̐��K������
 */
double NormalizedVectorDistance(const double *vec1, const double *vec2, size_t n) {
	size_t i;
	double unitDiff;
	double sum = 0.0;

	for (i = 0; i < n; i++) {
		unitDiff = (vec1[i] / (Norm(vec1, n) + DBL_MIN)) - (vec2[i] / (Norm(vec2, n) + DBL_MIN));
		sum += unitDiff * unitDiff;
	}

	return sqrt(sum);
}

/**
 * 3������ԏ�̓_p1�Ɠ_p2�����Ԓ�����ɓ_p�����݂��邩�𔻒肷��
 * @param[in] p1 �_1
 * @param[in] p2 �_2
 * @param[in] p 3������ԏ�̂���_
 * @retval 1 p1��p2�����Ԓ������p�����݂���
 * @rarval 0 p1��p2�����Ԓ������p�����݂��Ȃ�
 */
int IsDepthOnLine3D(const double *p1, const double *p2, const double *p) {
	double x = p[0];
	double y = p[1];
	double z = p[2];

	double v1 = (x - p1[0]) / (p2[0] - p1[0] + (1.0e-8));
	double v2 = (y - p1[1]) / (p2[1] - p1[1] + (1.0e-8));
	double v3 = (z - p1[2]) / (p2[2] - p1[2] + (1.0e-8));

	if ( (fabs(v1 - v2) <= 1.0e-8) && (fabs(v1 - v3) <= 1.0e-8)) {
		return 1;
	} else {
		return 0;
	}
}

/**
 * 3�������ʂ��쐬����
 * @param[in] AB 3�����x�N�g��1
 * @param[in] AC 3�����x�N�g��2
 * @param[in] A ���W1
 * @param[out] a ���ʂ̕����� ax + by + cz + d = 0�̌W���� a
 * @param[out] b ���ʂ̕����� ax + by + cz + d = 0 �̌W���� b 
 * @param[out] c ���ʂ̕����� ax + by + cz + d = 0 �̌W���� c
 * @param[out] d ���ʂ̕����� ax + by + cz + d = 0 �̌W���� d
 */
void CalculatePlane3D(const double *AB, const double *AC, const double *A, double *a, double *b, double *c, double *d) {
	double norm;

	/* �@���x�N�g�� */
	*a = AB[1] * AC[2] - AB[2] * AC[1];
	*b = AB[2] * AC[0] - AB[0] * AC[2];
	*c = AB[0] * AC[1] - AB[1] * AC[0];

	norm = sqrt(*a * *a + *b * *b + *c * *c);
	*a = *a / (norm + DBL_MIN);
	*b = *b / (norm + DBL_MIN);
	*c = *c / (norm + DBL_MIN);

	/* ���_���畽�ʂɐL�΂��������̒��� */
	*d = -1.0 * *a * A[0] + -1.0 * *b * A[1] + -1.0 * *c * A[2];
}

/**
 * 3�_����3�������ʂ��쐬����
 * @param[in] P1 �_1
 * @param[in] P2 �_2
 * @param[in] P3 �_3
 * @param[out] a ���ʂ̕����� ax + by + cz + d = 0�̌W���� a
 * @param[out] b ���ʂ̕����� ax + by + cz + d = 0 �̌W���� b 
 * @param[out] c ���ʂ̕����� ax + by + cz + d = 0 �̌W���� c
 * @param[out] d ���ʂ̕����� ax + by + cz + d = 0 �̌W���� d
 */
void CalculatePlane3D_Using3Depths(double P1[3], double P2[3], double P3[3], double *a, double *b, double *c, double *d) {
	double AB[3], AC[3];
	double norm;

	/* �x�N�g���ɕϊ� */
	ToVector(P1, P2, 3, AB);
	ToVector(P1, P3, 3, AC);

	/* �@���x�N�g�� */
	*a = AB[1] * AC[2] - AB[2] * AC[1];
	*b = AB[2] * AC[0] - AB[0] * AC[2];
	*c = AB[0] * AC[1] - AB[1] * AC[0];

	norm = sqrt(*a * *a + *b * *b + *c * *c);
	*a = *a / (norm + DBL_MIN);
	*b = *b / (norm + DBL_MIN);
	*c = *c / (norm + DBL_MIN);

	/* ���_���畽�ʂɐL�΂��������̒��� */
	*d = -1.0 * *a * P1[0] + -1.0 * *b * P1[1] + -1.0 * *c * P1[2];
}

/**
 * �����̕�������W���`�����ʌ`�ɕό`���� (y = ax + b --> ax + by + c = 0)
 * @param[in] a �W���`�̌W�� (�X��)
 * @param[in] b �W���`�̌W�� (�ؕ�)
 * @param[out] ga ��ʌ`�̌W����a
 * @param[out] gb ��ʌ`�̌W����b
 * @param[out] gc ��ʌ`�̌W����c
 * @return EXIT_SUCCESS
 */
int ToLineCanonicalForm(double a, double b, double *ga, double *gb, double *gc)
{
	double norm;

	/* �@���x�N�g�� */
	norm = sqrt(a * a + 1.0);
	*ga = a / (norm + DBL_MIN);
	*gb = -1.0 / (norm + DBL_MIN);

	/* ���_���畽�ʂɐL�΂��������̒��� */
	*gc = b;
	return EXIT_SUCCESS;
}


/**
 * ���ʂ̕�������W���`�����ʌ`�ɕό`���� (z = ax + by + c --> ax + by + cz + d = 0)
 * @param[in] a �W���`�̌W��a
 * @param[in] b �W���`�̌W��b
 * @param[in] c �W���`�̌W��c
 * @param[out] ga ��ʌ`�̌W����a
 * @param[out] gb ��ʌ`�̌W����b
 * @param[out] gc ��ʌ`�̌W����c
 * @param[out] gd ��ʌ`�̌W����d
 * @return EXIT_SUCCESS
 */
int ToPlaneCanonicalForm(double a, double b, double c, double *ga, double *gb, double *gc, double *gd) {

	/* �@���x�N�g�� */
	//*ga = a;  /* �W��a */
	//*gb = b;  /* �W��b */
	//*gc = -1; /* �W��c */

	double norm;

	norm = sqrt(a * a + b * b + 1.0);
	*ga = a / (norm + DBL_MIN);
	*gb = b / (norm + DBL_MIN);
	*gc = -1.0 / (norm + DBL_MIN);

	/* ���_���畽�ʂɐL�΂��������̒��� */
	*gd = c;  /* �W��d */

	return EXIT_SUCCESS;
}

/**
 * ���ʂ̕���������ʌ`����W���`�ɕό`���� (z = ax + by + c --> ax + by + cz + d = 0)
 * @param[in] a ��ʌ`�̌W����a
 * @param[in] b ��ʌ`�̌W����b
 * @param[in] c ��ʌ`�̌W����c
 * @param[in] d ��ʌ`�̌W����d
 * @param[out] ca �W���`�̌W��a
 * @param[out] cb �W���`�̌W��b
 * @param[out] cc �W���`�̌W��c
 * @retval EXIT_SUCCESS �ϊ�����
 * @retval EXIT_FAILURE �ϊ����s(0���Z)
 */
int ToPlaneGenericForm(double a, double b, double c, double d, double *ca, double *cb, double *cc) {
	if (c != 0) {
		/* �X�� */
		*ca = - a / c; /* �W��a */
		*cb = - b / c; /* �W��b */

		/* �ؕ� */
		*cc = - d / c; /* �W��c */
		return EXIT_SUCCESS;

	} else {
		*ca = 0;
		*cb = 0;
		*cc = 0;
		return EXIT_FAILURE;
	}
}

/**
 * �_�ƒ����Ƃ̋�����Ԃ�
 * @param[in] a ��ʌ`�̌W����a
 * @param[in] b ��ʌ`�̌W����b
 * @param[in] c ��ʌ`�̌W����c
 * @param[in] point �_
 * @return �_�ƒ����Ƃ̋���
 */
double DistanceLineToDepth(double a, double b, double c, double point[2]) {
	double x = point[0];
	double y = point[1];

	return fabs(a * x + b * y + c) / (sqrt(a * a + b * b) + DBL_MIN);
}

/**
 * �_�ƕ��ʂƂ̋�����Ԃ�
 * @param[in] a ��ʌ`�̌W����a
 * @param[in] b ��ʌ`�̌W����b
 * @param[in] c ��ʌ`�̌W����c
 * @param[in] d ��ʌ`�̌W����d
 * @param[in] point �_
 * @return �_�ƕ��ʂƂ̋���
 */
double DistancePlaneToDepth(double a, double b, double c, double d, double point[3]) {
	double x = point[0];
	double y = point[1];
	double z = point[2];

	return fabs(a * x + b * y + c * z + d) / (sqrt(a * a + b * b + c * c) + DBL_MIN);
}

/**
 * 2�̃x�N�g������Ȃ镽�ʂ��悢���ʂ��𔻒肷��
 * @param[in] vec1 �x�N�g��1
 * @param[in] vec2 �x�N�g��2
 * @param[in] n �x�N�g���̎�����
 * @retval 1 �悢���� (�@���x�N�g�����o�\)
 * @retval 0 �������� (�@���x�N�g�����o�s�\�ȋ��ꂪ����)
 */
int IsCorrectPlane(const double *vec1, const double *vec2, size_t n) {
	double angle;
	double t = pi() / 6.0;

	/* �Ȃ��p�����߂� */
	angle = IncludedAngle(vec1, vec2, n);

	/* �Ȃ��p���傫�����Ă����������Ă��悢���ʂƂ͂����Ȃ� */
	if ( angle > t && angle < pi() - t ) {
		return 1;
	} else {
		return 0;
	}
}

/**
 * ��1�ی������4�ی��܂ł̕Ίp�����߂�
 * @param[in] y y�������̍��W
 * @param[in] x x�������̍��W
 * @return xy���ʂ̕Ίp(0 <= Arg <= 2pi)
 */
double Arg(double y, double x) {
	double pi = atan(1.0) * 4.0;

	if (x >= 0 && y >= 0) { // ���ی�
		return atan2(y, x);

	} else if (x <= 0 && y >= 0) { // ���ی�
		return atan2(y, x);

	} else if (x <= 0 && y <= 0 ) { // ��O�ی�
		return 2.0 * pi + atan2(y, x);

	} else { // ��l�ی�
		return 2.0 * pi + atan2(y, x);

	}
}

/**
 * �V���p�����߂�
 * @param[in] x x�������̍��W
 * @param[in] y y�������̍��W
 * @param[in] z z�������̍��W
 * @return xyz��Ԃ̓V���p(0 <= Arg <= pi)
 */
double ZenithAngle(double x, double y, double z) {
	double pi = atan(1.0) * 4.0;
	double angle;

	angle = atan2(z, sqrt(x * x + y * y));
	return pi / 2.0 - angle;
}

/**
 * �V�O���C�h�֐����v�Z����
 * @param[in] x ����
 * @param[in] gain �Q�C��
 * @return �V�O���C�h�֐�
 */
double Sigmoid(double x, double gain) {
  return 1.0 / (1.0 + exp(-gain * x));
}

/**
 * �~���������߂�
 * @return �~����
 */
double pi(void) {
    return atan(1.0) * 4.0;
}

/**
 * �ʓx�@�\�L��x���@�\�L�ɕϊ�����
 * @param[in] r �p�x[rad]
 * @return �p�x[deg]
 */
double to_deg(double r) {
    return r * 180.0 / (atan(1.0) * 4.0);
}

/**
 * �x���@�\�L���ʓx�@�\�L�ɕϊ�����
 * @param[in] r �p�x[deg]
 * @return �p�x[rad]
 */
double to_rad(double a) {
    return a * atan(1.0) * 4.0 / 180.0;
}

/**
 * �C�ӂ̐������Ƃ���ΐ����v�Z����
 * @param[in] base ��
 * @param[in] antilog �^��
 */
double logn(double base, double antilog) {
    return log(antilog) / log(base);
}

/**
 * �j���[�g���@�ŗ��������ߎ�����
 * @param[in] a ����
 * @param[in] x ���ɋ߂����Ȓl
 * @return a �̗�����
 */
double cbrt_newton(double a, double x) {
    double e;

    do {
        e = (x * x * x - a) / (3.0 * x * x);
        x = x - e;
    } while ( fabs(e) > 1.0e-16 );

    return x;
}

double frac(double x, double y){
	return x / (y + DBL_MIN);
}