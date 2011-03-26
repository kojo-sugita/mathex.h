
/** 
 * @file    mathex.c
 * @brief   math.hの拡張
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
 * 総和を計算する
 * @param[in] a 総和を計算する配列
 * @param[in] n 配列の要素数
 * @return 総和
 */
double Summation(const double *a, size_t n) {
	double sum = 0.0;

	while (n--) 
		sum += a[n];

	return sum;
}

/**
 * 平方和を計算する
 * @param[in] a 平方和を計算する配列
 * @param[in] n 配列の要素数
 * @return 平方和
 */
double SumOfSquares(const double *a, size_t n) {
	double sum = 0.0;

	while (n--)
		sum += (a[n] * a[n]);

	return sum;
}

/**
 * 積和を計算する
 * @param[in] a1 積和を計算する配列1
 * @param[in] a2 積和を計算する配列2
 * @param[in] n 配列の要素数
 * @return 積和
 */
double SumOfProduct(const double *a1, const double *a2, size_t n) {
	double sum = 0.0;

	while (n--)
		sum += (a1[n] * a2[n]);

	return sum;
}

/**
 * 平均値を計算する
 * @param[in] a 入力配列
 * @param[in] n 配列の要素数
 * @return 平均値
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
 * 分散を計算する
 * @param[in] a 入力配列
 * @param[in] n 配列の要素数
 * @return 分散
 */
double Variance(const double *a, size_t n) {
	size_t i;
	double sum = 0.0;
	double sqsum = 0.0;
	double avg;

	if ( n < 1 ) return 0.0;

	for (i = 0; i < n; i++) {
		sum += a[i]; /* 和 */
		sqsum += a[i] * a[i]; /* 平方和 */
	}

	avg = sum / (double)n; /* 平均 */
	return sqsum / (double)n - avg * avg; /* 分散 */
}

/**
 * 不偏分散を計算する
 * @param[in] a 入力配列
 * @param[in] n 配列の要素数
 * @return 不偏分散
 */
double UnbiasedVariance(const double *a, size_t n) {
    size_t i;
    double sum = 0.0;
    double sqsum = 0.0;
    double avg;

    if ( n < 2 ) return 0.0;

    for (i = 0; i < n; i++) {
        sum += a[i]; /* 和 */
        sqsum += a[i] * a[i]; /* 平方和 */
    }

    avg = sum / (double)n; /* 平均 */
    return sqsum / (double)(n - 1) - avg * avg; /* 不偏分散 */
}

/**
 * 共分散を計算する
 * @param[in] x 入力配列1
 * @param[in] y 入力配列2
 * @param[in] n xとyの要素数
 * @return 共分散
 */
double Covariance(const double *x, const double *y, int n) {
    int i;
    double sum_x = 0.0, sum_y = 0.0;
    double avg_x, avg_y;
    double sumprod = 0.0;

    if ( n < 1 ) return 0.0;

    for (i = 0; i < n; i++) {
        sum_x += x[i]; /* xの和 */
        sum_y += y[i]; /* yの和 */
        sumprod += x[i] * y[i];
    }

    avg_x = sum_x / (double)n; /* xの平均 */
    avg_y = sum_y / (double)n; /* yの平均 */
    return sumprod / (double)n - avg_x * avg_y; /* 共分散 */
}

/**
 * 標準偏差を計算する
 * @param[in] a 入力配列
 * @param[in] n 配列の要素数
 * @return 標準偏差
 */
double StdDev(const double *a, size_t n) {
	size_t i;
	double sum = 0.0;
	double sqsum = 0.0;
	double avg, var;

	if ( n < 1 ) return 0.0;

	for (i = 0; i < n; i++) {
		sum += a[i]; /* 和 */
		sqsum += a[i] * a[i]; /* 平方和 */
	}

	avg = sum / (double)n; /* 平均 */
	var = sqsum / (double)n - avg * avg; /* 分散 */

	if ( var < 0.0 ) var = 0.0;
	return sqrt(var); /* 標準偏差 */

}

/**
 * 最大値を求める
 * @param[in] a 入力配列
 * @param[in] n 配列の要素数
 * @return 最大値
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
 * 最大値を求める(int型配列)
 * @param[in] a 入力配列
 * @param[in] n 配列の要素数
 * @return 最大値
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
 * 最小値を求める
 * @param[in] a 入力配列
 * @param[in] n 配列の要素数
 * @return 最小値
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
 * 最小値を求める(int型配列)
 * @param[in] a 入力配列
 * @param[in] n 配列の要素数
 * @return 最小値
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
 * 最大値と最小値を求める
 * @param[in] a 入力配列
 * @param[in] n 配列の要素数
 * @param[out] max 最大値
 * @param[out] min 最小値
 * @retval 1 最大値と最小値を取得
 * @retval 0 配列の要素が0
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
 * 階乗を計算する
 * @param[in] n 整数 
 * @return 階乗
 */
int Factorial(int n) {
	int fact = n;

	if (n <= 0) return 1;

	while (--n) 
		fact *= n;

	return fact;
}

/**
 * 順列組み合わせの総数を計算する
 * @param[in] n nPrのn 
 * @param[in] r nPrのr 
 * @return 順列組み合わせの総数
 */
int Permutation(int n, int r) {
	if (r <= 0) return 1;
	/* 0除算はありえない */
	return Factorial(n) / Factorial(n - r);
}

/**
 * 組み合わせの総数を計算する
 * @param[in] n nCrのn 
 * @param[in] r nCrのr 
 * @return 組み合わせの総数
 */
int Combination(int n, int r) {
	if (r <= 0) return 1;
	/* 0除算はありえない */
	return Factorial(n) / (Factorial(r) * Factorial(n - r));
}


/**
 * 2点間の距離を求める
 * @param[in] p1 点の座標1
 * @param[in] p2 点の座標2
 * @param[in] n 次元数
 * @return p1とp2間の距離
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
 * 2点p1,p2からなる線分をm:nに内分する点を求める (2次元)
 * @param[in] p1 点の座標1
 * @param[in] p2 点の座標2
 * @param[in] m m:nのm
 * @param[in] n m:nのn
 * @param[out] result 内分点
 * @retval 1 正常終了
 * @retval 0 異常終了
 */
int DivideInternally2D(const double p1[2], const double p2[2], int m, int n, double result[2]) {
	if (m < 1 && n < 1) return 0;
	
	result[0] = (n * p1[0] + m * p2[0]) / (m + n);
	result[1] = (n * p1[1] + m * p2[1]) / (m + n);

	return 1;
}

/**
 * 2点p1,p2からなる線分をm:nに外分する点を求める (2次元)
 * @param[in] p1 点の座標1
 * @param[in] p2 点の座標2
 * @param[in] m m:nのm
 * @param[in] n m:nのn
 * @param[out] result 外分点
 * @retval 1 正常終了
 * @retval 0 異常終了
 */
int DivideExternally2D(const double p1[2], const double p2[2], int m, int n, double result[2]) {
	if (m - n == 0) return 0;

	result[0] = (-n * p1[0] + m * p2[0]) / (m - n);
	result[1] = (-n * p1[1] + m * p2[1]) / (m - n);

	return 1;
}

/**
 * 3点p1,p2,p3からなる三角形の重心を求める (2次元)
 * @param[in] p1 点の座標1
 * @param[in] p2 点の座標2
 * @param[in] p3 点の座標3
 * @param[out] result 三角形の重心
 * @return resultの先頭アドレス
 */
double *Triangle_CenterOfGravity(const double p1[2], const double p2[2], const double p3[2], double result[2]) {
	result[0] = (p1[0] + p2[0] + p3[0]) / 3.0;
	result[1] = (p1[1] + p2[1] + p3[1]) / 3.0;
	return result;
}

/**
 * 2点を結ぶベクトルを求める
 * @param[in] p1 点の座標1
 * @param[in] p2 点の座標2
 * @param[in] n 次元数
 * @param[out] vec ベクトル
 * @return vec の先頭アドレス
 */
double *ToVector(const double *p1, const double *p2, size_t n, double *vec) {
	size_t i;
	
	for ( i = 0; i < n; i++ ) {
		vec[i] = p2[i] - p1[i];
	}
	return vec;
}

/**
 * ベクトルの和を求める
 * @param[in] vec1 ベクトル1
 * @param[in] vec2 ベクトル2
 * @param[out] sumVec ベクトルの和
 * @param[in] n 次元数
 * @return sumVecの先頭アドレス
 */
double *VectorSum(const double *vec1, const double *vec2, double *sumVec, size_t n) {
	while (n--) 
		sumVec[n] = vec1[n] + vec2[n];
	return sumVec;
}

/**
 * ベクトルの差を求める
 * @param[in] vec1 ベクトル1
 * @param[in] vec2 ベクトル2
 * @param[out] sumVec ベクトルの和
 * @param[in] n 次元数
 * @return diffVecの先頭アドレス
 */
double *VectorDifference(const double *vec1, const double *vec2, double *diffVec, size_t n) {
	while (n--) 
		diffVec[n] = vec1[n] - vec2[n];
	return diffVec;
}

/**
 * ベクトルの長さ(ノルム)を計算する
 * @param[in] vec ベクトル
 * @param[in] n ベクトルの次元数
 * @return vec の長さ
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
 * ベクトルの内積を計算する
 * @param[in] vec1 ベクトル1
 * @param[in] vec2 ベクトル2
 * @param[in] n ベクトルの次元数
 * @return vec1 と vec2 の内積
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
 * 2 つのベクトルのなす角を計算する
 * @param[in] vec1 ベクトル1
 * @param[in] vec2 ベクトル2
 * @param[in] n ベクトルの次元数
 * @return vec1 と vec2 のなす角
 */
double IncludedAngle(const double *vec1, const double *vec2, size_t n) {
    return acos(InnerProduct(vec1, vec2, n) / 
        (Norm(vec1, n) * Norm(vec2, n) + DBL_MIN));
}

/**
 * 単位ベクトルを求める
 * @param[in] vec ベクトル
 * @param[out] unitVec 単位ベクトル
 * @param[in] n ベクトルの次元数
 * @retval 1 単位ベクトルを正常に取得
 * @retval 0 ベクトルがゼロベクトル
 */
int UnitVector(const double *vec, double *unitVec, size_t n) {
	size_t i;
	double norm = Norm(vec, n);

	/* ベクトルのノルムが0でないか? */
	if (norm > 0.0) { 
		/* 単位ベクトルを取得 */
		for (i = 0; i < n; i++) {
			unitVec[i] = vec[i] / Norm(vec, n);
		}
		return 1;

	} else {
		/* 結果にすべて0を格納する */
		for (i = 0; i < n; i++) {
			unitVec[i] = 0;
		}
		return 0;
	}
}

/**
 * 3次元の法線ベクトルを求める
 * @param[in] p1 座標1
 * @param[in] p2 座標2
 * @param[in] p3 座標3
 * @param[in] result 法線ベクトル
 * @return resultの先頭アドレス
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
 * 3次元ベクトルの外積を計算する
 * @param[in] vec1 ベクトル1
 * @param[in] vec2 ベクトル2
 * @param[out] result ベクトルの外積
 * @return resultの先頭アドレス
 */
double *CrossProduct3D(const double *vec1, const double *vec2, double *result) {
	result[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1];
	result[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2];
	result[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];

	return result;
}

/**
 * 2 つのベクトルの類似度を計算する
 * @param[in] vec1 ベクトル1
 * @param[in] vec2 ベクトル2
 * @param[in] n ベクトルの次元数
 * @return vec1 と vec2 の類似度 (-1.0 <= return <= 1.0) 1.0に近いほど類似度大
 * @attention ゼロベクトルが引数に渡された場合0を返す
 */
double VectorSimilarity(const double *vec1, const double *vec2, size_t n) {
    return InnerProduct(vec1, vec2, n) / 
		((Norm(vec1, n) * Norm(vec2, n)) + DBL_MIN);
}

/**
 * 正規化距離を計算する
 * @param[in] vec1 ベクトル1
 * @param[in] vec2 ベクトル2
 * @param[in] n ベクトルの次元数
 * @return vec1 と vec2 の正規化距離
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
 * 3次元空間上の点p1と点p2を結ぶ直線上に点pが存在するかを判定する
 * @param[in] p1 点1
 * @param[in] p2 点2
 * @param[in] p 3次元空間上のある点
 * @retval 1 p1とp2を結ぶ直線上にpが存在する
 * @rarval 0 p1とp2を結ぶ直線上にpが存在しない
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
 * 3次元平面を作成する
 * @param[in] AB 3次元ベクトル1
 * @param[in] AC 3次元ベクトル2
 * @param[in] A 座標1
 * @param[out] a 平面の方程式 ax + by + cz + d = 0の係数の a
 * @param[out] b 平面の方程式 ax + by + cz + d = 0 の係数の b 
 * @param[out] c 平面の方程式 ax + by + cz + d = 0 の係数の c
 * @param[out] d 平面の方程式 ax + by + cz + d = 0 の係数の d
 */
void CalculatePlane3D(const double *AB, const double *AC, const double *A, double *a, double *b, double *c, double *d) {
	double norm;

	/* 法線ベクトル */
	*a = AB[1] * AC[2] - AB[2] * AC[1];
	*b = AB[2] * AC[0] - AB[0] * AC[2];
	*c = AB[0] * AC[1] - AB[1] * AC[0];

	norm = sqrt(*a * *a + *b * *b + *c * *c);
	*a = *a / (norm + DBL_MIN);
	*b = *b / (norm + DBL_MIN);
	*c = *c / (norm + DBL_MIN);

	/* 原点から平面に伸ばした垂線の長さ */
	*d = -1.0 * *a * A[0] + -1.0 * *b * A[1] + -1.0 * *c * A[2];
}

/**
 * 3点から3次元平面を作成する
 * @param[in] P1 点1
 * @param[in] P2 点2
 * @param[in] P3 点3
 * @param[out] a 平面の方程式 ax + by + cz + d = 0の係数の a
 * @param[out] b 平面の方程式 ax + by + cz + d = 0 の係数の b 
 * @param[out] c 平面の方程式 ax + by + cz + d = 0 の係数の c
 * @param[out] d 平面の方程式 ax + by + cz + d = 0 の係数の d
 */
void CalculatePlane3D_Using3Depths(double P1[3], double P2[3], double P3[3], double *a, double *b, double *c, double *d) {
	double AB[3], AC[3];
	double norm;

	/* ベクトルに変換 */
	ToVector(P1, P2, 3, AB);
	ToVector(P1, P3, 3, AC);

	/* 法線ベクトル */
	*a = AB[1] * AC[2] - AB[2] * AC[1];
	*b = AB[2] * AC[0] - AB[0] * AC[2];
	*c = AB[0] * AC[1] - AB[1] * AC[0];

	norm = sqrt(*a * *a + *b * *b + *c * *c);
	*a = *a / (norm + DBL_MIN);
	*b = *b / (norm + DBL_MIN);
	*c = *c / (norm + DBL_MIN);

	/* 原点から平面に伸ばした垂線の長さ */
	*d = -1.0 * *a * P1[0] + -1.0 * *b * P1[1] + -1.0 * *c * P1[2];
}

/**
 * 直線の方程式を標準形から一般形に変形する (y = ax + b --> ax + by + c = 0)
 * @param[in] a 標準形の係数 (傾き)
 * @param[in] b 標準形の係数 (切片)
 * @param[out] ga 一般形の係数のa
 * @param[out] gb 一般形の係数のb
 * @param[out] gc 一般形の係数のc
 * @return EXIT_SUCCESS
 */
int ToLineCanonicalForm(double a, double b, double *ga, double *gb, double *gc)
{
	double norm;

	/* 法線ベクトル */
	norm = sqrt(a * a + 1.0);
	*ga = a / (norm + DBL_MIN);
	*gb = -1.0 / (norm + DBL_MIN);

	/* 原点から平面に伸ばした垂線の長さ */
	*gc = b;
	return EXIT_SUCCESS;
}


/**
 * 平面の方程式を標準形から一般形に変形する (z = ax + by + c --> ax + by + cz + d = 0)
 * @param[in] a 標準形の係数a
 * @param[in] b 標準形の係数b
 * @param[in] c 標準形の係数c
 * @param[out] ga 一般形の係数のa
 * @param[out] gb 一般形の係数のb
 * @param[out] gc 一般形の係数のc
 * @param[out] gd 一般形の係数のd
 * @return EXIT_SUCCESS
 */
int ToPlaneCanonicalForm(double a, double b, double c, double *ga, double *gb, double *gc, double *gd) {

	/* 法線ベクトル */
	//*ga = a;  /* 係数a */
	//*gb = b;  /* 係数b */
	//*gc = -1; /* 係数c */

	double norm;

	norm = sqrt(a * a + b * b + 1.0);
	*ga = a / (norm + DBL_MIN);
	*gb = b / (norm + DBL_MIN);
	*gc = -1.0 / (norm + DBL_MIN);

	/* 原点から平面に伸ばした垂線の長さ */
	*gd = c;  /* 係数d */

	return EXIT_SUCCESS;
}

/**
 * 平面の方程式を一般形から標準形に変形する (z = ax + by + c --> ax + by + cz + d = 0)
 * @param[in] a 一般形の係数のa
 * @param[in] b 一般形の係数のb
 * @param[in] c 一般形の係数のc
 * @param[in] d 一般形の係数のd
 * @param[out] ca 標準形の係数a
 * @param[out] cb 標準形の係数b
 * @param[out] cc 標準形の係数c
 * @retval EXIT_SUCCESS 変換成功
 * @retval EXIT_FAILURE 変換失敗(0除算)
 */
int ToPlaneGenericForm(double a, double b, double c, double d, double *ca, double *cb, double *cc) {
	if (c != 0) {
		/* 傾き */
		*ca = - a / c; /* 係数a */
		*cb = - b / c; /* 係数b */

		/* 切片 */
		*cc = - d / c; /* 係数c */
		return EXIT_SUCCESS;

	} else {
		*ca = 0;
		*cb = 0;
		*cc = 0;
		return EXIT_FAILURE;
	}
}

/**
 * 点と直線との距離を返す
 * @param[in] a 一般形の係数のa
 * @param[in] b 一般形の係数のb
 * @param[in] c 一般形の係数のc
 * @param[in] point 点
 * @return 点と直線との距離
 */
double DistanceLineToDepth(double a, double b, double c, double point[2]) {
	double x = point[0];
	double y = point[1];

	return fabs(a * x + b * y + c) / (sqrt(a * a + b * b) + DBL_MIN);
}

/**
 * 点と平面との距離を返す
 * @param[in] a 一般形の係数のa
 * @param[in] b 一般形の係数のb
 * @param[in] c 一般形の係数のc
 * @param[in] d 一般形の係数のd
 * @param[in] point 点
 * @return 点と平面との距離
 */
double DistancePlaneToDepth(double a, double b, double c, double d, double point[3]) {
	double x = point[0];
	double y = point[1];
	double z = point[2];

	return fabs(a * x + b * y + c * z + d) / (sqrt(a * a + b * b + c * c) + DBL_MIN);
}

/**
 * 2つのベクトルからなる平面がよい平面かを判定する
 * @param[in] vec1 ベクトル1
 * @param[in] vec2 ベクトル2
 * @param[in] n ベクトルの次元数
 * @retval 1 よい平面 (法線ベクトル検出可能)
 * @retval 0 悪い平面 (法線ベクトル検出不可能な恐れがある)
 */
int IsCorrectPlane(const double *vec1, const double *vec2, size_t n) {
	double angle;
	double t = pi() / 6.0;

	/* なす角を求める */
	angle = IncludedAngle(vec1, vec2, n);

	/* なす角が大きすぎても小さすぎてもよい平面とはいえない */
	if ( angle > t && angle < pi() - t ) {
		return 1;
	} else {
		return 0;
	}
}

/**
 * 第1象限から第4象限までの偏角を求める
 * @param[in] y y軸方向の座標
 * @param[in] x x軸方向の座標
 * @return xy平面の偏角(0 <= Arg <= 2pi)
 */
double Arg(double y, double x) {
	double pi = atan(1.0) * 4.0;

	if (x >= 0 && y >= 0) { // 第一象限
		return atan2(y, x);

	} else if (x <= 0 && y >= 0) { // 第二象限
		return atan2(y, x);

	} else if (x <= 0 && y <= 0 ) { // 第三象限
		return 2.0 * pi + atan2(y, x);

	} else { // 第四象限
		return 2.0 * pi + atan2(y, x);

	}
}

/**
 * 天頂角を求める
 * @param[in] x x軸方向の座標
 * @param[in] y y軸方向の座標
 * @param[in] z z軸方向の座標
 * @return xyz空間の天頂角(0 <= Arg <= pi)
 */
double ZenithAngle(double x, double y, double z) {
	double pi = atan(1.0) * 4.0;
	double angle;

	angle = atan2(z, sqrt(x * x + y * y));
	return pi / 2.0 - angle;
}

/**
 * シグモイド関数を計算する
 * @param[in] x 引数
 * @param[in] gain ゲイン
 * @return シグモイド関数
 */
double Sigmoid(double x, double gain) {
  return 1.0 / (1.0 + exp(-gain * x));
}

/**
 * 円周率を求める
 * @return 円周率
 */
double pi(void) {
    return atan(1.0) * 4.0;
}

/**
 * 弧度法表記を度数法表記に変換する
 * @param[in] r 角度[rad]
 * @return 角度[deg]
 */
double to_deg(double r) {
    return r * 180.0 / (atan(1.0) * 4.0);
}

/**
 * 度数法表記を弧度法表記に変換する
 * @param[in] r 角度[deg]
 * @return 角度[rad]
 */
double to_rad(double a) {
    return a * atan(1.0) * 4.0 / 180.0;
}

/**
 * 任意の整数を底とする対数を計算する
 * @param[in] base 底
 * @param[in] antilog 真数
 */
double logn(double base, double antilog) {
    return log(antilog) / log(base);
}

/**
 * ニュートン法で立方根を近似する
 * @param[in] a 実数
 * @param[in] x 解に近そうな値
 * @return a の立方根
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