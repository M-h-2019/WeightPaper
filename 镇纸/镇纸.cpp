#include <iostream>
#include <cstdio>
#include <cstring>
#include <vector>
#include <math.h>
using namespace std;
#define PB push_back
double fabs(double num){		//计算|num| 
	return num < 0 ? -num : num; 
} 

double fmin(double num1, double num2 ){		//计算num1和num2的较小者
	return num1 <= num2 ? num1 : num2 ; 
}

double fmax(double num1, double num2){		//计算num1和num2的较大者
	return num1 >= num2 ? num1 : num2 ; 
}

double Sqr(double num){						//计算num^2 
	return num * num;
}

struct Point{								//空间向量点坐标的结构类型 
	double X;
	double Y;
	double Z;
};
bool operator == (Point A, Point B){ 		//返回空间向量A和B重合的标志
	return (A.X == B.X && A.Y == B.Y && A.Z == B.Z);		
} 
Point MakePoint(double X0, double Y0, double Z0){
											//构建坐标为（X0,Y0,Z0)的空间向量点
	Point Ret;
	Ret.X = X0;
	Ret.Y = Y0;
	Ret.Z = Z0; 
	return Ret; 
}
Point operator + (Point A, Point B){		//返回向量A+B
	return MakePoint(A.X + B.X, A.Y + B.Y, A.Z + B.Z ); 
}
Point operator - (Point A, Point B){		//返回向量A-B
	return MakePoint(A.X - B.X, A.Y - B.Y, A.Z - B.Z ); 
}
double operator * (Point A, Point B){		//返回数量积A·B 
	return A.X * B.X + A.Y * B.Y + A.Z * B.Z ; 
}
Point operator ^ (Point A,Point B){			//返回向量积A^B
	return MakePoint(A.Y * B.Z - A.Z * B.Y , A.Z * B.X - A.X * B.Z , A.X * B.Y - A.Y * B.X ); 
}
Point operator * (Point A, double x) {		//向量的数乘运算
	return MakePoint(A.X * x, A.Y * x, A.Z * x); 
}
Point operator / (Point A, double x){		//向量的数除运算
	return MakePoint(A.X / x, A.Y / x, A.Z / x);
}
double Volume(Point A, Point B, Point C, Point D){	//返回四面体体积VABCD
	return fabs(((B - A) ^ (C - A)) * (D - A) ) / 6.000 ; 
}
double Distance(Point A, Point B){			//返回空间向量点A和B间的欧式距离 
	return sqrt(Sqr(A.X - B.X) + Sqr(A.Y - B.Y ) + Sqr(A.Z - B.Z)); 
}
double Abs(Point p){						//返回向量p的模
	return sqrt(Sqr(p.X) + Sqr(p.Y) + Sqr(p.Z)); 
}
vector <Point> Flat;						//平面的点集
Point A , B , C , D , E , F , G , G_ABCD , G_ABCE , P[10];
double V_ABCD , V_ABCE ; 
void FindCentre(){							//计算镇纸的重心
	V_ABCD = Volume(A , B , C , D);			
	V_ABCE = Volume(A , B , C , E);
	
	G_ABCD = (A + B + C + D)/ 4.000;
	G_ABCE = (A + B + C + E)/ 4.000;
	
	G=G_ABCD + ((G_ABCE - G_ABCD) * V_ABCE) / (V_ABCD + V_ABCE) ; 
} 
bool SameSide(Point p1, Point p2 , Point p3){	//返回所有点位于平面（由p1,p2,p3构成）
												//同侧的标志 
	Point Nor = (p2 - p1) ^ (p3 - p1);			//计算平面的法向量
	bool plus = false , minus = false;
	for(int i = 1; i <= 5; i ++)				//判断平面两侧是否存在点
	{
		if(P[i] == p1 || P[i] == p2 || P[i] == p3)
		continue;								//排除重合的点
		double temp = (P[i] - p1) * Nor;
		if(temp < 0)
		minus = true;
		else if(temp > 0)
		plus = true;
	}
	if(plus && minus)
	return false;								//若平面两侧存在点，则失败退出
	
	Flat.clear() ;								//平面上的点集Flat初始化
	Flat.PB(p1) ; Flat.PB(p2) ; Flat.PB(p3);	//p1,p2,p3压入Flat
	for(int i = 1; i <= 5; i++)					//将位于平面上的点压入Flat
	{
		if(P[i] == p1 || P[i] == p2 || P[i] == p3)
		continue;
		double temp = (P[i] - p1) * Nor;
		if(temp == 0)
		Flat.PB(P[i]);
	}
	return true;
} 
Point Shadow(Point p , Point p1, Point p2, Point p3)	//返回点p在平面（由p1,p2,p3构成）上的射影点
{
	Point Nor = (p2 - p1) ^ (p3 - p1);			//计算平面度的法向量
	Point Nor0 = (Nor * (Nor * (p - p1) )/ Abs(Nor) / Abs(Nor));
	return (p - Nor0);
}
double ShadowHeight(Point p, Point p1, Point p2, Point p3)
												//返回p与p到平面（由p1,p2,p3构成）射影的距离
{
	Point Nor = (p2 - p1) ^ (p3 - p1);			//计算平面的法向量
	return fabs(Nor * (p - p1)) / Abs(Nor) ; 
}
bool InTheMiddle(Point A, Point B, Point C, Point P)	//返回空间向量AP夹在AB和AC间的标志
{
	Point AB = B - A;
	Point AC = C - A;
	Point AP = P - A;
	
	Point v1 = AB ^ AC;
	Point v2 = AB ^ AP;
	
	return (v1 * v2) >= 0 ;
}
bool PointinTriangle(Point A, Point B, Point C, Point P)
												//返回点p在三角形ABCn内的标志
{
	return InTheMiddle(A, B, C, P) && InTheMiddle(B, C, A, P) && InTheMiddle(C, A, B, P);
}
bool PointInArea(Point p)					//返回点p在平面区域Flat内的标志
{
	for (int i = 0 ; i < Flat.size() ; i ++)	//枚举平面区域Flat内3个点的所有可能组合
												//一旦发现点p在某三点构成的三角形内部，
												//则说明p在平面区域Flat内 
		for(int j = i + 1 ; j < Flat.size() ; j ++)
			for (int k = j + 1 ; k < Flat.size() ; k ++)
				if(PointinTriangle(Flat[i], Flat[j], Flat[k], p))
				return true;
			return false;					//点p在平面区域Flat外 
} 
double DistanceToLine(Point P, Point A, Point B)	//计算点P到直线AB的距离
{
	if((B - A) * (P - A) > 0.000 && (A - B) * (P - B) > 0.000)
		return sqrt(Sqr(Distance(A, P)) - Sqr((B - A) * (P - A) / Abs((B - A))));
	return fmin(Distance(A, P), Distance(B, P));
}
bool Balance(Point p)						//返回稳定标志（点p向任意方向移动0.2个单位后仍然在平面区域Flat内）
{
	if( !PointInArea(p))
	   return false;
	for(int i = 0; i < Flat.size() ; i ++)	//枚举Flat上的所有线段
		for(int j = i + 1; j < Flat.size() ; j ++){
			Point P_i = Flat[i], P_j = Flat[j] ;
			int k;
			for( k = 0 ; k < Flat.size() ; k ++ )	//寻找平面上线段P_i，P_j外的一个点P_k
				if(k != i && k != j)
				break;
			Point P_k = Flat[k];
			
			bool plus = false, minus = false;
			Point Nor = (P_j - P_i) ^ (P_k - P_i) ^ (P_j - P_i);	//计算叉积
			for (int l = 0 ; l < Flat.size(); l ++ ) if ( l != i && l != j)
			{
				if(Nor * (Flat[l] - P_i ) > 0 )
				plus = true;
				else if(Nor * (Flat[l] - P_i) < 0 )
				minus = true ;
			 } 
			 if (plus && minus)
			 continue;				//P_i,P_j不是凸包上的边
			 if (DistanceToLine(p , P_i, P_j) < 0.2)
			 						//若p到凸包边缘距离小于0.2，说明以p为圆心、0.2为半径的圆不在凸壳面内，返回不稳定标志
				return false;
		} 
		return true;
}

void FindHull()						//寻找三维凸包的所有面，计算和输出芯片至底面的最大距离和最小距离
{
	double Min_Answer = 2000.000, Max_Answer = 0.000;
	P[1] = A ; P[2] = B ; P [3] = C ; P [4] = D ; P[5] = E ; P [6] = F;
	
	for(int i = 1 ; i <= 5; i ++)	//枚举5个点中选3个点的所有组合
		for(int j = i + 1 ; j <= 5 ; j ++)
			for(int k = j + 1 ; k <= 5 ; k ++){
				if(SameSide(P[i], P[j], P[k]) )	//若所有点在平面（由点P[i]P[j]P[k]构成）的
												//同侧，则说明该平面可作为镇纸的底面
				{
					Point G0 = Shadow(G, P[i], P[j], P [k]);	//计算重心G在底面的射影点G0
					if(Balance(G0))						//若处于稳定状态（G0向任意方向移动0.2个单位
														//后仍然在底面范围内），则计算芯片点F与F底面
														//射影的距离，调整最大距离和最小距离
					{
						Min_Answer = fmin(Min_Answer, ShadowHeight(F, P[i], P[j], P[k]));
						Max_Answer = fmax(Max_Answer, ShadowHeight(F, P[i], P[j], P[k])); 
			        }
			    }
			}
		cout << fixed << Min_Answer << " " << Max_Answer << "\n" ;
														//输出最大距离和最小距离
}
int main(){
	cout.precision(5);			//定义精度为小数点后五位
	int Case = 0, deg ;		//测试用例编号初始化
	char ss[200];
	while( gets(ss))				//反复输入镇纸的5个坐标和芯片坐标，直至输入0为止
	{
		if(strlen(ss) == 1 && ss[0] == '0' )
		break;
		sscanf(ss,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
			&A.X , &A.Y , &A.Z , &B.X , &B.Y , &B.Z ,
			&C.X , &C.Y , &C.Z , &D.X , &D.Y , &D.Z ,
			&E.X , &E.Y , &E.Z , &F.X , &F.Y , &F.Z );
		cout << "Case " << ++ Case << ": " ;	//输出测试用例
		FindCentre();							//计算镇纸的重心
		FindHull();								//计算和输出芯片至底面的最大距离和最小距离
	}	 
}

