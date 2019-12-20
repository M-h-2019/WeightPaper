#include <iostream>
#include <cstdio>
#include <cstring>
#include <vector>
#include <math.h>
using namespace std;
#define PB push_back
double fabs(double num){		//����|num| 
	return num < 0 ? -num : num; 
} 

double fmin(double num1, double num2 ){		//����num1��num2�Ľ�С��
	return num1 <= num2 ? num1 : num2 ; 
}

double fmax(double num1, double num2){		//����num1��num2�Ľϴ���
	return num1 >= num2 ? num1 : num2 ; 
}

double Sqr(double num){						//����num^2 
	return num * num;
}

struct Point{								//�ռ�����������Ľṹ���� 
	double X;
	double Y;
	double Z;
};
bool operator == (Point A, Point B){ 		//���ؿռ�����A��B�غϵı�־
	return (A.X == B.X && A.Y == B.Y && A.Z == B.Z);		
} 
Point MakePoint(double X0, double Y0, double Z0){
											//��������Ϊ��X0,Y0,Z0)�Ŀռ�������
	Point Ret;
	Ret.X = X0;
	Ret.Y = Y0;
	Ret.Z = Z0; 
	return Ret; 
}
Point operator + (Point A, Point B){		//��������A+B
	return MakePoint(A.X + B.X, A.Y + B.Y, A.Z + B.Z ); 
}
Point operator - (Point A, Point B){		//��������A-B
	return MakePoint(A.X - B.X, A.Y - B.Y, A.Z - B.Z ); 
}
double operator * (Point A, Point B){		//����������A��B 
	return A.X * B.X + A.Y * B.Y + A.Z * B.Z ; 
}
Point operator ^ (Point A,Point B){			//����������A^B
	return MakePoint(A.Y * B.Z - A.Z * B.Y , A.Z * B.X - A.X * B.Z , A.X * B.Y - A.Y * B.X ); 
}
Point operator * (Point A, double x) {		//��������������
	return MakePoint(A.X * x, A.Y * x, A.Z * x); 
}
Point operator / (Point A, double x){		//��������������
	return MakePoint(A.X / x, A.Y / x, A.Z / x);
}
double Volume(Point A, Point B, Point C, Point D){	//�������������VABCD
	return fabs(((B - A) ^ (C - A)) * (D - A) ) / 6.000 ; 
}
double Distance(Point A, Point B){			//���ؿռ�������A��B���ŷʽ���� 
	return sqrt(Sqr(A.X - B.X) + Sqr(A.Y - B.Y ) + Sqr(A.Z - B.Z)); 
}
double Abs(Point p){						//��������p��ģ
	return sqrt(Sqr(p.X) + Sqr(p.Y) + Sqr(p.Z)); 
}
vector <Point> Flat;						//ƽ��ĵ㼯
Point A , B , C , D , E , F , G , G_ABCD , G_ABCE , P[10];
double V_ABCD , V_ABCE ; 
void FindCentre(){							//������ֽ������
	V_ABCD = Volume(A , B , C , D);			
	V_ABCE = Volume(A , B , C , E);
	
	G_ABCD = (A + B + C + D)/ 4.000;
	G_ABCE = (A + B + C + E)/ 4.000;
	
	G=G_ABCD + ((G_ABCE - G_ABCD) * V_ABCE) / (V_ABCD + V_ABCE) ; 
} 
bool SameSide(Point p1, Point p2 , Point p3){	//�������е�λ��ƽ�棨��p1,p2,p3���ɣ�
												//ͬ��ı�־ 
	Point Nor = (p2 - p1) ^ (p3 - p1);			//����ƽ��ķ�����
	bool plus = false , minus = false;
	for(int i = 1; i <= 5; i ++)				//�ж�ƽ�������Ƿ���ڵ�
	{
		if(P[i] == p1 || P[i] == p2 || P[i] == p3)
		continue;								//�ų��غϵĵ�
		double temp = (P[i] - p1) * Nor;
		if(temp < 0)
		minus = true;
		else if(temp > 0)
		plus = true;
	}
	if(plus && minus)
	return false;								//��ƽ��������ڵ㣬��ʧ���˳�
	
	Flat.clear() ;								//ƽ���ϵĵ㼯Flat��ʼ��
	Flat.PB(p1) ; Flat.PB(p2) ; Flat.PB(p3);	//p1,p2,p3ѹ��Flat
	for(int i = 1; i <= 5; i++)					//��λ��ƽ���ϵĵ�ѹ��Flat
	{
		if(P[i] == p1 || P[i] == p2 || P[i] == p3)
		continue;
		double temp = (P[i] - p1) * Nor;
		if(temp == 0)
		Flat.PB(P[i]);
	}
	return true;
} 
Point Shadow(Point p , Point p1, Point p2, Point p3)	//���ص�p��ƽ�棨��p1,p2,p3���ɣ��ϵ���Ӱ��
{
	Point Nor = (p2 - p1) ^ (p3 - p1);			//����ƽ��ȵķ�����
	Point Nor0 = (Nor * (Nor * (p - p1) )/ Abs(Nor) / Abs(Nor));
	return (p - Nor0);
}
double ShadowHeight(Point p, Point p1, Point p2, Point p3)
												//����p��p��ƽ�棨��p1,p2,p3���ɣ���Ӱ�ľ���
{
	Point Nor = (p2 - p1) ^ (p3 - p1);			//����ƽ��ķ�����
	return fabs(Nor * (p - p1)) / Abs(Nor) ; 
}
bool InTheMiddle(Point A, Point B, Point C, Point P)	//���ؿռ�����AP����AB��AC��ı�־
{
	Point AB = B - A;
	Point AC = C - A;
	Point AP = P - A;
	
	Point v1 = AB ^ AC;
	Point v2 = AB ^ AP;
	
	return (v1 * v2) >= 0 ;
}
bool PointinTriangle(Point A, Point B, Point C, Point P)
												//���ص�p��������ABCn�ڵı�־
{
	return InTheMiddle(A, B, C, P) && InTheMiddle(B, C, A, P) && InTheMiddle(C, A, B, P);
}
bool PointInArea(Point p)					//���ص�p��ƽ������Flat�ڵı�־
{
	for (int i = 0 ; i < Flat.size() ; i ++)	//ö��ƽ������Flat��3��������п������
												//һ�����ֵ�p��ĳ���㹹�ɵ��������ڲ���
												//��˵��p��ƽ������Flat�� 
		for(int j = i + 1 ; j < Flat.size() ; j ++)
			for (int k = j + 1 ; k < Flat.size() ; k ++)
				if(PointinTriangle(Flat[i], Flat[j], Flat[k], p))
				return true;
			return false;					//��p��ƽ������Flat�� 
} 
double DistanceToLine(Point P, Point A, Point B)	//�����P��ֱ��AB�ľ���
{
	if((B - A) * (P - A) > 0.000 && (A - B) * (P - B) > 0.000)
		return sqrt(Sqr(Distance(A, P)) - Sqr((B - A) * (P - A) / Abs((B - A))));
	return fmin(Distance(A, P), Distance(B, P));
}
bool Balance(Point p)						//�����ȶ���־����p�����ⷽ���ƶ�0.2����λ����Ȼ��ƽ������Flat�ڣ�
{
	if( !PointInArea(p))
	   return false;
	for(int i = 0; i < Flat.size() ; i ++)	//ö��Flat�ϵ������߶�
		for(int j = i + 1; j < Flat.size() ; j ++){
			Point P_i = Flat[i], P_j = Flat[j] ;
			int k;
			for( k = 0 ; k < Flat.size() ; k ++ )	//Ѱ��ƽ�����߶�P_i��P_j���һ����P_k
				if(k != i && k != j)
				break;
			Point P_k = Flat[k];
			
			bool plus = false, minus = false;
			Point Nor = (P_j - P_i) ^ (P_k - P_i) ^ (P_j - P_i);	//������
			for (int l = 0 ; l < Flat.size(); l ++ ) if ( l != i && l != j)
			{
				if(Nor * (Flat[l] - P_i ) > 0 )
				plus = true;
				else if(Nor * (Flat[l] - P_i) < 0 )
				minus = true ;
			 } 
			 if (plus && minus)
			 continue;				//P_i,P_j����͹���ϵı�
			 if (DistanceToLine(p , P_i, P_j) < 0.2)
			 						//��p��͹����Ե����С��0.2��˵����pΪԲ�ġ�0.2Ϊ�뾶��Բ����͹�����ڣ����ز��ȶ���־
				return false;
		} 
		return true;
}

void FindHull()						//Ѱ����ά͹���������棬��������оƬ����������������С����
{
	double Min_Answer = 2000.000, Max_Answer = 0.000;
	P[1] = A ; P[2] = B ; P [3] = C ; P [4] = D ; P[5] = E ; P [6] = F;
	
	for(int i = 1 ; i <= 5; i ++)	//ö��5������ѡ3������������
		for(int j = i + 1 ; j <= 5 ; j ++)
			for(int k = j + 1 ; k <= 5 ; k ++){
				if(SameSide(P[i], P[j], P[k]) )	//�����е���ƽ�棨�ɵ�P[i]P[j]P[k]���ɣ���
												//ͬ�࣬��˵����ƽ�����Ϊ��ֽ�ĵ���
				{
					Point G0 = Shadow(G, P[i], P[j], P [k]);	//��������G�ڵ������Ӱ��G0
					if(Balance(G0))						//�������ȶ�״̬��G0�����ⷽ���ƶ�0.2����λ
														//����Ȼ�ڵ��淶Χ�ڣ��������оƬ��F��F����
														//��Ӱ�ľ��룬�������������С����
					{
						Min_Answer = fmin(Min_Answer, ShadowHeight(F, P[i], P[j], P[k]));
						Max_Answer = fmax(Max_Answer, ShadowHeight(F, P[i], P[j], P[k])); 
			        }
			    }
			}
		cout << fixed << Min_Answer << " " << Max_Answer << "\n" ;
														//������������С����
}
int main(){
	cout.precision(5);			//���徫��ΪС�������λ
	int Case = 0, deg ;		//����������ų�ʼ��
	char ss[200];
	while( gets(ss))				//����������ֽ��5�������оƬ���ֱ꣬������0Ϊֹ
	{
		if(strlen(ss) == 1 && ss[0] == '0' )
		break;
		sscanf(ss,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
			&A.X , &A.Y , &A.Z , &B.X , &B.Y , &B.Z ,
			&C.X , &C.Y , &C.Z , &D.X , &D.Y , &D.Z ,
			&E.X , &E.Y , &E.Z , &F.X , &F.Y , &F.Z );
		cout << "Case " << ++ Case << ": " ;	//�����������
		FindCentre();							//������ֽ������
		FindHull();								//��������оƬ����������������С����
	}	 
}

