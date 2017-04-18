

#include <iostream>
#include <iomanip> 
#include <vector>
#include <cmath>

using namespace std;


struct tableau
{
	float value;
	int isBasic;
};

void print(int m, int n, float **A)
{
	int i, j;
	cout<<"----------------"<<endl;
	for (i=0; i<m; i++){
		for (j=0; j<(n); j++)
			cout<<setw(10)<<setprecision(5)<<A[i][j];
		cout<<endl;
	}
	// cout<<"----------------"<<endl;
}

void printT(int m, int n, tableau **A)
{
	int i, j;
	cout<<"----------------"<<endl;
	for (i=0; i<m; i++){
		for (j=0; j<(n); j++)
			cout<<setw(10)<<setprecision(5)<<A[i][j].value;
		cout<<endl;
	}
	// cout<<"----------------"<<endl;
}

void solution(tableau **T, float **cost, int nrow, int ncol)
{
	float value;
	for (int i=0; i<nrow; i++)
		for (int j=0; j<ncol; j++)
			if (T[i][j].isBasic!=0)
				value=value+cost[i][j]*T[i][j].value;

	cout<<"Value is "<<value;
}

void northWest(int nrow,int ncol,tableau** T)
{
	int curr_x, curr_y;
	curr_x=0;
	curr_y=0;

	while(1){
		T[curr_x][curr_y].isBasic=1;
		if (curr_x>nrow-1 || curr_y>ncol-1)
			break;
		if (T[nrow][curr_y].value<=T[curr_x][ncol].value){
			T[curr_x][curr_y].value=T[nrow][curr_y].value;
			T[curr_x][ncol].value=T[curr_x][ncol].value-T[nrow][curr_y].value;
			T[nrow][curr_y].value=0;
			curr_y++;
		}
		else
		{
			T[curr_x][curr_y].value=T[curr_x][ncol].value;
			T[nrow][curr_y].value=T[nrow][curr_y].value- T[curr_x][ncol].value;
			T[curr_x][ncol].value=0;
			curr_x++;
		}

		
	}
	printT(nrow+1, ncol+1, T);
}

void iterateop(tableau **T, float **cost, int s, int d)
{
	int i, j;
	for (i=0; i<s; i++)
		T[i][d+1].value=-1111;
	for (j=0; j<d; j++)
		T[s+1][j].value=-1111;
	// Multipliers calculation
	T[0][d+1].value=0;
	for (i=0; i<s; i++)
	{
		for (j=0; j<d; j++)
		{
			if (T[i][j].isBasic == 1)
			{
				if (T[i][d+1].value==-1111)
					T[i][d+1].value=cost[i][j]-T[s+1][j].value;
				if (T[s+1][j].value==-1111)
					T[s+1][j].value=cost[i][j]-T[i][d+1].value;
			}
		}
	}
	cout<<"dsdsdd"<<endl;
	printT(s+2, d+2, T);
}

int main()
{
	int s, d, i, j;
	float **cost;
	tableau **T;
	cout<<"Enter no of source and destination"<<endl;
	cin>>s>>d;

	cost = new float* [s+1]; 

    for (i=0; i<s+1; i++){
        cost[i]= new float [d+1];
    }

    cout<<"Enter the cost matrix with supply and demand"<<endl;
	for (i=0; i<s+1; i++){
		for (j=0; j<d+1; j++){
			cin>>cost[i][j];
		}
	}

    print(s+1, d+1, cost);

    T=new tableau* [s+1+1];
    for (i=0; i<s+1+1; i++){
        T[i]= new tableau [d+1+1];
    }

    for (i=0; i<s; i++)
    	T[i][d].value=cost[i][d];
    for (i=0; i<d; i++)
    	T[s][i].value=cost[s][i];

    printT(s+1, d+1, T);
    northWest(s, d, T);
    solution(T, cost, s, d);
    cout<<"&&&&&&"<<endl;
    iterateop(T, cost, s, d);



}