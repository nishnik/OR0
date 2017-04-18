				// TEST2.CPP

#include "gauss_elimination.cpp"
#include <iostream>
#include <iomanip> 
#include <vector>

using namespace std;

int gauss_elimination( int n,float **A, float *b, float *sol);

vector< vector<float> > solutions;


void optimize ()
{
    int var;
    float max=-10000000000, temp;
    cout << "Enter the no of variables in the objective equation :";
    cin>>var;
    float eqn[var];
    cout << "Enter the Equation to Maximise : "<<endl;
    for (int i=0; i<var; i++)
        cin>>eqn[i];
    bool bfs;
    for (int i=0; i<solutions.size(); i++)
    {
        bfs=true;
        // degenerate=0;
        temp=0;
        for (int j=0; j<solutions[i].size(); j++){
            if (solutions[i][j]<0)
                bfs=false;
            temp=temp+eqn[j]*solutions[i][j];
        }

        if (bfs == true && temp>max)
        {
            max=temp;
        }
    }

    cout<<endl<<"MAX Value is "<<max<<endl; 
}


void print_sol (int n, int m, float **A, float *b)
{

    cout<<"-----------------SYSTEM OF EQUATIONS ARE-----------------------------"<<endl;
    int i, j;
    for (i=0; i<m; i++){
        for (j=0; j<n; j++)
            cout<<setw(4)<<setprecision(3)<<A[i][j];
        cout<<" |  x"<<i+1<<" = "<<b[i];
        cout<<endl;
    }

    cout<<"-----------------SOLUTIONS-----------------------------"<<endl;
    for (int i=0; i< n; i++)
        cout<<setw(5)<<"x"<<i+1;
    cout<<endl;

    bool bfs;
    int degenerate;
    for (int i=0; i<solutions.size(); i++)
    {
        bfs=true;
        degenerate=0;
        for (int j=0; j<solutions[i].size(); j++){
            if (solutions[i][j]<0)
                bfs=false;
            if (solutions[i][j]==0)
                degenerate++;
            cout<<setw(6)<<solutions[i][j];
        }

        if (bfs == true )
            cout<<"   B.F.S ,, Extreme Point ,,";
        if (degenerate >(n-m))
            cout<<"  Degenerate Basic Solution";
        else
            cout<<"   Non-Degenrate Basic Solution";

        cout<<endl; 
    }
}

void combinationUtil(float **A, int arr[], int data[], int start, int end,
                     int index, int r, float* b)
{
    int i, j;
    // Current combination is ready to be printed, print it
    if (index == r)
    {
    	float **temp_A;
        temp_A= new float* [r]; 

        for (i=0; i<r; i++)
            temp_A[i]= new float [r];

	    float temp_b[r];
        float* sol = new float [r];

	    for (i=0; i<r; i++)
	    	temp_b[i]=b[i];

        for (j=0; j<r; j++)
        {
        	for(i=0; i<r; i++)
        	{
        		temp_A[j][i]=A[j][data[i]];
        	}
        }

        int res = gauss_elimination(r, temp_A, temp_b, sol);
        if (res == 1)
        {
            vector <float> full_sol;
            cout<<"FULL SOL : ";
            j=0;
            for (i=0; i<=end; i++)
            {
                if (i==data[j]){
                    full_sol.push_back(sol[j]);
                    j++;
                }
                else
                    full_sol.push_back(0);
            }

            solutions.push_back(full_sol);

            for (i=0; i<=end; i++)
                cout<<setw(4)<<solutions[solutions.size()-1][i];
            cout<<endl;
        }
        return;
    }
 
    for (i=start; i<=end && end-i+1 >= r-index; i++)
    {
        data[index] = arr[i];
        combinationUtil(A, arr, data, i+1, end, index+1, r, b);
    }
}

int main()
{
	int n, m, i, j, *arr;
	float **A, *b;

	cout<<"Enter No of Equations : ";
	cin>>m;

	cout<<"Enter No of Unknwons : ";
	cin>>n;

	// A = (float **)malloc(m * sizeof(float));
    A= new float* [m]; 

    for (i=0; i<m; i++)
         // A[i] = (float *)malloc(n * sizeof(float));
        A[i]= new float [n];

    // arr=(int *)malloc(n*sizeof(int));
    arr = new int [n];

	cout<<"Enter the matrix A : \n";
	for (i=0; i<m; i++){
		for (j=0; j<n; j++)
			cin>>A[i][j];
	}

	cout<<"Enter the column vector 'b' : \n";
	// b= (float *)malloc(m*sizeof(float));
    b=new float [m];
	for (i=0; i<m; i++)
		cin>>b[i];

	for (i=0; i<n; i++)
		arr[i]=i;

	int data[m]; 
	combinationUtil(A , arr, data, 0,n-1, 0, m, b);

    print_sol(n,m, A, b);

    optimize();
}
