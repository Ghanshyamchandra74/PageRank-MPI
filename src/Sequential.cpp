#include <bits/stdc++.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <time.h>
using namespace std;
 
void addEdge(vector< vector<pair<int,double> > > &adj, int u,
                                     int v, double wt)
{
    adj[u].push_back(make_pair(v, wt));
    adj[v].push_back(make_pair(u, wt));
}

double Randomdouble(double a, double b) {
    double random = ((double) rand()) / (double) RAND_MAX;
    double diff = b - a;
    double r = random * diff;
    return a + r;
}

double L2Norm(double *x, double *x_,int N){
    double L2Norm = 0.0;
    for (int i = 0; i < N; i++)
    {
        L2Norm += pow((x[i]-x_[i]),2);
    }
    
    return sqrt(L2Norm);
}
 

int main()
{
    int n_edges = 4194296; //No. of edges
    int nodes = 1048576;
    double d = 0.85;
    double epsilon = 1e-03;
    
    ifstream in_file("data/large_graph.dat");
    int *u = new int[n_edges];
    int *v = new int[n_edges];
    int *w = new int[n_edges];
    int i=0;
    
    string null,tmp,tmp_a , u_tmp,v_tmp, w_tmp;
    while (!in_file.eof())
    {
        in_file >> null;// >> u_tmp >> v_tmp >> w_tmp ;
       if (null=="a")
       {
           if (i==0)
           {
               //Skipping Headers
           }else if (i<=n_edges)
           {
                in_file >> u_tmp >> v_tmp >> w_tmp ;
                u[i-1] = stoi(u_tmp);
                v[i-1] = stoi(v_tmp);
                w[i-1] = stoi(w_tmp);
                //cout << u_tmp << "," << v_tmp << "," << w_tmp << endl;
           }
           i++;
       }  
    }
    in_file.close();
    
    vector< vector<pair<int,double> > > adj(nodes,vector<pair<int,double> >(0));
    for (size_t i = 0; i < n_edges; i++)
    {
        addEdge(adj, v[i], u[i], w[i]);
    }
    clock_t start = clock();
    for (size_t i = 0; i < nodes; i++)
    {
        double sum = 0.0;
        for (auto it = adj[i+1].begin(); it!=adj[i+1].end(); it++)
        {
            sum += it->second;
        }
        for (auto it = adj[i+1].begin(); it!=adj[i+1].end(); it++)
        {
            if (sum!=0)
            {
                it->second = (it->second)/sum;
            }
        }
    }
    for (size_t i = 0; i < nodes; i++)
    {
        for (auto it = adj[i+1].begin(); it!=adj[i+1].end(); it++)
        {
            it->second = double(((it->second)*d + ((1-d)/double(nodes))));
        }
    }
    double* x = new double[nodes]; *x = 0.0;
    double* x_ = new double[nodes]; *x_ = 0.0;
    double min  = 0.000, max = 0.9999;
    double norm1_x = 0.0;
    for (size_t i = 0; i < nodes; i++)
    {
        x[i] =  Randomdouble(min,max);
        norm1_x += x[i];
    }
    for (size_t i = 0; i < nodes; i++)
    {
         x[i] = x[i]/norm1_x;
    }
    double err = 1e+5;
    int z = 0;
    while (err>=epsilon)
    {
        for (int i = 0; i < nodes; i++)
        {
            x_[i] = x[i];
            x[i] = 0.0;
        }
        for (size_t i = 0; i < nodes; i++)
        {
            for (auto it = adj[i+1].begin(); it!=adj[i+1].end(); it++)
            {
                x[it->first] +=  double((it->second)*x_[i]);
            }
        }
        /*Normalising PageRank after Each iteration */
        double norm1_x = 0.0;
        for (size_t i = 0; i < n_edges; i++)
        {
            norm1_x += x[i];
        }
        for (size_t i = 1; i < nodes; i++)
        {
            x[i] = x[i]/norm1_x;
        }
        err = L2Norm(x,x_,nodes);
        cout << "L2 Error : " << err << endl;
    }
    clock_t end = clock();
    float seconds = (float)(end - start) / CLOCKS_PER_SEC;
    ofstream time_file("data/Sequential_time.csv");
    time_file << "t" << endl;
    time_file << seconds;
    /* Calculating Vector Sum */
    double sum = 0.0;
        for (size_t i = 1; i < nodes; i++)
        {
            sum+=x[i];
        }
        /* Sum of All Vector must be Equal to 1 */
        cout << "Vector Sum check : " << (sum <= 1 && sum >= 0.99) << endl;
    /* Writing Results */
    // x is final PageRank
    ofstream out_file("data/Sequential_PageRank.csv");
    for (size_t i = 0; i < nodes; i++)
    {
        out_file << x[i] << endl;
    }
    cout << " Finished Writing PageRank to a file :" << endl;
    return 0;
}