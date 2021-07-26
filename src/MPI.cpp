#include <bits/stdc++.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <mpi.h>
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
 

int main(int argc, char* argv[])
{
    int n_edges = 4194296; //No. of edges
    int nodes = 1048576;
    double d = 0.85;
    double epsilon = 5e-04;
    double t1,t2;
    /*String for Writing File Names */
    string str_time,file_name;
    
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
    /* Taking Values to Vector */
    vector< vector<pair<int,double> > > adj(nodes,vector<pair<int,double> >(0));
    for (size_t i = 0; i < n_edges; i++)
    {
        addEdge(adj, v[i], u[i], w[i]);
    }
    /*Begin MPI Region */
    int proc_count,my_proc_id;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD,&proc_count);
    MPI_Comm_rank(MPI_COMM_WORLD,&my_proc_id);
    file_name = "data/Page_rank_"+to_string(proc_count)+"_MPI.csv";
    ofstream page_rank_file(file_name);
    str_time = "data/Time_Exec_"+to_string(proc_count)+"_MPI.csv";
    ofstream time_file(str_time);
    MPI_Barrier(MPI_COMM_WORLD);
    time_file << "t" << endl; // Pushing Header
    t1 = MPI_Wtime();
    /* Calculating Chunk Size */
    int chunk_size_nodes = nodes/proc_count;
    int chunk_size_edges = n_edges/proc_count;
    double* time = new double[proc_count];
    /* Global Arrays */
    double* x = new double[nodes]; *x = 0.0;
    double* x_ = new double[nodes]; *x_ = 0.0;
    /* Local Arrays */
    double* x_loc = new double[chunk_size_nodes];
    double* x_loc_ = new double[chunk_size_nodes];
    //std::cout << "debug 1" << std::endl;
    
    for (size_t i = my_proc_id*chunk_size_nodes; i < (my_proc_id*chunk_size_nodes+chunk_size_nodes); i++)
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
    //std::cout << "debug 2" << std::endl;
    for (size_t i = my_proc_id*chunk_size_nodes; i < (my_proc_id*chunk_size_nodes+chunk_size_nodes); i++)
    {
        for (auto it = adj[i+1].begin(); it!=adj[i+1].end(); it++)
        {
            it->second = double(((it->second)*d + ((1-d)/double(nodes))));
        }
    }
    //std::cout << "debug 3" << std::endl;
    double min  = 0.000, max = 0.9999;
    double norm1_x = 0.0;
    double norm1_x_final = 0.0;
    for (size_t i = 0; i < chunk_size_nodes; i++)
    {
        x_loc[i] = Randomdouble(min,max);
        norm1_x += x_loc[i];
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allreduce(&norm1_x,&norm1_x_final,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    for (size_t i = 0; i < chunk_size_nodes; i++)
    {
        x_loc[i] = x_loc[i]/(norm1_x_final);
    }
    MPI_Allgather(x_loc,chunk_size_nodes,MPI_DOUBLE,x,chunk_size_nodes,MPI_DOUBLE,MPI_COMM_WORLD);
    //std::cout << "debug 4" << std::endl;
    //MPI_Allgather(x_loc,chunk_size_nodes,MPI_DOUBLE,x,chunk_size_nodes,MPI_DOUBLE,MPI_COMM_WORLD);
    if (my_proc_id==0)
    {
        double sum = 0.0;
        for (size_t i = 0; i < nodes; i++)
        {
            sum+=x[i];
        }
        std::cout << "Sum :" << sum << std::endl;
    }
    
    //std::cout << "debug 5" << std::endl;
    double err = 1e+5;
    int z = 0;
    //while (err>epsilon)
    while(z<50)
    {
        for (size_t i = 0; i < chunk_size_nodes; i++)
        {
            x_loc_[i] = x_loc[i];
            x_loc[i] = 0.0;
        }
        MPI_Allgather(x_loc_,chunk_size_nodes,MPI_DOUBLE,x_,chunk_size_nodes,MPI_DOUBLE,MPI_COMM_WORLD);
        //std::cout << "debug 6" << std::endl;
        for (size_t i = my_proc_id*chunk_size_nodes; i < (my_proc_id*chunk_size_nodes+chunk_size_nodes); i++)
        {
            for (auto it = adj[i+1].begin(); it!=adj[i+1].end(); it++)
            {
                x[it->first] +=  double((it->second)*x_[i]);
            }
        }
        //std::cout << "debug 7" << std::endl;
        /*Normalising PageRank after Each iteration */
        double norm1_x1 = 0.0;
        double norm1_x1_final = 0.0;
        for (int i = my_proc_id*chunk_size_nodes; i < (my_proc_id*chunk_size_nodes + chunk_size_nodes); i++)
        {
            x_loc[i-my_proc_id*chunk_size_nodes] = x[i];
            norm1_x1 += x[i];
        }
        //cout << "Rank : " << my_proc_id << "Norm : " << norm1_x1 << endl;
        //std::cout << "debug 8" << endl; 
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Allreduce(&norm1_x1,&norm1_x1_final,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        for (int i = 0; i < chunk_size_nodes; i++)
        {
            x_loc[i] = x_loc[i]/norm1_x1_final;
        }
        MPI_Allgather(x_loc,chunk_size_nodes,MPI_DOUBLE,x,chunk_size_nodes,MPI_DOUBLE,MPI_COMM_WORLD);
        // //cout << "Rank : " << my_proc_id << "Norm_Final : " << norm1_x1_final << endl;
        // //std::cout << "debug 8" << std::endl;
        double l2norm = 0.0;
        double l2norm_final = 0.0;
        for (int i = 0; i <  chunk_size_nodes; i++)
        {
            l2norm+=(pow((x_loc[i]-x_loc_[i]),2));
        }
        //std::cout << "debug 9" << endl;
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Reduce(&l2norm,&l2norm_final,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
        // //MPI_Scan(&l2norm,&l2norm_final,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        if (my_proc_id==0)
        {
            err = sqrt(l2norm_final);
            cout << " <=============== L2 Error : "<< err << endl;
        }
        //std::cout << "debug 10" << std::endl;
        //MPI_Barrier(MPI_COMM_WORLD);
        z++;   
    }
    t2 = MPI_Wtime();
    double diff = (t2 -t1),diff_final;
    MPI_Reduce(&diff,&diff_final,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    // time[my_proc_id] = t2 - t1;
    if (my_proc_id==0)
        {
            double sum = 0.0;
            for (size_t i = 1; i < nodes; i++)
            {
                sum+=x[i];
            }
            /* Sum of All Vector must be Equal to 1 */
            cout << "Vector Sum check : " << (sum <= 1 && sum >= 0.99) << endl;
            /* Writing Results */
            time_file << diff_final/proc_count << endl;
            for (size_t i = 0; i < nodes; i++)
            {
                page_rank_file << x[i] << endl;
            }
            cout << " Finished Writing PageRank to a file :" << endl;
            time_file.close();
            page_rank_file.close();
        };
    MPI_Finalize();
    return 0;
}