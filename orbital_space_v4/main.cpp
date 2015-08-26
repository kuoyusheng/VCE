//
//  main.cpp
//  class_trial
//
//  Created by kuoyusheng on 2015/8/12.
//  Copyright (c) 2015年 kuoyusheng. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include </usr/include/Eigen/Eigen> // works ok
#include "null_space.h"
using namespace std;
using Eigen::MatrixXd;
int power(int val, int pow){
    int result=1;
    if(pow==0)
        return result;
    
    for (int i=0; i<(pow); i++) {
        result*=val;
    }
    return result;
}

void digit_combination(int dim,int i,int arr[],int digit)
{
    if (dim / 3 != 0) {
        digit_combination((dim / 3),i+1,arr,digit);
    }
    arr[digit-1-i]=(dim % 3);
}


void tensor_ind(int **tensor,int tensor_dim, int no_of_ver){
    for (int i=0; i<tensor_dim; i++) {
        digit_combination(i,0,tensor[i],no_of_ver);
    }
}


int factorial(int x, int result = 1) {
    if (x == 1||x==0) return result; else return factorial(x - 1, x * result);
}


int** combination_fixed_total(int tot, int groups){
    string s;
    int tot_permu=factorial(tot+groups-1)/(factorial(tot)*factorial(groups-1));
    
    for (int i=0; i<groups-1; i++) {
        s+="0";
    }
    for (int i=0; i<tot; i++) {
        s+="1";
    }
    int **arr=new int*[tot_permu];
    for (int i=0; i<tot_permu; i++) {
        arr[i]=new int[tot+groups-1];
    }
    int i=0;
    do
    {
        //std::cout << s << '\n';
        for (int j=0; j<s.length(); j++) {
            int a=s[j]-48;
            arr[i][j]=a;
        }
        i++;
    }
    while (std::next_permutation(s.begin(), s.end()));
    return arr;
}
void swap(int *fir, int *sec)
{
    int temp = *fir;
    *fir = *sec;
    *sec = temp;
}
static int i=0;
void permutation(int **arr2,int * arr, int curr, int size,  int tot)
{
    if(curr == size-1)
    {
        for(int a=0; a<size; a++)
            arr2[i][a]=arr[a];
        i++;
        if (i==tot) {
            i=0;
        }
    }
    
    else
    {
        for(int i=curr; i<size; i++)
        {
            swap(&arr[curr], &arr[i]);
            permutation(arr2,arr, curr+1, size,tot);
            swap(&arr[curr], &arr[i]);
        }
    }
}
int **permu_fun(int no_of_ver){
    int tot_permu=factorial(no_of_ver);
    int**arr2=new int*[tot_permu];
    for (int i=0; i<tot_permu; i++) {
        arr2[i]=new int[no_of_ver];
    }
    int *arr=new int[no_of_ver];
    for (int i=0; i<no_of_ver; i++) {
        arr[i]=i;
    }
    permutation(arr2, arr, 0, no_of_ver, tot_permu);
    delete [] arr;
    return arr2;
}
class Cluster{
public:
    static Cluster* create(int num, int size);
    
    Cluster();
    Cluster(int no_of_ver);
    Cluster(int no_of_ver, float ** vertex);
    ~Cluster();
    
    void assign_val(Cluster a);
    void map_data(float *clu, int no);
    bool comare();
    float* read_vertex();
    int  read_vertex_no();
    bool compare_transl_dist(Cluster b);
    Cluster* delete_dup(Cluster* a);
    bool compare_sym_dist(Cluster a, float ***linear_trans, float ***transl, int no_of_sym);
    void permu_vertex(int *pi);
    void print_clu();
    void clu_sym_cal(Cluster a,int sym_no, float **linear_trans,float *transl);
    void clu_imp_enum(Cluster , int, int*);
    void imp_gamma_cal();
    void gamma_cal(float ***linear_trans);
    int no_of_sym;
    int isotromy=0;
    int no_of_os=0;
    float *clu_vertex;
    int *permu;
    int no_of_ver;
    bool sym_distinctive;
    bool imp=false;
    int *imp_idx;
    float clu_dist=0;
    float* gamma;
    float* imp_gamma;
    
private:
    static int dim;
    
};
void Cluster::imp_gamma_cal(){
    if (!imp) {
        return;
    }
    int tensor_dim= power(3, no_of_ver);
    int** arr =new int*[tensor_dim];
    for (int i=0; i<tensor_dim; i++) {
        arr[i]=new int[no_of_ver];
    }
    int** permu =new int*[tensor_dim];
    for (int i=0; i<tensor_dim; i++) {
        permu[i]=new int[no_of_ver];
    }
    int npt=no_of_ver;
    int k=0;
    for (int i=0; i<npt; i++) {
        for (int j=i+1; j<npt; j++) {
            if (imp_idx[i]==imp_idx[j]) {
                permu[k][0]=i;
                permu[k][1]=j;
                k++;
            }
        }
    }
//    if(npt==3)
//        cout<<k<<endl;
    imp_gamma=new float[(tensor_dim*k)*tensor_dim];
    for (int i=0; i<tensor_dim*k*tensor_dim; i++) {
        imp_gamma[i]=0;
    }

    float **temp=new float*[tensor_dim];
    for (int i=0; i<tensor_dim; i++) {
        temp[i]=new float[tensor_dim];
    }
    tensor_ind(arr, tensor_dim, 3);
    for (int u=0; u<k; u++) {
        //cout<<"after permu"<<endl;
        for (int i=0; i<tensor_dim;i++ ) {
            swap(arr[i][permu[u][0]],arr[i][permu[u][1]]);
        }
        
        for (int i=0; i<tensor_dim;i++ ) {
            int tmp=0;
            for (int j=0; j<npt; j++) {
                tmp+=(arr[i][j]*power(3, npt-j-1));
            }temp[i][tmp]=1;
        }
        for (int i=0; i<tensor_dim; i++) {
            for (int j=0; j<tensor_dim; j++) {
                imp_gamma[(u*tensor_dim+i)*tensor_dim+j]=temp[i][j];
            }
        }
    }
    delete []temp;
    delete []permu;
    delete []arr;
}
void Cluster::gamma_cal(float ***linear_trans){
    int tensor_dim=power(dim, no_of_ver);
    int **tensor =new int*[tensor_dim];
    for (int i=0; i<tensor_dim; i++) {
        tensor[i]=new int[no_of_ver];
    }
    gamma=new float[tensor_dim*tensor_dim];
    tensor_ind(tensor, tensor_dim, no_of_ver);
    for (int i=0; i<tensor_dim; i++) {
        for (int j=0; j<tensor_dim; j++) {
            float tmp=1;
            for (int k=0; k<no_of_ver; k++) {
                tmp*=linear_trans[no_of_sym][tensor[i][k]][tensor[j][permu[k]]];
                if (tmp==(-0)) {
                    tmp=0;
                }
            }
            gamma[i*tensor_dim+j]=tmp;
        }
    }
    delete []tensor;
}
int Cluster::dim=3;

bool Cluster::compare_sym_dist(Cluster a, float ***linear_tran, float ***transl, int no_of_sym){
    Cluster* tmp=create(no_of_sym, a.no_of_ver);
    for (int i=0;i<no_of_sym; i++) {
        tmp[i].clu_sym_cal(a, i, linear_tran[i], transl[i][0]);
    }
    for (int i=0; i<no_of_sym; i++) {
        if(this->compare_transl_dist(tmp[i]))
            return true;
    }
    return false;
}



void Cluster::clu_imp_enum(Cluster a, int no_of_imp, int *arr){
    imp=true;
    int size=a.no_of_ver+no_of_imp-1;
    clu_dist=a.clu_dist;
    no_of_ver=a.no_of_ver+no_of_imp;
    imp_idx=new int[no_of_ver];
    for (int i=0; i<a.no_of_ver; i++) {
        imp_idx[i]=i;
    }
    
    for (int i=0; i<a.no_of_ver; i++) {
        for (int j=0; j<dim; j++) {
            clu_vertex[i*dim+j]=a.clu_vertex[i*dim+j];
        }
    }
    
    
    for (int i=0,j=0,k=0; i<size;i++) {
        if (arr[i]==0) {
            j++;
        }
        if (arr[i]==1) {
            for (int pos=0; pos<dim; pos++) {
                clu_vertex[(a.no_of_ver+k)*dim+pos]=clu_vertex[j*dim+pos];
            }
            imp_idx[a.no_of_ver+k]=j;
            k++;
        }
    }
}

Cluster::Cluster(){
    clu_vertex=NULL;
    no_of_ver=0;
    no_of_sym=0;
}

void Cluster::assign_val(Cluster a){
    no_of_ver=a.no_of_ver;
    clu_dist=a.clu_dist;
    for (int i=0; i<no_of_ver; i++) {
        for (int j=0; j<dim; j++) {
            clu_vertex[i*dim+j]=a.clu_vertex[i*dim+j];
        }
    }
}


Cluster* Cluster::create(int num, int size)
{
    Cluster* clusters = new Cluster[num];
    for (int i=0; i < num; i++) {
        clusters[i].no_of_ver = size;
        clusters[i].clu_vertex=new float[size*dim];
        clusters[i].permu=new int[size];
    }
    return clusters;
};


void Cluster::print_clu(){
    for (int i=0; i<no_of_ver; i++) {
        for (int j=0; j<dim; j++) {
            cout<<clu_vertex[i*dim+j]<<"\t";
        }cout<<endl;
    }cout<<endl;
}
void Cluster::map_data(float *clu, int no){
    no_of_ver=no;
    clu_vertex=new float[no_of_ver*dim];
    for (int i=0; i<no; i++) {
        for (int j=0; j<dim;j++) {
            clu_vertex[i*dim+j]=clu[i*dim+j];
        }
    }
}

Cluster::~Cluster(){
}
//how to deal with mallac
void Cluster::clu_sym_cal(Cluster a,int sym_idx, float** linear_trans, float *transl){
    no_of_sym=sym_idx;
    isotromy=sym_idx;
    no_of_ver=a.no_of_ver;
    imp_idx=new int[no_of_ver];
    if (a.imp) {
        imp=a.imp;
        for (int i=0; i<a.no_of_ver; i++) {
            imp_idx[i]=a.imp_idx[i];
        }
    }
    permu=new int[no_of_ver];
    for (int i=0; i<no_of_ver; i++) {
        permu[i]=i;
    }
    clu_vertex=new float[no_of_ver*dim];
    Eigen::MatrixXd coord(3,1);
    Eigen::MatrixXd lin_trans(3,3);
    Eigen::MatrixXd tran(3,1);
    Eigen::MatrixXd result(3,1);
    for (int i=0; i<no_of_ver; i++) {
        for (int j=0; j<dim; j++) {
            coord(j,0)=a.clu_vertex[i*dim+j];
            tran(j,0)=transl[j];
            for (int k=0; k<dim; k++) {
                lin_trans(j,k)=linear_trans[j][k];
            }
        }
        result=lin_trans*coord+tran;
        for (int j=0; j<dim; j++) {
            clu_vertex[i*dim+j]=result(j,0);
        }
    }
}

Cluster::Cluster(int no_ver, float** vertex){
    no_of_ver=no_ver;
    clu_vertex=new float [no_ver*3];
    for (int i=0; i<no_of_ver; i++) {
        for (int j=0; j<dim;j++) {
            clu_vertex[(dim*i)+j]=vertex[i][j];
        }
    }
}

Cluster::Cluster(int no){
    no_of_ver=no;
    clu_vertex=new float[no*dim];
}


float* Cluster::read_vertex(){
    return clu_vertex;
}

bool Cluster::compare_transl_dist(Cluster b){
    int tot_permu=factorial(no_of_ver);
    int a=b.no_of_ver;
    int **pi=permu_fun(a);
    
    
    float** distance=new float*[no_of_ver];
    for (int i=0;i<no_of_ver ; i++) {
        distance[i]=new float[3];
    }
    bool t=true;
    for(int idx=0;idx<tot_permu;idx++,t=1){
        for (int i=0; i<no_of_ver; i++) {
            for (int j=0; j<dim; j++) {
                distance[i][j]=clu_vertex[dim*pi[idx][i]+j]-b.clu_vertex[dim*i+j];
                int test=distance[i][j];
                if (distance[i][j]!=test)
                    t=0;
            }
        }
        if(!t){
            continue;
        }
        for (int i=0; i<no_of_ver-1; i++) {
            for (int j=0; j<3; j++) {
                if (distance[i][j]!=distance[i+1][j]){
                    t=0;
                    break;
                }
            }
            if (!t) {
                break;
            }
        }
        if(t){
            for(int i=0;i<no_of_ver;i++)
                permu[i]=pi[idx][i];
            delete[] distance;
            return t;
        }
    }
    
    delete []distance;
    return t=0;
    
    return t;
};

void find_orbital_space(Cluster* clu, int no_of_ver, int no_of_rep, int no_of_sym){
    
    for(int idx=0;idx<no_of_rep;idx++){
        for (int i=0; i<no_of_sym-1; i++) {
            for (int j=i+1; j<no_of_sym; j++) {
                if(clu[idx*no_of_sym+j].compare_transl_dist(clu[idx*no_of_sym+i])){
                    if (clu[idx*no_of_sym+i].no_of_os) {
                        clu[idx*no_of_sym+j].isotromy=clu[idx*no_of_sym+i].isotromy;
                    }
                    else
                        clu[idx*no_of_sym+j].isotromy=clu[idx*no_of_sym+i].no_of_sym;
                    clu[idx*no_of_sym+j].no_of_os++;
                }
            }
        }
    }
}

Cluster* imp_clu_enu(int no_of_imp,int no_of_ver,int no_of_rep, Cluster* clu_rep,int & tot){
    int **com=combination_fixed_total(no_of_imp,no_of_ver);
    tot=factorial(no_of_imp+no_of_ver-1)/(factorial(no_of_imp)*factorial(no_of_ver-1));
    
    Cluster* imp_pair=Cluster::create(no_of_rep*tot, no_of_ver+no_of_imp);
    
    for (int idx=0; idx<no_of_rep; idx++) {
        for (int i=0; i<tot; i++) {
            imp_pair[idx*tot+i].clu_imp_enum(clu_rep[idx], no_of_imp, com[i]);
        }
    }
    delete []com;
    return imp_pair;
}
Cluster* delete_dup(Cluster* a, int no_of_rep,int& no_of_dist, float ***linear_tran, float ***transl){
    int count=0;
    for (int i=0; i<no_of_rep; i++) {
        for (int j=i+1; j<no_of_rep; j++) {
            if(a[j].compare_sym_dist(a[i], linear_tran, transl, 48)){
                a[j].no_of_os++;
            }
        }
        if (!a[i].no_of_os) {
            count++;
        }
    }
    Cluster* tmp=Cluster::create(count, a[0].no_of_ver);
    int k=0;
    for (int i=0; i<no_of_rep; i++) {
        if (!a[i].no_of_os) {
            tmp[k]=a[i];
            k++;
        }
    }
    no_of_dist=count;
    return tmp;
}

void cluster_rep(int** arr_rep_no,Cluster **clusters_rep,int no_of_ver,float***linear_tran, float***transl)
{
    for(int ver=1;ver<=no_of_ver;ver++){
        int prop_no=arr_rep_no[ver][0];
        int imp_prop_no=prop_no;
        int tot=0;
        for (int i=1; i<ver; i++) {
            cout<<"shit"<<endl;
            int no_of_imp=i;
            int no_ver=ver-no_of_imp;
            int multip=0;
            int no_of_dist=0;
            tot+=arr_rep_no[no_ver][0];
            
            Cluster* imp=imp_clu_enu(no_of_imp, no_ver, tot,clusters_rep[no_ver], multip);
            Cluster* imp_clu_rep=delete_dup(imp, tot*multip,no_of_dist,linear_tran, transl);
            arr_rep_no[ver][i]=no_of_dist;
            imp_prop_no+=no_of_dist;
            for (int k=prop_no,j=0; k<imp_prop_no; k++,j++) {
                clusters_rep[ver][k]=imp_clu_rep[j];
            }
            prop_no=imp_prop_no;
        }
    }
}
Cluster** clu_symmety(int no_of_ver,int no_of_sym, Cluster** clusters_rep, int **arr_rep_no, float***linear_tran, float***transl){
    int clu_num[4];
    
    for (int i=1; i<=no_of_ver; i++) {
        for (int j=0; j<4; j++) {
            clu_num[i]+=arr_rep_no[i][j];
        }
    }
    Cluster** clu_sym = new Cluster*[no_of_ver];
    for (int i=1; i<=no_of_ver; i++) {
        clu_sym[i]=Cluster::create(5000,i);
    }

    for (int ver=1; ver<=no_of_ver; ver++) {
        for (int no=0; no<clu_num[ver]; no++) {
            for (int sym=0; sym<no_of_sym; sym++) {
                clu_sym[ver][no*no_of_sym+sym].clu_sym_cal(clusters_rep[ver][no],sym, linear_tran[sym],transl[sym][0]);
            }
        }
    }
    return clu_sym;
}
void isotromy_gamma(Cluster* clu_sym, int* arr, int no_of_os){
    int num=clu_sym[arr[0]].no_of_ver;
    int tensor_dim=power(3, num);
    int tot=(no_of_os*tensor_dim)*tensor_dim;
    float* temp=new float[tot];
    for (int i=0,u=0; i<no_of_os; i++) {
        for(int j=0;j<tensor_dim;j++){
            for(int k=0;k<tensor_dim;k++,u++){
                temp[u]=clu_sym[arr[i]].gamma[j*tensor_dim+k];
                if (j==k) {
                    temp[u]--;
                }
            }
        }
    }
    clu_sym[arr[0]].gamma=new float[tot];
    for (int i=0; i<tot; i++) {
        clu_sym[arr[0]].gamma[i]=temp[i];
    }
    delete []temp;
}


int main() {
    ifstream sym("/Users/kuoyusheng/Desktop/sym1.out");
    ifstream clu("/Users/kuoyusheng/Desktop/clusters1.out");
    float a=0;
    int no_of_sym=0;
    sym>>a;
    no_of_sym=a;
    //create the matrix to read in the symmetry group
    //one linear transformation group, one translation group
    float ***linear_tran=new float**[no_of_sym];
    float ***transl=new float**[no_of_sym];
    for(int i=0;i<a;i++){
        linear_tran[i]=new float*[3];
        for(int j=0;j<3;j++){
            linear_tran[i][j]=new float[3];
            
        }
    }
    
    for(int i=0;i<a;i++){
        transl[i]=new float*[1];
        for(int j=0;j<1;j++){
            transl[i][j]=new float[3];
            
        }
    }
    for(int u=0;u<no_of_sym;u++){
        for(int i=0;i<3;i++){
            for(int j=0;j<3;j++){
                sym>>a;
                linear_tran[u][i][j]=a;
            }
        }
        for(int i=0;i<3;i++){
            sym>>a;
            transl[u][0][i]=a;
        }
        
    }
    
    
    //create 1,2,3 clusters 3-d matrix[ind][row][column]
    float **single=new float*[2];
    for(int i=0;i<2;i++){
        single[i]=new float[3];
    }
    
    float **pair=new float*[24];
    for(int i=0;i<24;i++){
        pair[i]=new float[6];
    }
    float **triplet=new float*[6];
    for(int i=0;i<6;i++){
        triplet[i]=new float[9];
    }
    float **data=new float*[100];
    for (int i=0; i<100; i++) {
        data[i]=new float[3];
    }
    string cluster;
    int single_ind=0;
    int pair_ind=0;
    int tri_ind=0;
    float distance[100];
    int k=0;
//    float dist=1.5;
//    float dist_tri=1.0;
    //acquire rep. cluster
    while(getline(clu,cluster)){
        for(int i=0;i<2;i++){
            getline(clu,cluster);
            istringstream s(cluster);
            if (i==0) {
                float a;
                s>>a;
                distance[k]=a;
                k++;
            }
        }
        if (cluster=="0"){
            getline(clu,cluster);
        }
        
        else if(cluster=="1"){
            getline(clu,cluster);
            istringstream s(cluster);
            for(int i=0; i<3;i++){
                float a;
                s>>a;
                single[single_ind][i]=a;
            }
            single_ind++;
            getline(clu,cluster);
        }
        else if (cluster=="2"){
            for(int i=0;i<2;i++){
                getline(clu,cluster);
                istringstream s(cluster);
                for(int j=0; j<3;j++){
                    float a;
                    s>>a;
                    pair[pair_ind][i*3+j]=a;
                }
            }
            pair_ind++;
            getline(clu,cluster);
        }
        else if (cluster=="3"){
            for(int i=0;i<3;i++){
                getline(clu,cluster);
                istringstream s(cluster);
                for(int j=0; j<3;j++){
                    float a;
                    s>>a;
                    triplet[tri_ind][i*3+j]=a;
                }
            }
            tri_ind++;
            getline(clu,cluster);
            //cout<<endl;
        }
    }
    int** arr_rep_no=new int*[4];
    for (int i=0; i<4; i++) {
        arr_rep_no[i]=new int[4];
    }
    //first ind rep. no_of_cluster, second ind rep no_of_imp;
    arr_rep_no[0][0]=1;
    arr_rep_no[1][0]=single_ind;
    arr_rep_no[2][0]=pair_ind;
    arr_rep_no[3][0]=tri_ind;
    
    Cluster** clusters_rep=new Cluster*[4];
    for (int i=0; i<4; i++) {
        clusters_rep[i]=new Cluster[100];
    }
    //int no_of_imp=0;
    //Cluster* pro_clu_rep= new Cluster[index];
    for (int i=0; i<arr_rep_no[2][0]; i++) {
        clusters_rep[2][i].map_data(pair[i], 2);
        clusters_rep[2][i].clu_dist=distance[1+single_ind+i];
    }
    for (int i=0; i<arr_rep_no[1][0]; i++) {
        clusters_rep[1][i].map_data(single[i], 1);
    }
    for(int i=0;i<arr_rep_no[3][0];i++){
        clusters_rep[3][i].map_data(triplet[i], 3);
        clusters_rep[3][i].clu_dist=distance[1+single_ind+pair_ind+i];
    }
//    
//    for (int i=0; i<arr_rep_no[3][0]; i++) {
//            clusters_rep[3][i].print_clu();
//    }
//    for (int i=0; i<arr_rep_no[2][0]; i++) {
//            clusters_rep[2][i].print_clu();
//    }
//
    int no_of_ver=3;

    cout<<"Enumerate imp clusters"<<endl;
    cluster_rep(arr_rep_no, clusters_rep, 3, linear_tran, transl);
//    int prop_no=arr_rep_no[no_of_ver][0];
//    int tot=0;
//    int imp_prop_no=prop_no;
//    for (int i=1; i<no_of_ver; i++) {
//        int no_of_imp=i;
//        int no_ver=no_of_ver-no_of_imp;
//        int multip=0;
//        int no_of_dist=arr_rep_no[no_ver][0];
//        tot+=arr_rep_no[no_ver][0];
////        for (int j=0; j<4; j++) {
////            tot+=arr_rep_no[no_ver][j];
////        }
//        Cluster* imp=imp_clu_enu(no_of_imp, no_ver, tot, clusters_rep[no_ver], multip);
//        Cluster* imp_clu_rep=delete_dup(imp, tot*multip,no_of_dist,linear_tran, transl);
//        arr_rep_no[no_of_ver][i]=no_of_dist;
//        imp_prop_no+=no_of_dist;
//        for (int k=prop_no,j=0; k<imp_prop_no; k++,j++) {
//            clusters_rep[no_of_ver][k]=imp_clu_rep[j];
//        }
//        prop_no=imp_prop_no;
//    }
//    for (int i=0; i<no_of_ver; i++) {
//        cout<<arr_rep_no[no_of_ver][i]<<"\t";
//    }
    int clu_num[4];
    
    for (int i=1; i<=3; i++) {
        for (int j=0; j<4; j++) {
            clu_num[i]+=arr_rep_no[i][j];
        }
    }
    //cout<<imp_prop;
//    for (int i=0; i<clu_num[3]; i++) {
//        if(clusters_rep[no_of_ver][i].clu_dist<dist_tri)
//            clusters_rep[no_of_ver][i].print_clu();
//    }
//    for (int i=0; i<clu_num[3]; i++) {
//        if (clusters_rep[no_of_ver][i].imp) {
//            for (int j=0; j<no_of_ver; j++) {
//                cout<<clusters_rep[no_of_ver][i].imp_idx[j]<<"\t";
//            }cout<<endl;
//        }
//    }


//    for (int i=0; i<no_of_dist; i++) {
//        imp_clu_rep[i].print_clu();
//    }
    Cluster** clu_sym = new Cluster*[no_of_ver];
    for (int i=1; i<=no_of_ver; i++) {
        clu_sym[i]=Cluster::create(5000,i);
    }

    for (int ver=1; ver<=no_of_ver; ver++) {
        for (int no=0; no<clu_num[ver]; no++) {
            for (int sym=0; sym<no_of_sym; sym++) {
                clu_sym[ver][no*no_of_sym+sym].clu_sym_cal(clusters_rep[ver][no],sym, linear_tran[sym],transl[sym][0]);
            }
        }
    }
    
//    cout<<"Calculate"<<endl;
//    
//        for (int i=0; i<no_of_sym; i++) {
//            clu_sym[i].clu_sym_cal(clusters_rep[2][0], i, linear_tran[i], transl[i][0]);
//        }
//
    cout<<"Find Orbital Space"<<endl;
    for (int ver=1; ver<=no_of_ver; ver++) {
        find_orbital_space(clu_sym[ver], ver, clu_num[ver], no_of_sym);
    }
    
    
//
//
    int array[4]={};
    for (int i=0,u=0; i<no_of_sym; i++) {
        if (clu_sym[3][i].isotromy==10) {
            array[u]=i;
            u++;
        }
    }cout<<endl;
    cout<<clu_num[3]<<endl;
    cout<<"Find gamma conversion"<<endl;
    for (int ver=1; ver<=no_of_ver; ver++) {
        for (int i=0; i<clu_num[ver]*48; i++) {
            clu_sym[ver][i].gamma_cal(linear_tran);
            clu_sym[ver][i].imp_gamma_cal();
        }
    }
    isotromy_gamma(clu_sym[3], array, 4);
//    for (int i=0; i<4*27; i++) {
//        for (int j=0; j<27; j++) {
//            cout<<clu_sym[3][10].gamma[i*27+j]<<"\t";
//        }cout<<endl;
//    }
    
//    for (int i=0; i<27; i++) {
//        for (int j=0; j<27; j++) {
//            cout<<clu_sym[3][12].gamma[27*i+j]<<"\t";
//        }cout<<endl;
//    }
//    cout<<endl;
    double* b_out=new double[27*27*4];
    int col=nullspace(27, 27*4, clu_sym[3][10].gamma, b_out);
    cout<<col<<endl;
    for (int i=0; i<27*27*4; i++) {
        for (int j=0;j<col; j++) {
            cout<<b_out[i*col+j]<<"\t";
        }cout<<endl;
    }
////
//    for (int i=0; i<clu_num[3]*48; i++) {
//        if(clu_sym[3][i].imp){
//            for (int k=0; k<27; k++) {
//                for (int j=0; j<27; j++) {
//                    cout<<clu_sym[3][i].imp_gamma[k*27+j]<<"\t";
//                }cout<<endl;
//            }
//            break;
//        }
//    }
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    return 0;
}
