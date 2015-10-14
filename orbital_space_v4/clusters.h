
//
//  clusters.h
//  orbital_space_v4
//
//  Created by kuoyusheng on 2015/9/21.
//  Copyright (c) 2015年 kuoyusheng. All rights reserved.
//

#ifndef orbital_space_v4_clusters_h
#define orbital_space_v4_clusters_h
#include "util.h"
class Cluster{
public:
    static Cluster* create(int num, int size);
    
    Cluster();
    Cluster(int no_of_ver);
    Cluster(int no_of_ver, float ** vertex);
    //~Cluster();
    
    void assign_val(Cluster a);
    void map_data(float *clu, int no);
    void map_data_atom(int*clu, int no);
    bool comare();
    float* read_vertex();
    int  read_vertex_no();
    bool compare_transl_dist(Cluster b);
    Cluster* delete_dup(Cluster* a);
    bool compare_sym_dist(Cluster a, float ***linear_trans, float ***transl, int no_of_sym);
    void permu_vertex(int *pi);
    void index_elements(Cluster prim);
    void print_clu();
    void print_clu_atom();
    void clu_sym_cal(Cluster a,int sym_no, float **linear_trans,float *transl);
    void clu_imp_enum(Cluster , int, int*);
    float* find_red_fct(int& col);
    void imp_gamma_cal();
    void gamma_cal(float ***linear_trans);
    int no_of_sym;
    int isotromy=0;
    int no_of_os=0;
    int num_in_os=0;
    float *clu_vertex=NULL;
    int *clu_vertex_atom=NULL;
    int *permu=NULL;
    int no_of_ver=0;
    bool sym_distinctive=false;
    bool imp=false;
    bool rep=false;
    int *imp_idx=NULL;
    float clu_dist=0;
    float* gamma=NULL;
    float* imp_gamma=NULL;
    float* bmat=NULL;
    float* red_fct=NULL;
    int imp_ga_mulip=0;
    
private:
    static int dim;
    
};
void Cluster::index_elements(Cluster prim){
    clu_vertex_atom=new int[no_of_ver*(dim+1)];
    int no_ver_prim= prim.no_of_ver;
    for (int i=0; i<no_of_ver; i++) {
        for (int j=0; j<no_ver_prim; j++) {
            int* tmp_dist= new int[dim];
            for (int k=0,l=0; k<dim; k++) {
                double test=clu_vertex[i*dim+k]-prim.clu_vertex[j*dim+k];
                tmp_dist[k]=test;
                int test2=test;
                if (test==test2) {
                    l++;
                }
                if (l==dim) {
                    for(int m=0;m<dim;m++){
                        clu_vertex_atom[i*(dim+1)+m]=tmp_dist[m];
                    }
                    clu_vertex_atom[i*(dim+1)+dim]=j;
                    j=no_ver_prim;
                }
            }
            delete []tmp_dist;
        }
    }    
}

void Cluster::imp_gamma_cal(){
    if (!imp) {
        return;
    }
    int tensor_dim= power(3, no_of_ver);
    int** arr =new int*[tensor_dim];
    for (int i=0; i<tensor_dim; i++) {
        arr[i]=new int[no_of_ver];
        for(int j=0;j<no_of_ver;j++){
            arr[i][j]=0;
        }
    }
//    cout<<"array"
//    for(int i=0;i<tensor_dim;i++){
//        for (int j=0; j<no_of_ver; j++) {
//            cout<<arr[i][j]<<" ";
//        }
//    }
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
//    for(int i=0;i<k;i++){
//        cout<<permu[i][0]<<permu[i][1]<<endl;
//    }
    //    if(npt==3)
    //        cout<<k<<endl;
    imp_ga_mulip=k;
    imp_gamma=new float[(tensor_dim*k)*tensor_dim];
    for (int i=0; i<tensor_dim*k*tensor_dim; i++) {
        imp_gamma[i]=0;
    }
    
    
    tensor_ind(arr, tensor_dim, no_of_ver);
//        cout<<"array";
//        for(int i=0;i<tensor_dim;i++){
//            for (int j=0; j<no_of_ver; j++) {
//                cout<<arr[i][j]<<" ";
//            }cout<<endl;
//        }cout<<endl;
    
    for (int u=0; u<k; u++) {
        //cout<<"after permu"<<endl;
        for (int i=0; i<tensor_dim;i++ ) {
            swap(arr[i][permu[u][0]],arr[i][permu[u][1]]);
        }
        float **temp=new float*[tensor_dim];
        for (int i=0; i<tensor_dim; i++) {
            temp[i]=new float[tensor_dim];
            for(int j=0;j<tensor_dim;j++){
                temp[i][j]=0;
            }
        }
        
        for (int i=0; i<tensor_dim;i++ ) {
            int tmp=0;
            for (int j=0; j<npt; j++) {
                tmp+=(arr[i][j]*power(3, npt-j-1));
            }
            temp[i][tmp]=1;
        }
        for (int i=0; i<tensor_dim; i++) {
            for (int j=0; j<tensor_dim; j++) {
                imp_gamma[(u*tensor_dim+i)*tensor_dim+j]=temp[i][j];
            }
        }
        for(int i=0;i<tensor_dim;i++){
            delete [] temp[i];
        }
        delete []temp;
    }
    for(int i=0;i<tensor_dim;i++){
        delete[] permu[i];
    }
    delete []permu;
    for(int i=0;i<tensor_dim;i++){
        delete[] arr[i];
    }
    delete []arr;
}

void Cluster::gamma_cal(float ***linear_trans){
    int tensor_dim=power(dim, no_of_ver);
    int **tensor =new int*[tensor_dim];
    for (int i=0; i<tensor_dim; i++) {
        tensor[i]=new int[no_of_ver];
        for(int j=0;j<no_of_ver;j++){
            tensor[i][j]=0;
        }
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
    for(int i=0;i<tensor_dim;i++){
        delete[] tensor[i];
    }
    delete[] tensor;
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
void Cluster::print_clu_atom(){
    for (int i=0; i<no_of_ver; i++) {
        for (int j=0; j<(dim+1); j++) {
            cout<<clu_vertex_atom[i*(dim+1)+j]<<"\t";
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
void Cluster::map_data_atom(int*clu, int no){
    no_of_ver=no;
    clu_vertex_atom=new int[no_of_ver*(dim+1)];
    for (int i=0; i<no; i++) {
        for (int j=0; j<(dim+1);j++) {
            clu_vertex_atom[i*(dim+1)+j]=clu[i*(dim+1)+j];
        }
    }
}
//Cluster::~Cluster(){
//    delete []
//}
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
    //permu=new int[no_of_ver];
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
                distance[i][j]=clu_vertex[dim*i+j]-b.clu_vertex[dim*pi[idx][i]+j];
                int int_test=distance[i][j];
                if (distance[i][j]!=int_test)
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
    for(int i=0;i<no_of_ver;i++){
        delete[] distance[i];
    }
    delete []distance;
    return t=0;
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
                    //                    for (int k=0; k<no_of_ver; k++) {
                    //                        cout<<clu[idx*no_of_sym+j].permu[k]<<"\t";
                    //                    }cout<<endl;
                }
            }
        }
    }
}
void find_isotromy_pi(Cluster* clu_rep,Cluster* clu_sym,int no_of_ver, int no_of_rep, int no_of_sym){
    for (int idx=0,k=0; idx<no_of_rep; idx++,k=0) {
        for (int i=0; i<no_of_sym; i++) {
            if (clu_sym[idx*no_of_sym+i].compare_transl_dist(clu_rep[idx])) {
//                                for (int k=0; k<no_of_ver; k++) {
//                                    cout<<clu_sym[idx*no_of_sym+i].permu[k]<<'\t';
//                                }cout<<endl;
                //                k++;
            }
        }//cout<<k<<endl;
    }
}

Cluster* imp_clu_enu(int no_of_imp,int no_of_ver,int no_of_rep, Cluster* clu_rep,int & tot
                     ){
    int **com=combination_fixed_total(no_of_imp,no_of_ver);
    tot=factorial(no_of_imp+no_of_ver-1)/(factorial(no_of_imp)*factorial(no_of_ver-1));
    
    Cluster* imp_pair=Cluster::create(no_of_rep*tot, no_of_ver+no_of_imp);
    
    for (int idx=0; idx<no_of_rep; idx++) {
        for (int i=0; i<tot;i++) {
            imp_pair[idx*tot+i].clu_imp_enum(clu_rep[idx], no_of_imp, com[i]);
        }
    }
    for(int i=0;i<tot;i++){
        delete[] com[i];
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

void cluster_rep(int** arr_rep_no,Cluster **clusters_rep,int no_of_ver,float***linear_tran, float***transl, float cut_off[])
{
    for(int ver=1;ver<=no_of_ver;ver++)
    {
        int prop_no=arr_rep_no[ver][0];
        int imp_prop_no=prop_no;
        int tot=0;
        for (int i=1; i<ver; i++) {
            int no_of_imp=i;
            int no_ver=ver-no_of_imp;
            int multip=0;
            int no_of_dist=0;
            for(int j=0;j<arr_rep_no[no_ver][0];j++){
                if (clusters_rep[no_ver][j].clu_dist<=cut_off[no_of_ver]) {
                    tot++;
                }
            }
            //tot+=arr_rep_no[no_ver][0];
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
int* find_isotromy(Cluster* clu_sym,int st_p, int &no_os, int num_of_sym, int id_idx){
    int id=clu_sym[st_p+id_idx].isotromy;
    no_os=0;
    for (int i=0; i<num_of_sym; i++) {
        if (clu_sym[st_p+i].isotromy==id) {
            no_os++;
//            cout<<"this is "<<i<<'\t';
//            for(int k=0;k<clu_sym[st_p+i].no_of_ver;k++)
//                cout<<clu_sym[st_p+i].permu[k]<<'\t';
//            cout<<endl;
//            for(int k=0;k<power(3,clu_sym[st_p+i].no_of_ver);k++){
//                for(int j=0;j<power(3,clu_sym[st_p+i].no_of_ver);j++)
//                    cout<<clu_sym[st_p+i].gamma[(k*power(3,clu_sym[st_p+i].no_of_ver))+j]<<'\t';
//                cout<<endl;
//            }cout<<endl;
        }
    }
    int *arr =new int[no_os];
    for (int i=0,u=0;i<num_of_sym ; i++) {
        if (clu_sym[st_p+i].isotromy==id) {
            arr[u]=i;
            u++;
        }
    }
    
    clu_sym[arr[0]+st_p].rep=true;
    clu_sym[arr[0]+st_p].num_in_os=no_os;
    return arr;
}

void isotromy_gamma(Cluster* clu_sym, int no_of_sym, int tot_no, int id_idx){
    for (int clu=0; clu<tot_no; clu+=no_of_sym) {
        int no_os=0;
        int* arr=find_isotromy(clu_sym,clu, no_os, no_of_sym, id_idx);
        int num=clu_sym[clu+arr[0]].no_of_ver;
        int tensor_dim=power(3, num);
        int tot=(no_os*tensor_dim)*tensor_dim;
        float* temp=new float[tot];
        for (int i=0,u=0; i<no_os; i++) {
            for(int j=0;j<tensor_dim;j++){
                for(int k=0;k<tensor_dim;k++,u++){
                    temp[u]=clu_sym[clu+arr[i]].gamma[j*tensor_dim+k];
                    if (j==k) {
                        temp[u]--;
                    }
                }
            }
        }
//        for(int u=0;u<no_os*tensor_dim;u++){
//            for(int i=0;i<tensor_dim;i++){
//                cout<<temp[u*tensor_dim+i]<<'\t';
//            }cout<<endl;
//        }
        //for(int i=0;i<no_os;i++)
        if (clu_sym[clu+arr[0]].imp) {
            int k=clu_sym[clu+arr[0]].imp_ga_mulip;
            int row_imp=tensor_dim*(k+no_os);
            int tot_imp=row_imp*tensor_dim;
            clu_sym[clu+arr[0]].bmat=new float[tot_imp];
            for (int i=0; i<tot; i++) {
                clu_sym[clu+arr[0]].bmat[i]=temp[i];
            }
            int* temp_imp=new int[k*tensor_dim*tensor_dim];
            for (int u=0; u<k; u++) {
                for (int i=0; i<tensor_dim; i++) {
                    for (int j=0;j<tensor_dim ; j++) {
                        temp_imp[(u*tensor_dim+i)*tensor_dim+j]=clu_sym[clu+arr[0]].imp_gamma[(u*tensor_dim+i)*tensor_dim+j];
                        if (i==j) {
                            temp_imp[(u*tensor_dim+i)*tensor_dim+j]--;
                        }
                    }
                }
            }
            for (int i=tot; i<tot_imp; i++) {
                clu_sym[clu+arr[0]].bmat[i]=temp_imp[i-tot];
            }
            delete[] temp_imp;
        }
        else{
            //cout<<no_os<<endl;;
            //cout<<"proper"<<endl;
            clu_sym[clu+arr[0]].bmat=new float[tot];
            for (int i=0; i<tot; i++) {
                clu_sym[clu+arr[0]].bmat[i]=temp[i];
            }
//            if (clu_sym[clu+arr[0]].no_of_ver==3) {
//                for(int i=0;i<(tot/tensor_dim);i++){
//                    for(int j=0;j<tensor_dim;j++){
//                        cout<<clu_sym[clu+arr[0]].bmat[i*tensor_dim+j]<<'\t';
//                    }cout<<endl;
//                }cout<<endl;
//            }
        }
        delete []temp;
    }
}


void print_max(float* v, int row, int col){
    for (int i=0; i<row; i++) {
        for (int j=0; j<col; j++) {
            cout<<v[i*col+j]<<"\t";
        }cout<<endl;
    }cout<<endl;
}
float* Cluster::find_red_fct(int& col){
    int ten_dim=power(3, no_of_ver);
    if (imp) {
        cout<<"imp"<<endl;
        cout<<imp_ga_mulip<<endl;
        int tot=(ten_dim*(num_in_os+imp_ga_mulip))*ten_dim;
        double* b_out_tmp=new double[tot];
        double*v= new double[tot];
        for (int i=0; i<tot; i++) {
            v[i]=bmat[i];
        }
//        for (int i=0; i<ten_dim*(1+imp_ga_mulip); i++) {
//            for (int j=0; j<ten_dim; j++) {
//                cout<<bmat[i]<<'\t';
//            }cout<<endl;
//        }
        int row=(num_in_os+imp_ga_mulip)*ten_dim;
        //cout<<row<<endl;
        col=nullspace(row, ten_dim, v, b_out_tmp);
        cout<<col<<endl;
        red_fct=new float[ten_dim*col];
        for (int i=0; i<row; i++) {
            for (int j=0; j<col; j++) {
                red_fct[i*col+j]=b_out_tmp[i*col+j];
            }
        }
        print_max(red_fct, ten_dim, col);
        //cout<<endl;

        delete []b_out_tmp;
        delete []v;
        return red_fct;

    }
    else{
        cout<<"prop"<<endl;
        int tot=(ten_dim*num_in_os)*ten_dim;
        if(no_of_ver==3)
            cout<<num_in_os<<endl;
        double* b_out_tmp=new double[tot];
        double*v= new double[tot];
        for (int i=0; i<tot; i++) {
            v[i]=bmat[i];
        }
//        if (no_of_ver==3) {
//            for(int i=0;i<ten_dim*num_in_os;i++){
//                for(int j=0;j<ten_dim;j++)
//                    cout<<v[ten_dim*i+j]<<'\t';
//                cout<<endl;
//            }cout<<endl<<endl;;
//        }
        col=nullspace(ten_dim*num_in_os, ten_dim, v, b_out_tmp);
        red_fct=new float[ten_dim*col];
        int row=ten_dim;
        cout<<col<<endl;
        
        for (int i=0; i<row; i++) {
            for (int j=0; j<col; j++) {
                red_fct[i*col+j]=b_out_tmp[i*col+j];
            }
        }
        print_max(red_fct, row, col);
        delete []b_out_tmp;
        delete []v;
        return red_fct;

    }
}
float* gamma_cal(float **linear_trans, int* permu, int no_of_ver){
    int tensor_dim=power(3, no_of_ver);
    int **tensor =new int*[tensor_dim];
    for (int i=0; i<tensor_dim; i++) {
        tensor[i]=new int[no_of_ver];
    }
    float* gamma;
    gamma=new float[tensor_dim*tensor_dim];
    tensor_ind(tensor, tensor_dim, no_of_ver);
    for (int i=0; i<tensor_dim; i++) {
        for (int j=0; j<tensor_dim; j++) {
            float tmp=1;
            for (int k=0; k<no_of_ver; k++) {
                tmp*=linear_trans[tensor[i][k]][tensor[j][permu[k]]];
                if (tmp==(-0)) {
                    tmp=0;
                }
            }
            gamma[i*tensor_dim+j]=tmp;
        }
    }
    for(int i=0;i<tensor_dim;i++){
        delete [] tensor[i];
    }
    delete []tensor;
    return gamma;
}
float* transl_invariant_conversion(Cluster rep, Cluster& aa, float***linear_trans, int& igp){
    int no_of_ver=rep.no_of_ver;
    int tot_permu=factorial(no_of_ver);
    int **pi=permu_fun(no_of_ver);
    int dim=3;
    
    float** distance=new float*[no_of_ver];
    for (int i=0;i<no_of_ver ; i++) {
        distance[i]=new float[3];
    }
    bool t=true;
    for(int idx=0;idx<tot_permu;idx++,t=1){
        for (int i=0; i<no_of_ver; i++) {
            if(rep.clu_vertex_atom[pi[idx][i]*dim+dim]!=aa.clu_vertex_atom[i*dim+dim]){
                t=0;
                break;
            }
            else for(int j=0;j<dim;j++)
                distance[i][j]
                =rep.clu_vertex_atom[pi[idx][i]*(dim+1)+j]-aa.clu_vertex_atom[i*(dim+1)+j];
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
            igp=rep.no_of_sym;
            aa.no_of_sym=igp;
            aa.no_of_ver=no_of_ver;
            for(int i=0;i<no_of_ver;i++){
                aa.permu[i]=pi[idx][i];
            }
            //aa.permu=pi[idx];
            aa.gamma_cal(linear_trans);
//            //cout<<igp<<endl;
//            if(igp==10)
//            for (int i=0; i<9; i++) {
//                for(int j=0;j<9;j++)
//                    cout<<aa.gamma[i*9+j]<<'\t';
//                cout<<endl;
//            }
            for(int i=0;i<no_of_ver;i++){
                delete distance[i];
            }
            delete [] distance;
            for(int i=0;i<tot_permu;i++){
                delete [] pi[i];
            }
            delete [] pi;
            return aa.gamma;
        }
    }
    for(int i=0;i<no_of_ver;i++){
        delete distance[i];
    }
    delete [] distance;
    for(int i=0;i<tot_permu;i++){
        delete [] pi[i];
    }
    delete [] pi;
    return NULL;
}

Cluster* find_all_clu(Cluster prim,int** all_atom,int num_atom, int no_of_ver, int atom_type,int&tot){
    tot=power(num_atom,no_of_ver-1);
    Cluster* all_clu=Cluster::create(tot,no_of_ver);
    int **tensor=new int*[tot];
    for (int i=0; i<tot; i++) {
        tensor[i]=new int[no_of_ver-1];
    }
    find_np(tensor, tot, num_atom, no_of_ver-1);
    for(int i=0;i<tot;i++){
        for (int j=0; j<no_of_ver-1; j++) {
            cout<<tensor[i][j]<<'\t';
        }cout<<endl;
    }
    for (int i=0; i<tot; i++) {
        int* clu=new int[no_of_ver*4];
        for (int j=0; j<no_of_ver; j++) {
            if(j==0){
                for (int k=0; k<4; k++) {
                    clu[k]=prim.clu_vertex_atom[atom_type*4+k];
                }
            }
            else {
                for (int k=0; k<4; k++) {
                    clu[j*4+k]=all_atom[tensor[i][j-1]][k];
                }
            }
        }
        all_clu[i].map_data_atom(clu,no_of_ver);
        delete[] clu;
    }
    return all_clu;
}























#endif
