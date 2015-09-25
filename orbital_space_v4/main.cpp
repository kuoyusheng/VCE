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
#include "clusters.h"
#include "util.h"
using namespace std;
using Eigen::MatrixXd;


int** find_all_atoms(int cut_off, int no_of_atom, int& no){
    int a=power(cut_off, 3);
    int **arr1=new int*[power(cut_off, 3)];
    for (int i=0; i<power(cut_off, 3); i++) {
        arr1[i]=new int[3];
    }
    int b=power(no_of_atom, 3);
    int **arr2=new int*[b];
    for (int i=0; i<b; i++) {
        arr2[i]=new int[3];
    }
    int c=power(cut_off, 3)*power(no_of_atom, 3)*no_of_atom;
    int** arr3=new int*[c];
    for (int i=0; i<c; i++) {
        arr3[i]=new int[4];
    }
    find_np(arr1,a,cut_off);
    find_np(arr2,b, 2);
    int l=0;
    for (int i=0; i<a; i++) {
        for (int j=0; j<b; j++) {
            for (int k=0,u=0; k<3; k++) {
                if ((arr1[i][k]||arr2[j][k])) {
                    u++;
                }
                if (u==3) {
                    for (int atom=0; atom<no_of_atom; atom++) {
                        for (int m=0; m<3; m++) {
                            arr3[l][3]=atom;
                            if (arr2[j][m]==1) {
                                arr3[l][m]=arr1[i][m];
                            }
                            else arr3[l][m]=-1*arr1[i][m];
                        }
                        l++;
                    }
                }
            }
        }
    }
    no=l;
    int **arr=new int*[l];
    for (int i=0; i<l; i++) {
        arr[i]=new int[4];
    }
    for (int i=0; i<l; i++) {
        for (int j=0; j<4; j++) {
            arr[i][j]=arr3[i][j];
        }
    }
    delete[] arr1;
    delete[] arr2;
    delete[] arr3;
    return arr;
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
    float ***linear_tran_c=new float**[no_of_sym];
    for (int i=0; i<no_of_sym; i++) {
        linear_tran_c[i]=new float*[3];
        for (int j=0; j<3; j++) {
            linear_tran_c[i][j]=new float[3];
        }
        
    }
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
    Eigen::MatrixXd trans(3,3);
    trans<<0,1,1,1,0,1,1,1,0;
    symop_DtoC(linear_tran, linear_tran_c, no_of_sym, trans);
//    for (int i=0; i<no_of_sym; i++) {
//        for (int j=0; j<3; j++) {
//            for (int k=0; k<3; k++) {
//                cout<<linear_tran_c[i][j][k]<<"\t";
//            }cout<<endl;
//        }cout<<endl<<endl;
//    }
    
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
//    for (int i=0; i<no_of_ver; i++) {
//        cout<<i<<"th order clusters"<<endl;
//        for (int j=0; j<clu_num[i]*no_of_sym; j++) {
//            clu_sym[i][j].print_clu();
//        }
//    }
    
//    cout<<"Calculate"<<endl;
//    
//        for (int i=0; i<no_of_sym; i++) {
//            clu_sym[i].clu_sym_cal(clusters_rep[2][0], i, linear_tran[i], transl[i][0]);
//        }
//
    cout<<"Find Orbital Space"<<endl;
    for (int ver=1; ver<=no_of_ver; ver++) {
        find_orbital_space(clu_sym[ver], ver, clu_num[ver], no_of_sym);
        find_isotromy_pi(clusters_rep[ver], clu_sym[ver], ver, clu_num[ver], no_of_sym);
    }
    
    
//    cout<<"Find isotromy permu"<<endl;
//    for (int ver=1; ver<=no_of_ver; ver++) {
//
//    }
    
    
//
//
//    int array[12]={};
//    for (int i=0,u=0; i<no_of_sym; i++) {
//        if (clu_sym[2][i].isotromy==2) {
//            array[u]=i;
//            cout<<i<<"\t";
//            u++;
//        }
//    }cout<<endl;
//    int *arr=find_isotromy(clu_sym[3], no_of_sym, 10);
//    for (int i=0; i<4; i++) {
//        cout<<arr[i]<<"\t";
//    }
    
//    cout<<clu_num[3]<<endl;
    cout<<"Find gamma conversion"<<endl;
    for (int ver=2; ver<=no_of_ver; ver++) {
        for (int i=0; i<clu_num[ver]*48; i++) {
            clu_sym[ver][i].gamma_cal(linear_tran_c);
            //clu_sym[ver][i].gamma_cal(linear_tran);
            clu_sym[ver][i].imp_gamma_cal();
        }
    }

    //cout<<clu_sym[3][34*48].imp_ga_mulip<<endl;
//    for (int i=0; i<27*clu_sym[3][34*48].imp_ga_mulip; i++) {
//        for (int j=0; j<27; j++) {
//            cout<<clu_sym[3][34*48].imp_gamma[i*27+j]<<"\t";
//        }cout<<endl;
//    }
//    for (int i=24*48,k=0; i<clu_num[2]*48; i++,k++) {
////        if (clu_sym[2][i].no_of_os==0) {
////            cout<<i%48;
////        }
//        if (clu_sym[2][i].isotromy==0&&clu_sym[2][i].imp) {
////            for (int j=0; j<2; j++) {
////                cout<<clu_sym[2][i].permu[j]<<"\t";
////            }cout<<endl;
//            for (int j=0; j<9&&k==0; j++) {
//                for(int k=0;k<9;k++)
//                    cout<<clu_sym[2][i].imp_gamma[j*9+k]<<"\t";
//                cout<<endl;
//            }
//        }
//    }
    cout<<"find reduced C"<<endl;;
    isotromy_gamma(clu_sym[3], no_of_sym, clu_num[3]*no_of_sym, 10);
    isotromy_gamma(clu_sym[2], no_of_sym, clu_num[2]*no_of_sym, 10);

    
//        for (int i=0; i<27*clu_sym[3][34*48].imp_ga_mulip; i++) {
//            for (int j=0; j<27; j++) {
//                cout<<clu_sym[3][34*48].imp_gamma[i*27+j]<<"\t";
//            }cout<<endl;
//        }cout<<endl<<endl;
//    for (int i=0; i<9*25; i++) {
//                for (int j=0; j<9; j++) {
//                    cout<<clu_sym[2][24*48].gamma[i*9+j]<<"\t";
//                }cout<<endl;
//            }
    float** Cmat=new float*[9*clu_num[2]];
    for (int i=0; i<9*clu_num[2]; i++) {
        Cmat[i]=new float[9*clu_num[2]];
    }
    int col_tot=0;
    int num=0;
    for (int i=0; i<clu_num[2]*48; i++) {
        if (clu_sym[2][i].rep) {
            //cout<<i%48<<endl;
            int col=0;
            float*fct = clu_sym[2][i].find_red_fct(col);
                for (int j=0; j<9; j++) {
                    for (int k=0; k<9; k++) {
                        if(fct[j*9+k]<1e-9)
                            fct[j*9+k]=0;
                        Cmat[num*9+j][col_tot+k]=fct[j*9+k];
                    }
                }
            num++,col_tot+=col;
        }
    }cout<<num<<endl;
    cout<<col_tot<<endl;
    cout<<"Cmat: "<<endl;
    float** Cmat_new=new float*[9*clu_num[2]];
    for (int i=0; i<9*clu_num[2]; i++) {
        Cmat_new[i]=new float[col_tot];
        for (int j=0; j<col_tot; j++) {
            Cmat_new[i][j]=Cmat[i][j];
        }cout<<endl;
    }
    return 1;

    int no=2;
    Cluster prim(no);
    float prim_ver[]={1,1,1,0.25,0.25,0.25};
    for (int i=0; i<2*3; i++) {
        prim.clu_vertex[i]=prim_ver[i];
    }
    prim.index_elements(prim);
   
    for (int i=0; i<no_of_ver; i++) {
        for (int j=0;j<clu_num[i]*no_of_sym; j++) {
            clu_sym[i][j].index_elements(prim);
        }
    }
    num=0;
    int no_of_ver_1=2;
    int **all_atom=find_all_atoms(4, 2, num);
    Cluster* all_np=Cluster::create(power(num, no_of_ver_1-1), no_of_ver_1);
    for (int prim_ind=0; prim_ind<no; prim_ind++) {
        for (int i=0; i<num; i++) {
            int* clu=new int[no_of_ver_1*4];
            int j=0;
            for (j=0; j<4; j++) {
                clu[j]=prim.clu_vertex_atom[prim_ind*4+j];
            }
            for (int k=0; k<4; k++,j++) {
                clu[j]=all_atom[i][k];
            }
            all_np[i].map_data_atom(clu, no_of_ver_1);
            delete[] clu;
        }
    }
    int igp=0;
    int no_rep=clu_num[2]*no_of_sym;
    int ten_dim=power(3, no_of_ver_1);
    float** Bmat=new float*[ten_dim];
    for (int i=0; i<ten_dim; i++) {
        Bmat[i]=new float[ten_dim*no_rep];
    }
    for (int ip=0; ip<num; ip++) {
        for (int idx=0; idx<no_rep; idx++) {
            float* gamma=transl_invariant_conversion(clu_sym[no_of_ver_1][idx], all_np[ip], linear_tran_c,igp);
            if (gamma==NULL)
                continue;
            int s_col=idx/no_of_sym;
            for (int i=0;i<ten_dim; i++) {
                for (int j=ten_dim*s_col,l=0; l<ten_dim; j++,l++) {
                    Bmat[i][j]+=gamma[i*ten_dim+l];
                }
            }
        }
    }
//    for (int i=0; i<ten_dim; i++) {
//        for (int j=0; j<ten_dim*no_rep; j++) {
//            cout<<Bmat[i][j]<<'\t';
//        }cout<<endl;
//    }
    double* bmat=new double[ten_dim*ten_dim*clu_num[2]];
    double* b_out=new double[ten_dim*ten_dim*clu_num[2]];
    for (int i=0; i<ten_dim; i++) {
        for (int j=0; j<ten_dim*clu_num[2]; j++) {
            bmat[i*ten_dim*clu_num[2]+j]=Bmat[i][j];
        }
    }
    int col=nullspace(ten_dim, ten_dim*clu_num[2], bmat, b_out);
    cout<<col<<'\t'<<ten_dim*clu_num[2];
    
    
    
   
    
    
    
    
    
    
    
//    for (int i=48; i<2*48; i++) {
//        if (clu_sym[3][i].no_of_os==0) {
//            cout<<i-48<<"\t";
//        }
//    }
//    for (int i=48; i<2*48; i++) {
//        if (clu_sym[3][i].isotromy==6) {
//            cout<<i-48<<"\t";
//        }
//    }
//    for (int i=0; i<6*27; i++) {
//        for (int j=0; j<27; j++) {
//            cout<<clu_sym[3][54].gamma[i*27+j]<<"\t";
//        }cout<<endl;
//    }
    //isotromy_gamma(clu_sym[2], array, 12);
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
    
//
    
//    double* ver=new double[6*27*27];
//    int tot=6*27*27;
//    for (int i=0; i<tot; i++) {
//        ver[i]=clu_sym[3][54].gamma[i];
//    }
//    double* b_out=new double[6*27*27];
//    int col=nullspace(6*27, 27, ver, b_out);
//    
//    cout<<col<<endl;
//    for (int i=0; i<27*6; i++) {
//        for (int j=0;j<col; j++) {
//            cout<<b_out[i*col+j]<<"\t";
//        }cout<<endl;
//    }
    
    
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
