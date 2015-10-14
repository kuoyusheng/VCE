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
#include </Users/kuoyusheng/Desktop/Eigen/Eigen>// works ok
#include "null_space.h"
#include "clusters.h"
#include "util.h"
using namespace std;
using Eigen::MatrixXd;

bool examine_cut_off(int* arr,float* coord, float cut_off_radius){
    int j=arr[3];
    Eigen::MatrixXd tran(3,3);
    Eigen::MatrixXd dir(3,1);
    Eigen::MatrixXd cart(3,1);
    for (int i=0; i<3; i++) {
        dir(i,0)=arr[i]+coord[9+j*3+i];
        for (int k=0; k<3; k++) {
            tran(i,k)=coord[i*3+k];
        }
    }
    cart=tran*dir;
    float length=cart.squaredNorm();
    
    if (length<=cut_off_radius) {
        return true;
    }
    else return false;
}
int** find_all_atom(double *l_para,float cut_off_r, int no_of_atom, int& tot, float* lat_in){
    int cut_off[3];
    for (int i=0; i<3; i++) {
        cut_off[i]=cut_off_r/l_para[i];
    }
    int pos_tot=1;
    for (int i=0; i<3; i++) {
        if (cut_off[i]) {
            pos_tot*=(cut_off[i]*2+1);
        }
    }
    int** arr1=new int*[pos_tot*no_of_atom];
    for (int i=0; i<pos_tot*no_of_atom; i++) {
        arr1[i]=new int[4];
    }
    for (int i=(-cut_off[0]),l=0; i<=cut_off[0]; i++) {
        for (int j=(-cut_off[1]); j<=cut_off[1]; j++) {
            for (int k=(-cut_off[2]); k<=cut_off[2]; k++) {
                for (int at=0; at<no_of_atom; at++) {
                    arr1[l][0]=i;
                    arr1[l][1]=j;
                    arr1[l][2]=k;
                    arr1[l][3]=at;
                    l++;
                }
            }
        }
    }
    
    tot=0;
    int **temp=new int*[pos_tot*no_of_atom];
    for (int i=0; i<pos_tot*no_of_atom;i++) {
        temp[i]=new int[4];
        if (examine_cut_off(arr1[i], lat_in, cut_off_r)) {
            for (int j=0; j<4; j++) {
                temp[tot][j]=arr1[i][j];
            }
            tot++;
        }
    }
    int **arr=new int*[tot];
    for (int i=0; i<tot; i++) {
        arr[i]=new int[4];
        for (int j=0; j<4; j++) {
            arr[i][j]=temp[i][j];
        }
    }
    for (int i=0; i<tot; i++) {
        delete [] temp[i];
    }
    delete[] temp;
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
//    for (int i=0; i<3; i++) {
//        for (int j=0; j<3; j++) {
//            cout<<linear_tran_c[6][i][j];
//        }
//    }cout<<endl;
    //return 1;
    
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
    int max_ver=(3+1);//no of ver plus zero order;
    int single_ind=0;
    int pair_ind=0;
    int tri_ind=0;
    float cut_off[4]={0,0,1.45,0.8};
    float** distance=new float*[max_ver];
    for (int i=0; i<max_ver; i++) {
        distance[i]=new float[100];
    }
    //float distance[100][;
    int k=0;
    int ver=1;
    //    float dist=1.5;
    //    float dist_tri=1.0;
    //acquire rep. cluster
    while(getline(clu,cluster)){
        float dist=0;
        for(int i=0;i<2;i++){
            getline(clu,cluster);
            istringstream s(cluster);
            if (i==0) {
                s>>dist;
            }
            else{
                int tmp=ver;
                s>>ver;
                if (tmp==ver)
                    k++;
                else
                    k=0;
                distance[ver][k]=dist;
            }
        }
        if(dist>cut_off[ver]){
            //getline(clu,cluster);
            for (int i=0; i<ver; i++) {
                getline(clu,cluster);
            }
            getline(clu,cluster);
            continue;
        }
        else if (cluster=="0"){
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
        clusters_rep[2][i].clu_dist=distance[2][i];
    }
    for (int i=0; i<arr_rep_no[1][0]; i++) {
        clusters_rep[1][i].map_data(single[i], 1);
    }
    for(int i=0;i<arr_rep_no[3][0];i++){
        clusters_rep[3][i].map_data(triplet[i], 3);
        clusters_rep[3][i].clu_dist=distance[3][i];
    }
//    
//    for (int i=0; i<arr_rep_no[3][0]; i++) {
//            clusters_rep[3][i].print_clu();
//    }

//
    int no_of_ver=3;

    cout<<"Enumerate imp clusters"<<endl;
    cluster_rep(arr_rep_no, clusters_rep, 3, linear_tran, transl,cut_off);
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
    cout<<"#of improper"<<endl;
    for (int i=1; i<=3; i++) {
        for (int j=0; j<4; j++) {
            clu_num[i]+=arr_rep_no[i][j];
            if (i==3) {
                cout<<arr_rep_no[i][j]<<'\t';
            }
        }
    }cout<<endl;
        for (int i=0; i<clu_num[3]; i++) {
            clusters_rep[3][i].print_clu();
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
        clu_sym[i]=Cluster::create(2000,i);
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
    cout<<"find reduced Cmat"<<endl;
    cout<<clu_num[3]<<endl;
    isotromy_gamma(clu_sym[3], no_of_sym, clu_num[3]*no_of_sym, 10);
//    for (int i=0; i<clu_num[3]*48; i++) {
//        if (clu_sym[3][i].rep) {
//            cout<<clu_sym[3][i].num_in_os;
//            cout<<endl;
//        }
//    }
    //cout<<"this is second";
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
    int tensor_dim=power(3, no_of_ver);
    float** Cmat=new float*[tensor_dim*clu_num[no_of_ver]];
    for (int i=0; i<tensor_dim*clu_num[no_of_ver]; i++) {
        Cmat[i]=new float[tensor_dim*clu_num[no_of_ver]];
        for (int j=0; j<tensor_dim*clu_num[no_of_ver]; j++) {
            Cmat[i][j]=0;
        }
    }
    int no=2;
    Cluster prim(no);
    float prim_ver[]={0,0,0,0.25,0.25,0.25};
    for (int i=0; i<2*3; i++) {
        prim.clu_vertex[i]=prim_ver[i];
    }
    prim.index_elements(prim);
    
    for (int i=0; i<no_of_ver; i++) {
        for (int j=0;j<clu_num[i]*no_of_sym; j++) {
            clu_sym[i][j].index_elements(prim);
            //cout<<"dist: "<<clu_sym[i][j].clu_dist<<endl;
            //clu_sym[i][j].print_clu_atom();
        }
    }
//    for (int i=10; i<clu_num[2]*no_of_sym; i+=48) {
//        clu_sym[2][i].print_clu_atom();
//    }
    int col_tot=0;
    int num=0;
    for (int i=0; i<clu_num[no_of_ver]*no_of_sym; i++) {
        if (clu_sym[no_of_ver][i].rep) {
            //clu_sym[2][i].print_clu_atom();
            //cout<<i%48<<endl;
            int s_row=num*tensor_dim;
            int s_col=col_tot;
            cout<<"s_row:"<<num*tensor_dim<<endl;
            cout<<"s_col:"<<col_tot<<endl;

            int col=0;
            float*fct = clu_sym[no_of_ver][i].find_red_fct(col);
                for (int j=0; j<tensor_dim; j++) {
                    for (int k=0; k<col; k++) {
                        if(abs(fct[j*col+k])<1e-9)
                            fct[j*col+k]=0;
                        Cmat[s_row+j][s_col+k]=fct[j*col+k];
                    }
                }
            num++,col_tot+=col;
        }
    }
    cout<<endl;
    //cout<<num<<endl;
    cout<<"tot col of Cmat:"<<col_tot<<endl;
    cout<<"Cmat: "<<endl;
    //return 0;
   // Eigen::MatrixXd Cmat_iso(9*clu_num[2],col_tot);
    float* Cmat_new=new float[9*clu_num[2]*col_tot];
    for (int i=0; i<9*clu_num[2]; i++) {
        for (int j=0; j<col_tot; j++) {
            Cmat_new[i*col_tot+j]=Cmat[i][j];
        }
    }
    
    //return 1;

//    int no=2;
//    Cluster prim(no);
//    float prim_ver[]={1,1,1,0.25,0.25,0.25};
//    for (int i=0; i<2*3; i++) {
//        prim.clu_vertex[i]=prim_ver[i];
//    }
//    prim.index_elements(prim);
//   
//    for (int i=0; i<no_of_ver; i++) {
//        for (int j=0;j<clu_num[i]*no_of_sym; j++) {
//            clu_sym[i][j].index_elements(prim);
//            //clu_sym[i][j].print_clu_atom();
//        }
//    }
    num=0;
    int no_of_ver_1=2;
    int no_of_atom=2;
    float lat_in[]={0,0.5,0.5,0.5,0,0.5,0.5,0.5,0,0,0,0,0.25,0.25,0.25};
    double l_para[3]={sqrt(0.25),sqrt(0.25),sqrt(0.25)};
    //find out all the atom by iterating from 0000 to cut_off atom
    int **all_atom=find_all_atom(l_para,cut_off[2],no_of_atom,num,lat_in);
    for (int i=0; i<num; i++) {
            for (int k=0; k<4; k++) {
                cout<<all_atom[i][k]<<'\t';
            }cout<<endl;
    }
    //return 1;
   // int no=power(num, no_of_ver_1-1);
    Cluster* all_np=Cluster::create(power(num, no_of_ver_1-1), no_of_ver_1);
    //prim_ind choose what atom in prim cell
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
        for (int j=0; j<ten_dim*no_rep; j++) {
            Bmat[i][j]=0;
        }
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
//            Bmat_tran(i,j)=Bmat[i][j];
//        }
//    }
    float* bmat=new float[ten_dim*ten_dim*clu_num[2]];
    double* b_out=new double[ten_dim*ten_dim*col_tot];
    for (int i=0; i<ten_dim; i++) {
        for (int j=0; j<ten_dim*clu_num[2]; j++) {
            bmat[i*ten_dim*clu_num[2]+j]=Bmat[i][j];
        }
    }
//    for (int i=0; i<ten_dim; i++) {
//        for (int j=0; j<ten_dim*clu_num[2]; j++) {
//            cout<<bmat[i*ten_dim*clu_num[2]+j]<<" ";
//        }cout<<endl;
//    }cout<<endl;
    cout<<"bmat[1]"<<bmat[1]<<endl;
    cout<<"Cmat[1]"<<Cmat_new[1]<<endl;
    float* _Bmat=Matrix_product(ten_dim, ten_dim*clu_num[2], ten_dim*clu_num[2],col_tot, bmat, Cmat_new);
    double* _bmat=new double[ten_dim*col_tot];
    for (int i=0;i<ten_dim; i++) {
        for (int j=0; j<col_tot; j++) {
            _bmat[i*col_tot+j]=_Bmat[i*col_tot+j];
        }
    }
    
    //cout<<_bmat[1];
    {
    ofstream fout("/Users/kuoyusheng/Desktop/bmat2.txt");
    if (!fout) {
        cerr << "Could not open file." << endl;
        return 1;
    }
    for (int i=0; i<ten_dim; i++) {
        for (int j=0; j<ten_dim*clu_num[2]; j++) {
//            if (_bmat!=0) {
//                //fout<<"n";
//            }
            fout<<Bmat[i][j]<<" ";
        }fout<<endl;
    }
    fout.close();
    }
    {
        ofstream fout("/Users/kuoyusheng/Desktop/Cmat3.txt");
        if (!fout) {
            cerr << "Could not open file." << endl;
            return 1;
        }
        for (int i=0; i<ten_dim*clu_num[2]; i++) {
            for (int j=0; j<col_tot; j++) {
                //            if (_bmat!=0) {
                //                //fout<<"n";
                //            }
                fout<<Cmat[i][j]<<" ";
            }fout<<endl;
        }
        fout.close();
    }
    ofstream fout("/Users/kuoyusheng/Desktop/result4.txt");
    if (!fout) {
        cerr << "Could not open file." << endl;
        return 1;
    }
    for (int i=0; i<ten_dim; i++) {
        for (int j=0; j<col_tot; j++) {
            if (_bmat!=0) {
                //fout<<"n";
            }
            fout<<_Bmat[i*col_tot+j]<<" ";
        }fout<<endl;
    }
    fout.close();
    
    int col=nullspace(ten_dim, col_tot, _bmat, b_out);
    cout<<col<<'\t'<<ten_dim*clu_num[2]<<endl;;
    
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
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    sym.close();
    clu.close();
    
    
    
    
    
    
    
    return 0;
}
