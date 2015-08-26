//
//  Header.h
//  orbital_space_v4
//
//  Created by kuoyusheng on 2015/8/19.
//  Copyright (c) 2015年 kuoyusheng. All rights reserved.
//

#ifndef orbital_space_v4_Header_h
#define orbital_space_v4_Header_h
void FCT_conversion(int** tensor,int***permu,float***linear_tran,float**gamma,int**o_s,int ind_of_rep, int ind_of_iso,
                    int dim,int no_of_ver){
    
    int n_s=o_s[ind_of_rep][ind_of_iso];
    for (int j=0; j<dim; j++) {
        for (int k=0; k<dim; k++) {
            float tmp=1;
            for (int u=0; u<no_of_ver; u++) {
                tmp*=linear_tran[n_s] [tensor[j][u]][tensor[k][permu[ind_of_rep][n_s][u]]];
                if (tmp==(-0))
                    tmp=0;
                
            }
            gamma[j][k]=tmp;
        }
    }
}
void tensor_ind(int **tensor,int dim, int no_of_ver){
    for (int i=0; i<dim; i++) {
        combination(i,0,tensor[i],no_of_ver);
    }
}

void combination(int dim,int i,int arr[],int digit)//
{
    if (dim / 3 != 0) {
        combination((dim / 3),i+1,arr,digit);
    }
    arr[digit-1-i]=(dim % 3);
}




cout<<"Enumerate imp clusters"<<endl;
int prop_no=arr_rep_no[no_of_ver][0];
int imp_prop_no=prop_no;
for (int i=1; i<no_of_ver; i++) {
    int no_of_imp=i;
    int no_ver=no_of_ver-no_of_imp;
    int multip=0;
    int no_of_dist=0;
    int tot=0;
    for (int j=0; j<4; j++) {
        tot+=arr_rep_no[no_ver][j];
    }
    Cluster* imp=imp_clu_enu(no_of_imp, no_ver, tot, clusters_rep[no_ver], multip);
    Cluster* imp_clu_rep=delete_dup(imp, tot*multip,no_of_dist,linear_tran, transl);
    imp_prop_no+=no_of_dist;
    for (int k=prop_no,j=0; k<imp_prop_no; k++,j++) {
        clusters_rep[no_of_ver][k]=imp_clu_rep[j];
    }
    prop_no=imp_prop_no;
}
for (int i=0; i<imp_prop_no; i++) {
    if(clusters_rep[no_of_ver][i].clu_dist<dist_tri)
        clusters_rep[no_of_ver][i].print_clu();
}


0.25	0.25	0.25
0	0	0
0.25	0.25	-0.75

0.25	0.25	0.25
1.25	0.25	-0.75
0.25	0.25	-0.75

0.25	0.25	0.25
0.25	-0.75	0.25
0.25	0.25	-0.75

0.25	0.25	0.25
0	0	0
1	0	-1

0.25	0.25	0.25
0.25	1.25	-0.75
1	0	-1

0.25	0.25	0.25
0.25	-0.75	0.25
1	0	-1

1	1	1
1.25	1.25	1.25
1.25	1.25	1.25

0.25	0.25	0.25
0.25	0.25	-0.75
0.25	0.25	-0.75

0.25	0.25	0.25
1	0	-1
1	0	-1

1	1	1
1	1	1
1	1	1




































#endif
