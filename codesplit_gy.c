#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "sys/time.h"

#define pie 3.1415926

void analysis (char *infile, char * outfile, char * logfile, char **argv)
{
	strcpy(infile,  argv[1]);
    strcpy(outfile, argv[2]);
    strcpy(logfile, argv[3]);
    printf("123\n");
}

void log_start (char *logfile, int lenght)
{
	char tmp[80];
	FILE  *flog;

	strcpy(tmp, "date ");
    strncat(tmp, ">> ",3);
    strncat(tmp, logfile, lenght);

    flog = fopen(logfile, "w");
    fprintf(flog, "------------start time------------\n");
    fclose(flog);

    system(tmp);  
}
void read_data(char *infile, int *nx, int *ny, int *nz, int *lt, int *nedge, int *ncx_shot1, int *ncy_shot1, int *ncz_shot, int *nxshot, int *nyshot, float *frequency, float *velmax, float *dt, float *unit,int *dxshot,int *dyshot)
{
	FILE *fin;

	fin = fopen(infile,"r");
    if(fin == NULL)
    {
        printf("file %s is  not exist\n", infile);
        exit(0);
    }
    fscanf(fin,"nx=%d\n",nx);
    fscanf(fin,"ny=%d\n",ny);
    fscanf(fin,"nz=%d\n",nz);
    fscanf(fin,"lt=%d\n",lt);
    fscanf(fin,"nedge=%d\n",nedge);
    fscanf(fin,"ncx_shot1=%d\n",ncx_shot1);
    fscanf(fin,"ncy_shot1=%d\n",ncy_shot1);
    fscanf(fin,"ncz_shot=%d\n",ncz_shot);
    fscanf(fin,"nxshot=%d\n",nxshot);
    fscanf(fin,"nyshot=%d\n",nyshot);
    fscanf(fin,"frequency=%f\n",frequency);
    fscanf(fin,"velmax=%f\n",velmax);
    fscanf(fin,"dt=%f\n",dt);
    fscanf(fin,"unit=%f\n",unit);
    fscanf(fin,"dxshot=%d\n",dxshot);
    fscanf(fin,"dyshot=%d\n",dyshot);
    fclose(fin);
    printf("wanl\n");
}
void print_out(int *nx, int *ny, int *nz, int *lt, int *nedge, int *ncx_shot1, int *ncy_shot1, int *ncz_shot, int *nxshot, int *nyshot, float *frequency, float *velmax, float *dt, float *unit,int *dxshot,int *dyshot)
{
	printf("\n--------workload parameter--------\n");
    printf("nx = %d\n", *nx);
    printf("ny = %d\n", *ny);
    printf("nz = %d\n", *nz);
    printf("lt = %d\n", *lt);
    printf("nedge = %d\n", *nedge);
    printf("ncx_shot1 = %d\n", *ncx_shot1);
    printf("ncy_shot1 = %d\n", *ncy_shot1);
    printf("ncz_shot = %d\n", *ncz_shot);
    printf("nxshot = %d\n", *nxshot);
    printf("nyshot = %d\n", *nyshot);
    printf("frequency = %f\n", *frequency);
    printf("velmax = %f\n", *velmax);
    printf("dt = %f\n", *dt);
    printf("unit = %f\n", *unit);
    printf("dxshot = %d\n", *dxshot);
    printf("dyshot = %d\n\n", *dyshot);
}
void log_end(char *logfile, int *nx, int *ny, int *nz, int *lt, int *nedge, int *ncx_shot1, int *ncy_shot1, int *ncz_shot, int *nxshot, int *nyshot, float *frequency, float *velmax, float *dt, float *unit,int *dxshot,int *dyshot)
{
    FILE *flog;

	flog = fopen(logfile, "a");
    fprintf(flog, "\n--------workload parameter--------\n");
    fprintf(flog, "nx = %d\n", nx);
    fprintf(flog, "ny = %d\n", ny);
    fprintf(flog, "nz = %d\n", nz);
    fprintf(flog, "lt = %d\n", lt);
    fprintf(flog, "nedge = %d\n", nedge);
    fprintf(flog, "ncx_shot1 = %d\n", ncx_shot1);
    fprintf(flog, "ncy_shot1 = %d\n", ncy_shot1);
    fprintf(flog, "ncz_shot = %d\n", ncz_shot);
    fprintf(flog, "nxshot = %d\n", nxshot);
    fprintf(flog, "nyshot = %d\n", nyshot);
    fprintf(flog, "frequency = %f\n", frequency);
    fprintf(flog, "velmax = %f\n",velmax);
    fprintf(flog, "dt = %f\n", dt);
    fprintf(flog, "unit = %f\n", unit);
    fprintf(flog, "dxshot = %d\n", dxshot);
    fprintf(flog, "dyshot = %d\n\n", dyshot);
    fclose(flog);
}

void init (int *nz, int *ny, int *nx, float *vpp, float *vss, float *density)
{
	int i, j, k;

	for ( i = 0; i < *nz; i++ )
		for ( j = 0; j < *ny; j++ )
			for ( k = 0; k < *nx; k++ )
			{
				int yuan;
                yuan = i * *ny * *nx + j * *nx + k;
                
                if ( i < 210 )
                {
                	vpp[ yuan ] = 2300.;
                    vss[ yuan]=1232.;
                    density[ yuan]=1.;    //密度ρ

                   
                }
                else if(i >= 210 && i < 260)
                {
                    vpp[ yuan ] = 2800.;
                    vss[ yuan] = 1500.;
                    density[ yuan] = 2.;   //密度ρ

                   

                }
                else
                {

                    vpp[yuan] = 3500.;
                    vss[yuan] = 1909.;
                    density[yuan] = 2.5;           //密度

                   
                }
            }

}

void formula (int *lt, float *dt, float *tt, float *t0, float *wave, float *frequency)
{
	int l;

	for ( l = 0; l < (*lt); l++ )
    {
        *tt = l * *dt;
        *tt = *tt - *t0;
        float sp = pie * (*frequency) * (*tt);
        float fx = 100000. * exp ( -sp * sp ) * ( 1. - 2. * sp * sp);
        wave[l] = fx;
    }
}
void init_c(int *mm, float c[][7], float *c0)
{
	int i, j;

	if( *mm == 5)
    {
        *c0 = -2.927222164;
        c[0][0] = 1.66666665;
        c[1][0] = -0.23809525;
        c[2][0] = 0.03968254;
        c[3][0] = -0.004960318;
        c[4][0] = 0.0003174603;
    }

    c[0][1] = 0.83333;
    c[1][1] = -0.2381;
    c[2][1] = 0.0595;
    c[3][1] = -0.0099;
    c[4][1] = 0.0008;

    
}
void log_insert(char *logfile, int *ishot)
{
	FILE *flog;

	printf("shot = %d\n", *ishot);
    flog = fopen(logfile, "a");
    fprintf(flog, "shot = %d\n", ishot);
    
    fclose(flog);

}

void set_in (int yuan, float *u, float *v, float *w, float *up, float *up1, float *up2, float *vp, float *vp1, float *vp2, float *wp, float *wp1, float *wp2, float *us, float *us1, float *us2, float *vs, float *vs1, float *vs2, float *ws, float *ws1, float *ws2)
{
    *(u + 1) = 2;

	u[yuan]=0.0f;
 
    v[yuan]=0.0f;
 
    w[yuan]=0.0f;
  
    up[yuan]=0.0f;

    up1[yuan]=0.0f;
   
    up2[yuan]=0.0f;
 
    vp[yuan]=0.0f;
  
    vp1[yuan]=0.0f;
   
    vp2[yuan]=0.0f;
  
    wp[yuan]=0.0f;

    wp1[yuan]=0.0f;
  
    wp2[yuan]=0.0f;
 
    us[yuan]=0.0f;
  
    us1[yuan]=0.0f;

    us2[yuan]=0.0f;
  
    vs[yuan]=0.0f;
  
    vs1[yuan]=0.0f;
 
    vs2[yuan]=0.0f;
 
    ws[yuan]=0.0f;
  
    ws1[yuan]=0.0f;
    
    ws2[yuan]=0.0f;
}
void set_out (int *nx, int *ny, int *nz, float *u, float *v, float *w, float *up, float *up1, float *up2, float *vp, float *vp1, float *vp2, float *wp, float *wp1, float *wp2, float *us, float *us1, float *us2, float *vs, float *vs1, float *vs2, float *ws, float *ws1, float *ws2)
{
	int i, j, k;

	for ( i = 0; i < (*nz); i++ )
		for ( j = 0; j < (*ny); j++ )
			for ( k = 0; k < (*nx); k++ )
			{
				int yuan;

                yuan = i * (*ny) * (*nx) + j * (*nx) + k;
                set_in(yuan, u, v, w, up, up1, up2, vp, vp1, vp2, wp, wp1, wp2, us, us1, us2, vs, vs1, vs2, ws, ws1, ws2);
            }
}
void judge_modify(int *l, float *dt, float *velmax, float *xmax, int *nleft, int *nright, int *nfront, int *nback, int *ntop, int *nbottom, float *unit, int *ncx_shot, int *ncy_shot, int *ncz_shot, int *nx, int *ny, int *nz)
{
	*xmax    = *l * *dt * *velmax;
    *nleft   = *ncx_shot - *xmax / *unit - 10;
    *nright  = *ncx_shot + *xmax / *unit + 10;
    *nfront  = *ncy_shot - *xmax / *unit - 10;
    *nback   = *ncy_shot + *xmax / *unit + 10;
    *ntop    = *ncz_shot - *xmax / *unit - 10;
    *nbottom = *ncz_shot + *xmax / *unit + 10;
            
    if(*nleft < 5)
    	*nleft = 5;
    if(*nright > *nx - 5)
        *nright = *nx - 5;
    if(*nfront < 5) 
        *nfront = 5;
    if(*nback > *ny - 5) 
        *nback = *ny - 5;
    if(*ntop < 5) 
        *ntop = 5;
    if(*nbottom > *nz - 5) 
        *nbottom = *nz - 5;

}

void set_temp2(float *vvp2, int *mm, float c[][7], float *u, float *v, float *w, float *tempux2, float *tempuy2, float *tempuz2, float *tempvx2, float *tempvy2, float *tempvz2, float *tempwx2, float *tempwy2, float *tempwz2, int *nx, int *ny, int *i, int *j, int *k, float *c0, float *vpp2, float *vvs2, float *dtx, float *dtz)
{
	int kk;
	float tmp_d2z;
	float tmp_d2x;
	int yuan = *k * *ny * *nx + *j * *nx + *i;

	for ( kk = 1; kk <= *mm; kk++ )
    {
    	(*tempux2) = *tempux2 + c[kk - 1][0] * (u[yuan + kk] + u[yuan - kk] );
        (*tempuy2) = *tempuy2 + c[kk - 1][0] * (u[yuan + kk * *nx] + u[yuan - kk * *nx] );
        (*tempuz2) = *tempuz2 + c[kk - 1][0] * (u[yuan + kk * *ny] + u[yuan - kk * *ny] );

        (*tempvx2) = *tempvx2 + c[kk - 1][0] * (v[yuan + kk] + v[yuan - kk] );
        (*tempvy2) = *tempvy2 + c[kk - 1][0] * (v[yuan + kk * *nx] + v[yuan - kk * *nx] );
        (*tempvz2) = *tempvz2 + c[kk - 1][0] * (v[yuan + kk * *ny] + v[yuan - kk * *ny] );

        (*tempwx2) = *tempwx2 + c[kk - 1][0] * (w[yuan + kk] + w[yuan - kk] );
        (*tempwy2) = *tempwy2 + c[kk - 1][0] * (w[yuan + kk * *nx] + w[yuan - kk * *nx] );
        (*tempwz2) = *tempwz2 + c[kk - 1][0] * (w[yuan + kk * *ny] + w[yuan - kk * *ny] );

        //temp = temp + 某系数 * [(x + △x) + (x - △x)]

    } //for(kk=1;kk<=mm;kk++) end


    tmp_d2x = *dtx * *dtx;
    tmp_d2z = *dtz * *dtz;

    *tempux2 = (*tempux2 + *c0 * u[yuan]) * *vvp2 * tmp_d2x;
    *tempuy2 = (*tempuy2 + *c0 * u[yuan]) * *vvs2 * tmp_d2x;
    *tempuz2 = (*tempuz2 + *c0 * u[yuan]) * *vvs2 * tmp_d2z;

    *tempvx2 = (*tempvx2 + *c0 * v[yuan]) * *vvs2 * tmp_d2x;
    *tempvy2 = (*tempvy2 + *c0 * v[yuan]) * *vvp2 * tmp_d2x;
    *tempvz2 = (*tempvz2 + *c0 * v[yuan]) * *vvs2 * tmp_d2z;

    *tempwx2 = (*tempwx2 + *c0 * w[yuan]) * *vvs2 * tmp_d2x;
    *tempwy2 = (*tempwy2 + *c0 * w[yuan]) * *vvs2 * tmp_d2x;
    *tempwz2 = (*tempwz2 + *c0 * w[yuan]) * *vvp2 * tmp_d2z;
}
void set_temp(int *mm, float *tempwz2, float *tempuxz, float *tempuxy, float *tempvyz, float *tempvxy, float *tempwxz, float *tempwyz, float *u, float *v, float *w, int *nx, int *ny, int *i, int *j, int *k, float c[][7])
{
	int kk;
	int kkk;
	int yuan = *k * *ny * *nx + *j * *nx + *i;
	//int yuan2;

	for ( kk = 1; kk <= *mm; kk++)
	{
		for( kkk = 1; kkk <= *mm; kkk++)
		{
			*tempuxz = *tempuxz + c[kkk - 1][1 + kk] * (u[yuan + kkk  * *ny * *nx + kk ]
				- u [yuan - kkk  * *ny * *nx + kk ]
				+ u [yuan - kkk  * *ny * *nx - kk]
				- u [yuan + kkk  * *ny * *nx - kk]);
			*tempuxy = *tempuxy + c[kkk - 1][1 + kk] * (u[yuan + kkk * *nx + kk]
				- u[yuan - kkk * *nx + kk ]
				+ u[yuan - kkk * *nx - kk]
				- u[yuan + kkk * *nx - kk]);

			*tempvyz = *tempvyz + c[kkk - 1][1 + kk] * (v[yuan + kkk * *nx * *ny + kk * *nx]
				- v[yuan - kkk * *nx * *ny + kk * *nx]
				+ v[yuan - kkk * *nx * *ny - kk * *nx]
				- v[yuan + kkk * *nx * *ny - kk * *nx]);
			*tempvxy = *tempvxy + c[kkk - 1][1 + kk] * (v[yuan + kkk * *nx + kk]
				- v[yuan - kkk * *nx + kk]
				+ v[yuan - kkk * *nx - kk]
				- v[yuan + kkk * *nx - kk]);

			*tempwyz = *tempwyz + c[kkk - 1][1 + kk] * (w[yuan + kkk * *nx * *ny + kk * *nx]
				- w[yuan - kkk * *nx * *ny + kk * *nx]
				+ w[yuan - kkk * *nx * *ny - kk * *nx]
				- w[yuan + kkk * *nx * *ny - kk * *nx]);
			*tempwxz = *tempwxz+c[kkk - 1][1 + kk]*(w[yuan + kkk  * *ny * *nx + kk ]
				- w[yuan + kkk  * *ny * *nx + kk]
				+ w[yuan - kkk  * *ny * *nx - kk]
				- w[yuan + kkk  * *ny * *nx - kk]);
		} // for(kkk=1;kkk<=mm;kkk++) end
    } //for(kk=1;kk<=mm;kk++) end
}

void set_ps( float *dtx, float *tempwyz, float *px, float *up, float *up1, float *up2, float *vp, float *vp1, float *vp2, float *wp, float *wp1, float *wp2, float *us, float *us1, float *us2, float *vs, float *vs1, float *vs2, float *ws, float *ws1, float *ws2, int *i, int *j, int *k, int *nx, int *ny, float *tempux2, float *tempuy2, float *tempuz2, float *tempvx2, float *tempvy2, float *tempvz2, float *tempwx2, float *tempwy2, float *tempwz2, float *tempuxz, float *tempuxy, float *tempvyz, float *tempvxy, float *tempwxz, float *dtz, float *vvp2, float *vvs2, float *wave, int *l)
{
	int yuan = *k * *ny * *nx + *j * *nx + *i;

	up[yuan] = 2. * up1[yuan] - up2[yuan] + *tempux2 + *tempwxz * *vvp2 * *dtz * *dtx + *tempvxy * *vvp2 * *dtz * *dtx;
	vp[yuan] = 2. * vp1[yuan] - vp2[yuan] + *tempvy2 + *tempuxy * *vvp2 * *dtz * *dtx + *tempwyz * *vvp2 * *dtz * *dtx;
	wp[yuan] = 2. * wp1[yuan] - wp2[yuan] + *tempwz2 + *tempuxz * *vvp2 * *dtz * *dtx + *tempvyz * *vvp2 * *dtz * *dtx + *px * wave[*l-1];
	us[yuan] = 2. * us1[yuan] - us2[yuan] + *tempuy2 + *tempuz2 - *tempvxy * *vvs2 * *dtz * *dtx - *tempwxz * *vvs2 * *dtz * *dtx;
	vs[yuan] = 2. * vs1[yuan] - vs2[yuan] + *tempvx2 + *tempvz2 - *tempuxy * *vvs2 * *dtz * *dtx - *tempwyz * *vvs2 * *dtz * *dtx;
	ws[yuan] = 2. * ws1[yuan] - ws2[yuan] + *tempwx2 + *tempwy2 - *tempuxz * *vvs2 * *dtz * *dtx - *tempvyz * *vvs2 * *dtz * *dtx;
}
void set_ps2(int *ntop, int *nbottom, int *nfront, int *nback, int *nleft, int *nright, int *nx, int *ny, float *u, float *v, float *w, float *up, float *up1, float *up2, float *vp, float *vp1, float *vp2, float *wp, float *wp1, float *wp2, float *us, float *us1, float *us2, float *vs, float *vs1, float *vs2, float *ws, float *ws1, float *ws2)
{
	int i, j, k;
	int yuan;

	for(k = *ntop; k < *nbottom; k++)
		for(j = *nfront; j < *nback; j++)
			for( i = *nleft; i < *nright; i++)
			{
				yuan = k * (*ny) * (*nx) + j * (*nx) + i;

				u[yuan] = up[yuan] + us[yuan];
				v[yuan] = vp[yuan] + vs[yuan];
				w[yuan] = wp[yuan] + ws[yuan];
				up2[yuan] = up1[yuan];
				up1[yuan] = up[yuan];
				us2[yuan] = us1[yuan];
				us1[yuan] = us[yuan];
				vp2[yuan] = vp1[yuan];
				vp1[yuan] = vp[yuan];
				vs2[yuan] = vs1[yuan];
				vs1[yuan] = vs[yuan];
				wp2[yuan] = wp1[yuan];
				wp1[yuan] = wp[yuan];
				ws2[yuan] = ws1[yuan];
				ws1[yuan] = ws[yuan];
            }//for(i=nleft;i<nright;i++) end
}
void get_result(float *vpp2, float *c0, float c[][7], float *tempwz2, float *vpp, float *drd1, float *dr1, float *drd2, float * dr2, float *vss, char *outfile, char *logfile, float *nshot, int *ncy_shot, int *ncy_shot1, int *nxshot, int *dyshot, int *ncx_shot, int *ncx_shot1, int *dxshot, int *nx, int *ny, int *nz, float *u, float *v, float *w, float *up,float *up1, float *up2, float *vp, float *vp1, float *vp2, float *wp, float *wp1, float *wp2, float *us, float *us1, float *us2, float *vs, float *vs1, float *vs2, float *ws, float *ws1, float *ws2, float *tempux2, float *tempuy2, float *tempuz2, float *tempvx2, float *tempvy2, float *tempvz2,
           float *tempwx2, float *tempwy2, float *tempuxz, float *tempuxy, float *tempvyz, float *tempvxy, float *tempwxz, float *tempwyz, float *dtx, float *dtz, float *vvp2, float *vvs2, float *wave, int *lt, int *mm,  float *dt, float *velmax, float *xmax, int *nleft, int *nright, int *nfront, int *nback, int *ntop, int*nbottom, float *unit, int *ncz_shot, float *px, float *sx)
{
	FILE *fout;

	int l;
	int ishot;
    int i, j , k;

	fout = fopen(outfile, "wb");

    for ( ishot = 1; ishot <= *nshot; ishot++)
    {
        log_insert(logfile, &ishot);   //根据log_insert可看出下面两句的时间
        
        (*ncy_shot) = (*ncy_shot1) + ( ishot / (*nxshot)) * (*dyshot);
        (*ncx_shot) = (*ncx_shot1) + ( ishot % (*nxshot)) * (*dxshot);

        //set_out (int *nx, int *ny, int *nz, float *u, float *v, float *w, float *up,float *up1, float *up2, float *vp, float *vp1, float *vp2, float *wp, float *wp1, float *wp2, float *us, float *us1, float *us2, float *vs, float *vs1, float *vs2, float *ws, float *ws1, float *ws2)

        set_out(nx, ny, nz, u, v, w, up, up1, up2, vp, vp1, vp2, wp, wp1, wp2, us, us1, us2, vs, vs1, vs2, ws, ws1, ws2);

        for( l = 1; l <= *lt; l++)
        {
            judge_modify(&l, dt, velmax, xmax, nleft, nright, nfront, nback, ntop, nbottom, unit, ncx_shot, ncy_shot, ncz_shot, nx, ny, nz);
           // printf("chu judge_modify\n");
            (*ntop) = *ntop - 1;
            *nfront = *nfront - 1;
            *nleft = *nleft - 1;

            for ( k = *ntop; k < *nbottom; k++ )
                for ( j = *nfront; j < *nback; j++ )
                    for ( i = *nleft; i < *nright; i++ )
                    {
                        int yuan;

                        yuan = k * *ny * *nx + j * *nx + i;

                        if( i == *ncx_shot - 1 && j == *ncy_shot - 1 && k == *ncz_shot - 1)
                        {
                            *px = 1.;
                            *sx = 0.;
                        }
                        else
                        {
                            *px = 0.;
                            *sx = 0.;
                        }
                        *vvp2 = vpp[ yuan ] * vpp[ yuan ];
                        *drd1 = (*dr1) * *vvp2;
                        *drd2 = (*dr2) * *vvp2;

                        *vvs2 = vss[ yuan ] * vss[ yuan ];
                        *drd1 = (*dr1) * *vvs2;
                        *drd2 = (*dr2) * *vvs2;

                        *tempux2 = 0.0f;
                        *tempuy2 = 0.0f;
                        *tempuz2 = 0.0f;
                        *tempvx2 = 0.0f;
                        *tempvy2 = 0.0f;
                        *tempvz2 = 0.0f;
                        *tempwx2 = 0.0f;
                        *tempwy2 = 0.0f;
                        *tempwz2 = 0.0f;
                        *tempuxz = 0.0f;
                        *tempuxy = 0.0f;
                        *tempvyz = 0.0f;
                        *tempvxy = 0.0f;
                        *tempwxz = 0.0f;
                        *tempwyz = 0.0f;
                        set_temp2(vvp2, mm, c , u, v, w, tempux2, tempuy2, tempuz2, tempvx2, tempvy2, tempvz2, tempwx2, tempwy2, tempwz2, nx, ny, &i, &j, &k, c0, vpp2, vvs2,  dtx, dtz);
                       // printf("chu set_temp2\n");
                        set_temp(mm, tempwz2, tempuxz, tempuxy, tempvyz, tempvxy, tempwxz, tempwyz, u, v, w, nx, ny, &i, &j, &k, c );
                     //   printf("chu set_temp\n");
                        set_ps(dtx, tempwyz, px, up, up1, up2, vp, vp1, vp2, wp, wp1, wp2, us, us1, us2, vs, vs1, vs2, ws, ws1, ws2, &i, &j, &k, nx, ny, tempux2, tempuy2, tempuz2, tempvx2, tempvy2, tempvz2, tempwx2,  tempwy2, tempwz2, tempuxz, tempuxy, tempvyz, tempvxy, tempwxz, dtz, vvp2, vvs2, wave, &l);
                   //     printf("chu set_ps\n");
                        
                    }
            set_ps2(ntop, nbottom, nfront, nback, nleft, nright, nx, ny, u, v, w, up, up1, up2, vp, vp1, vp2, wp, wp1, wp2, us, us1, us2, vs, vs1, vs2, ws, ws1, ws2);
        }
        
        fwrite(up + 169 * (*ny) * (*nx),sizeof(float),*ny * *nx, fout);

    }
    fclose(fout);
}
int main(int argc, char **argv)
{
	int mm = 5;
	int nx,ny,nz,lt,nedge;
    int nleft,nright,nfront,nback,ntop,nbottom;    
    int ishot,ncy_shot,ncx_shot;
    int nxshot,nyshot,dxshot,dyshot;
    int ncx_shot1,ncy_shot1,ncz_shot;
    float unit, dt;
    float dtx,dtz,dtxz,dr1, drd1, dr2,drd2, dtx4,dtz4,dtxz4;
    float *u, *v, *w, *up, *up1, *up2, 
          *vp, *vp1, *vp2, *wp, *wp1, *wp2, 
          *us, *us1, *us2, *vs, *vs1, *vs2,
          *ws, *ws1, *ws2, *vpp, *density, *vss;
    float *vpp2;
    float vvp2,vvs2,tempux2,tempuy2,tempuz2,tempvx2,tempvy2,tempvz2,
          tempwx2,tempwy2,tempwz2,tempuxz,tempuxy,tempvyz,tempvxy,tempwxz,tempwyz;
	struct timeval start, end;
    float all_time;
    float c[5][7];
    float *wave;
    float velmax;
    float frequency;
      float xmax,px,sx;
    float nshot,t0,tt,c0;
	char infile[80],outfile[80],logfile[80];
    FILE *flog;
    int tmp_len;

	if(argc < 4)
    {
        printf("please add 3 parameter: inpurfile, outfile, logfile\n");
        exit(0);
    }
    analysis(infile, outfile, logfile, argv);
    log_start(logfile, strlen(logfile));
    gettimeofday(&start, NULL);

    read_data(infile, &nx, &ny, &nz, &lt, &nedge, &ncx_shot1, &ncy_shot1, &ncz_shot, &nxshot, &nyshot, &frequency, &velmax, &dt, &unit, &dxshot, &dyshot);
    print_out(&nx, &ny, &nz, &lt, &nedge, &ncx_shot1, &ncy_shot1, &ncz_shot, &nxshot, &nyshot, &frequency, &velmax, &dt, &unit,&dxshot,&dyshot);

printf("chu print_out\n");
 //   apply(u, v, w, up,up1, up2, vp, vp1, vp2, wp, wp1, wp2, us, us1, us2, vs, vs1, vs2, ws, ws1, ws2, vpp, density, vss, wave, , &lt);


 tmp_len= sizeof(float) * (nz*ny*nx);

    u       = (float*)malloc(tmp_len);
    v       = (float*)malloc(tmp_len);
    w       = (float*)malloc(tmp_len);
    up      = (float*)malloc(tmp_len);
    up1     = (float*)malloc(tmp_len);
    up2     = (float*)malloc(tmp_len);
    vp      = (float*)malloc(tmp_len);
    vp1     = (float*)malloc(tmp_len);
    vp2     = (float*)malloc(tmp_len);
    wp      = (float*)malloc(tmp_len);
    wp1     = (float*)malloc(tmp_len);
    wp2     = (float*)malloc(tmp_len);
    us      = (float*)malloc(tmp_len);
    us1     = (float*)malloc(tmp_len);
    us2     = (float*)malloc(tmp_len);
    vs      = (float*)malloc(tmp_len);
    vs1     = (float*)malloc(tmp_len);
    vs2     = (float*)malloc(tmp_len);
    ws      = (float*)malloc(tmp_len);
    ws1     = (float*)malloc(tmp_len);
    ws2     = (float*)malloc(tmp_len);
    vpp     = (float*)malloc(tmp_len);
    density = (float*)malloc(tmp_len);
    vss     = (float*)malloc(tmp_len);
    wave    = (float*)malloc(sizeof(float) * (lt));

    nshot = nxshot * nyshot;

    t0 = 1.0 / frequency;    //周期?

    formula (&lt, &dt, &tt, &t0, wave, &frequency);



    init_c(&mm, c, &c0);

for ( i = 0; i < 5; i++)
        for ( j = 0; j < 5; j++)
            c[j][2 + i] = c[i][1] * c[j][1];

    dtx = dt / unit;
    dtz = dt / unit;
    dtxz = dtx * dtz;

    dr1 = dtx * dtx / 2. ;
    dr2 = dtz * dtz / 2. ;

    dtx4 = dtx * dtx * dtx * dtx;
    dtz4 = dtz * dtz * dtz * dtz;
    dtxz4 = dtx * dtx * dtz * dtz;

  
    get_result(vpp2, &c0, c, &tempwz2, vpp, &drd1, &dr1, &drd2, &dr2, vss, outfile, logfile, &nshot, &ncy_shot, &ncy_shot1, &nxshot , &dyshot, &ncx_shot, &ncx_shot1,  &dxshot, &nx, &ny,  &nz,  u, v, w, up, up1, up2, vp, vp1, vp2, wp, wp1, wp2,  us, us1, us2,  vs, vs1, vs2, ws, ws1, ws2, &tempux2, &tempuy2, &tempuz2, &tempvx2,  &tempvy2, &tempvz2,
            &tempwx2, &tempwy2, &tempuxz, &tempuxy, &tempvyz, &tempvxy, &tempwxz, &tempwyz, &dtx, &dtz, &vvp2, &vvs2, wave, &lt, &mm, &dt,  &velmax, &xmax, &nleft, &nright, &nfront, &nback, &ntop, &nbottom, &unit, &ncz_shot, &px, &sx);


    free(u);
    free(v);
    free(w);
    free(up);
    free(up1);
    free(up2);
    free(vp);
    free(vp1);
    free(vp2);
    free(wp);
    free(wp1);
    free(wp2);
    free(us);
    free(us1);
    free(us2);
    free(vs);
    free(vs1);
    free(vs2);
    free(ws);
    free(ws1);
    free(ws2);
    free(vpp);
    free(density);
    free(vss);
    free(wave);

    gettimeofday(&end,NULL);
    all_time = (end.tv_sec-start.tv_sec)+(float)(end.tv_usec-start.tv_usec)/1000000.0;
    printf("run time:\t%f s\n",all_time);
    flog = fopen(logfile,"a");
    fprintf(flog,"\nrun time:\t%f s\n\n",all_time);
    fclose(flog);
    flog = fopen(logfile,"a");
    fprintf(flog,"------------end time------------\n");
    fclose(flog);
//    system(tmp);
    return 1;
}