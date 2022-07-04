// Adapted from Prof. Koyama's MPF code in his textbook
// Author: Chuanqi Zhu
// Created on: 2022/2/16
// Updated on 2022/06/20

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "mpi.h"

#define NDX 50 //差分計算における計算領域一辺の分割数
#define NDY 50 //差分計算における計算領域一辺の分割数
#define NDZ 50
#define N 3 //考慮する結晶方位の数＋１(MPF0.cppと比較して、この値を大きくしている)
#define BEGIN 1
#define UTAG 2
#define DTAG 3
#define NONE 0
#define DONE 4
#define MASTER 0

int ndx = NDX;
int ndy = NDY;
int ndz = NDZ;
int ndmx = NDX - 1;
int ndmy = NDY - 1; //計算領域の一辺の差分分割数(差分ブロック数), ND-1を定義
int ndmz = NDZ - 1;
int nm = N - 1, nmm = N - 2; //考慮する結晶方位の数、N-2（考慮する結晶方位の数－１）を定義
double PI = 3.141592;        //π、計算カウント数
double RR = 8.3145;          //ガス定数

double aij[N][N]; //勾配エネルギー係数
double wij[N][N];
double mij[N][N];
double fij[N][N];
double thij[N][N];
double vpij[N][N];
double etaij[N][N];
double face[16][3];
int phinum;
double th, vp, eta;

int i, j, k, l, ii, jj, kk, ll, it; //整数
int ip, im, jp, jm, kp, km;         //整数
int n1, n2, n3;                     //整数

int istep = 0;
// int n000;		//位置(i,j)において、pが０ではない方位の個数（n00>=n000）
int nstep;               //計算カウント数の最大値（計算終了カウント）
double dtime, L, dx;     // L計算領域の一辺の長さ(nm), 差分プロック１辺の長さ(m)
double M0;               //粒界の易動度
double W0;               //ペナルティー項の係数
double A0;               //勾配エネルギー係数
double F0;               //粒界移動の駆動力
double temp;             //温度
double sum1, sum2, sum3; //各種の和の作業変数
double pddtt;            //フェーズフィールドの時間変化率
double phidxii, phidyii, phidzii, phiabsii;
double nxii, nyii, nzii, alii;
double al111, alm111, al1m11, al11m1;
double alm511, al1m51, al11m5;
double al511, al15m1, al1m15, al51m1, al151, alm115, al5m11, alm151, al115;
double miijj;
double min_val, rp0, rp1;
int min_idx;
double zeta1, zeta2, zeta3;

double phidx, phidy, phidz, phiabs, phiabs2;
double phidxx, phidxy, phidxz;
double phidyx, phidyy, phidyz;
double phidzx, phidzy, phidzz;

double xxp, xyp, xzp;
double yxp, yyp, yzp;
double zxp, zyp, zzp;

double dphiabs2dx, dphiabs2dy, dphiabs2dz;
double del, al0, al, alm, am;
double nx, ny, nz;
double nxx, nxy, nxz;
double nyx, nyy, nyz;
double nzx, nzy, nzz;
double ux, uy, uz, uu;

double nxphix, nyphix, nzphix;
double nxphiy, nyphiy, nzphiy;
double nxphiz, nyphiz, nzphiz;
double nxphixdx, nyphixdx, nzphixdx;
double nxphiydy, nyphiydy, nzphiydy;
double nxphizdz, nyphizdz, nzphizdz;

double PP, QQ, CC, SS;
double dPdx, dPdy, dPdz;
double dQdx, dQdy, dQdz;
double dPdphix, dPdphiy, dPdphiz;
double dQdphix, dQdphiy, dQdphiz;
double dPdphixdx, dPdphiydy, dPdphizdz;
double dQdphixdx, dQdphiydy, dQdphizdz;

double epsilon0, ep;
double epdx, epdy, epdz;
double epdphix, epdphiy, epdphiz;
double epdphixdx, epdphiydy, epdphizdz;

double termx, termy, termz;
double termiikk, termjjkk;

double gamma0; //粒界エネルギ密度
double delta;  //粒界幅（差分ブロック数にて表現）
double mobi;   //粒界の易動度
double vm0;    //モル体積

int x11, y11, z11, x1h[10], y1h[10], z1h[10]; //初期核の座標
double t, r0, r;

//******* メインプログラム ******************************************
int main(int argc, char *argv[])
{
    nstep = 2000;
    dtime = 1.0;
    temp = 1000.0;
    L = 2000.0;
    vm0 = 7.0e-6;
    delta = 7.0;
    mobi = 1.0;
    zeta1 = 0.001;
    zeta2 = 0.6;
    zeta3 = 0.8;
    rp0 = 0.02;
    rp1 = 0.05;
    al0 = 35.0 / 180.0 * PI;

    dx = L / 100 * 1.0e-9;               //差分プロック１辺の長さ(m)
    gamma0 = 0.5 * vm0 / RR / temp / dx; //粒界エネルギ密度（0.5J/m^2）を無次元化
    A0 = 8.0 * delta * gamma0 / PI / PI; //勾配エネルギー係数[式(4.40)]
    W0 = 4.0 * gamma0 / delta;           //ペナルティー項の係数[式(4.40)]
    M0 = mobi * PI * PI / (8.0 * delta); //粒界の易動度[式(4.40)]
    F0 = 100.0 / RR / temp;              //粒界移動の駆動力

    for (ii = 1; ii <= nm; ii++)
    {
        for (jj = 1; jj <= nm; jj++)
        {
            wij[ii][jj] = W0;
            aij[ii][jj] = A0;
            mij[ii][jj] = M0;
            fij[ii][jj] = 0.0;
            thij[ii][jj] = 0.0;
            vpij[ii][jj] = 0.0;
            etaij[ii][jj] = 0.0;
            if ((ii == nm) || (jj == nm))
            {
                fij[ii][jj] = F0;
            }
            if (ii > jj)
            {
                fij[ii][jj] = -fij[ii][jj];
            }
            if (ii == jj)
            {
                wij[ii][jj] = 0.0;
                aij[ii][jj] = 0.0;
                mij[ii][jj] = 0.0;
                fij[ii][jj] = 0.0;
            }
        }
    }

    // thij[1][2] = PI / 4.0;
    // thij[2][1] = PI / 4.0;

    face[0][0] = 1.0;
    face[0][1] = 1.0;
    face[0][2] = 1.0;

    face[1][0] = -1.0;
    face[1][1] = 1.0;
    face[1][2] = 1.0;

    face[2][0] = 1.0;
    face[2][1] = -1.0;
    face[2][2] = 1.0;

    face[3][0] = 1.0;
    face[3][1] = 1.0;
    face[3][2] = -1.0;

    face[4][0] = -5.0;
    face[4][1] = 1.0;
    face[4][2] = 1.0;

    face[5][0] = 1.0;
    face[5][1] = -5.0;
    face[5][2] = 1.0;

    face[6][0] = 1.0;
    face[6][1] = 1.0;
    face[6][2] = -5.0;

    face[7][0] = 5.0;
    face[7][1] = 1.0;
    face[7][2] = 1.0;

    face[8][0] = 1.0;
    face[8][1] = 5.0;
    face[8][2] = -1.0;

    face[9][0] = 1.0;
    face[9][1] = -1.0;
    face[9][2] = 5.0;

    face[10][0] = 5.0;
    face[10][1] = 1.0;
    face[10][2] = -1.0;

    face[11][0] = 1.0;
    face[11][1] = 5.0;
    face[11][2] = 1.0;

    face[12][0] = 1.0;
    face[12][1] = 1.0;
    face[12][2] = 5.0;

    face[13][0] = 5.0;
    face[13][1] = -1.0;
    face[13][2] = 1.0;

    face[14][0] = -1.0;
    face[14][1] = 5.0;
    face[14][2] = 1.0;

    face[15][0] = 1.0;
    face[15][1] = 1.0;
    face[15][2] = 5.0;

    int taskid,
        numworkers,
        numtasks,
        rows, offset,
        dest, source,
        up, down,
        msgtype,
        rc, start, end,
        ix, iy, iz, it;

    MPI_Status status;

    // Allocate taskid to each core (cores = tasks = master + workers)
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);

    numworkers = numtasks - 1;
    rows = NDX / numworkers;

    /************************* master code *******************************/
    if (taskid == MASTER)
    {
        clock_t start_t, end_t, total_t;
        start_t = clock();
        if (NDX % numworkers != 0)
        {
            MPI_Abort(MPI_COMM_WORLD, rc);
            exit(1);
        }

        double(*phi)[N][NDX][NDY][NDZ] = malloc(sizeof(*phi));
        double(*intphi)[NDX][NDY][NDZ] = malloc(sizeof(*intphi));

        for (i = 0; i <= ndmx; i++)
        {
            for (j = 0; j <= ndmy; j++)
            {
                for (k = 0; k <= ndmz; k++)
                {
                    for (ii = 1; ii <= nm - 1; ii++)
                    {
                        (*phi)[ii][i][j][k] = 0.0;
                    }
                    (*phi)[nm][i][j][k] = 1.0; // nm番目のフェーズフィールドを１に初期化
                }
            }
        }

        for (i = 0; i <= ndmx; i++)
        {
            for (j = 0; j <= ndmy; j++)
            {
                for (k = 0; k <= ndmy; k++)
                {
                    r = sqrt(((i - NDX / 2)) * (i - NDX / 2) + (j - NDY / 2) * (j - NDY / 2) + (k - NDZ / 2) * (k - NDZ / 2));
                    if (r <= NDX / 10)
                    {
                        (*phi)[1][i][j][k] = 1.0;
                        (*phi)[2][i][j][k] = 0.0;
                    } //初期核位置のフェーズフィールドを設定
                }
            }
        }

        offset = 0;
        // Send to workers
        for (i = 1; i <= numworkers; i++)
        {
            dest = i;
            if (dest == 1)
                up = NONE;
            else
                up = dest - 1;
            if (dest == numworkers)
                down = NONE;
            else
                down = dest + 1;

            MPI_Send(&offset, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
            MPI_Send(&rows, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
            MPI_Send(&up, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
            MPI_Send(&down, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
            //// send phase fields
            for (ii = 1; ii <= nm; ii++)
            {
                MPI_Send(&(*phi)[ii][offset], rows * NDY * NDZ, MPI_DOUBLE, dest, BEGIN, MPI_COMM_WORLD);
            }

            offset = offset + rows;
        }

        // Receive from workers
        for (i = 1; i <= numworkers; i++)
        {
            source = i;
            msgtype = DONE;
            MPI_Recv(&offset, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD,
                     &status);
            MPI_Recv(&rows, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
            //// receive phase fields
            MPI_Recv(&(*intphi)[offset], rows * NDY * NDZ, MPI_DOUBLE, source,
                     msgtype, MPI_COMM_WORLD, &status);
        }

        FILE *stream;
        char buffer[30];
        sprintf(buffer, "data/2d%d.vtk", nstep);
        stream = fopen(buffer, "a");

        fprintf(stream, "# vtk DataFile Version 1.0\n");
        fprintf(stream, "phi_%d.vtk\n", nstep);
        fprintf(stream, "ASCII\n");
        fprintf(stream, "DATASET STRUCTURED_POINTS\n");
        fprintf(stream, "DIMENSIONS %d %d %d\n", NDX, NDY, NDZ);
        fprintf(stream, "ORIGIN 0.0 0.0 0.0\n");
        fprintf(stream, "ASPECT_RATIO 1.0 1.0 1.0\n");
        fprintf(stream, "\n");
        fprintf(stream, "POINT_DATA %d\n", NDX * NDY * NDZ);
        fprintf(stream, "SCALARS scalars float\n");
        fprintf(stream, "LOOKUP_TABLE default\n");

        for (k = 0; k <= ndmz; k++)
        {
            for (j = 0; j <= ndmy; j++)
            {
                for (i = 0; i <= ndmx; i++)
                {
                    // fprintf(streamc0, "%e\n", phi[1][i][j][k]);
                    fprintf(stream, "%e\n", (*intphi)[i][j][k]);
                }
            }
        }
        fclose(stream);

        end_t = clock();
        total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
        printf("Total time taken: %lu secs\n", total_t);
        printf("Exiting of the program...\n");

        MPI_Finalize();
    }

    /************************* workers code **********************************/
    if (taskid != MASTER)
    {
        double(*phi)[N][rows + 2][NDY][NDZ] = malloc(sizeof(*phi));
        double(*phi2)[N][rows + 2][NDY][NDZ] = malloc(sizeof(*phi2));
        int(*phiNum)[rows + 2][NDY][NDZ] = malloc(sizeof(*phiNum));
        int(*phiIdx)[N + 1][rows + 2][NDY][NDZ] = malloc(sizeof(*phiIdx));
        double(*intphi)[rows + 2][NDY][NDZ] = malloc(sizeof(*intphi));

        // Receive from master
        source = MASTER;
        msgtype = BEGIN;
        MPI_Recv(&offset, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&rows, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&up, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&down, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
        //// receive phase fields
        for (ii = 1; ii <= nm; ii++)
        {
            MPI_Recv(&(*phi)[ii][1], rows * NDY * NDZ, MPI_DOUBLE, source, msgtype, MPI_COMM_WORLD, &status);
        }

    start:;

        // Communicate with neighor workers before computation
        // (The size of each message should be limited under 8500 floating point number for Mac M1 chip)
        if (up != NONE)
        {
            //// send up boundaries of phase fields
            for (ii = 1; ii <= nm; ii++)
            {
                MPI_Send(&(*phi)[ii][1][0], NDY * NDZ / 2, MPI_DOUBLE, up, DTAG, MPI_COMM_WORLD);
                MPI_Send(&(*phi)[ii][1][NDY / 2], NDY * NDZ / 2, MPI_DOUBLE, up, DTAG, MPI_COMM_WORLD);
            }

            source = up;
            msgtype = UTAG;
            //// receive up boundaries of phase fields
            for (ii = 1; ii <= nm; ii++)
            {
                MPI_Recv(&(*phi)[ii][0][0], NDY * NDZ / 2, MPI_DOUBLE, source,
                         msgtype, MPI_COMM_WORLD, &status);
                MPI_Recv(&(*phi)[ii][0][NDY / 2], NDY * NDZ / 2, MPI_DOUBLE, source,
                         msgtype, MPI_COMM_WORLD, &status);
            }
        }
        if (down != NONE)
        {
            //// send down boundaries of phase fields
            for (ii = 1; ii <= nm; ii++)
            {
                MPI_Send(&(*phi)[ii][rows][0], NDY * NDZ / 2, MPI_DOUBLE, down,
                         UTAG, MPI_COMM_WORLD);
                MPI_Send(&(*phi)[ii][rows][NDY / 2], NDY * NDZ / 2, MPI_DOUBLE, down,
                         UTAG, MPI_COMM_WORLD);
            }

            source = down;
            msgtype = DTAG;
            //// receive down boundaries of phase fields
            for (ii = 1; ii <= nm; ii++)
            {
                MPI_Recv(&(*phi)[ii][rows + 1][0], NDY * NDZ / 2, MPI_DOUBLE, source, msgtype,
                         MPI_COMM_WORLD, &status);
                MPI_Recv(&(*phi)[ii][rows + 1][NDY / 2], NDY * NDZ / 2, MPI_DOUBLE, source, msgtype,
                         MPI_COMM_WORLD, &status);
            }
        }

        // Compute after sending and receiving data
        start = 1;
        end = rows;

        for (i = start; i <= end; i++)
        {
            for (j = 0; j <= ndmy; j++)
            {
                for (k = 0; k <= ndmz; k++)
                {
                    ip = i + 1;
                    im = i - 1;
                    jp = j + 1;
                    jm = j - 1;
                    kp = k + 1;
                    km = k - 1;
                    if (up == NONE && i == 1)
                    {
                        im = 1;
                    }
                    if (down == NONE && i == rows)
                    {
                        ip = rows;
                    }
                    if (j == ndmy)
                    {
                        jp = 0;
                    }
                    if (j == 0)
                    {
                        jm = ndmy;
                    }
                    if (k == ndmz)
                    {
                        kp = 0;
                    }
                    if (k == 0)
                    {
                        km = ndmz;
                    }

                    //--- 位置(i,j)およびその周囲(i±1,j±1)において、pが０ではない方位の個数---
                    phinum = 0;
                    for (ii = 1; ii <= nm; ii++)
                    {
                        if (((*phi)[ii][i][j][k] > 0.0) ||
                            (((*phi)[ii][i][j][k] == 0.0) && ((*phi)[ii][ip][j][k] > 0.0) ||
                             ((*phi)[ii][im][j][k] > 0.0) ||
                             ((*phi)[ii][i][jp][k] > 0.0) ||
                             ((*phi)[ii][i][jm][k] > 0.0) ||
                             ((*phi)[ii][i][j][kp] > 0.0) ||
                             ((*phi)[ii][i][j][km] > 0.0)))
                        {
                            phinum++;
                            (*phiIdx)[phinum][i][j][k] = ii;
                            // printf("%d  ", n00);
                        }
                    }
                    (*phiNum)[i][j][k] = phinum;
                }
            }
        }

        // Evolution Equations
        for (i = start; i <= end; i++)
        {
            for (j = 0; j <= ndmy; j++)
            {
                for (k = 0; k <= ndmz; k++)
                {
                    ip = i + 1;
                    im = i - 1;
                    jp = j + 1;
                    jm = j - 1;
                    kp = k + 1;
                    km = k - 1;
                    if (up == NONE && i == 1)
                    {
                        im = 1;
                    }
                    if (down == NONE && i == rows)
                    {
                        ip = rows;
                    }
                    if (j == ndmy)
                    {
                        jp = 0;
                    }
                    if (j == 0)
                    {
                        jm = ndmy;
                    }
                    if (k == ndmz)
                    {
                        kp = 0;
                    }
                    if (k == 0)
                    {
                        km = ndmz;
                    }

                    for (n1 = 1; n1 <= (*phiNum)[i][j][k]; n1++)
                    {
                        ii = (*phiIdx)[n1][i][j][k];
                        pddtt = 0.0;

                        phidxii = ((*phi)[ii][ip][j][k] - (*phi)[ii][im][j][k]) / 2.0 / dx;
                        phidyii = ((*phi)[ii][i][jp][k] - (*phi)[ii][i][jm][k]) / 2.0 / dx;
                        phidzii = ((*phi)[ii][i][j][kp] - (*phi)[ii][i][j][km]) / 2.0 / dx;
                        phiabsii = phidxii * phidxii + phidyii * phidyii + phidzii * phidzii;

                        // thii = thij[ii][jj];
                        // vpii = vpij[ii][jj];
                        // etaii = etaij[ii][jj];

                        // xxpii = cos(th) * cos(vp);
                        // yxpii = sin(th) * cos(vp);
                        // zxpii = sin(vp);
                        // xypii = -sin(th) * cos(eta) - cos(th) * sin(vp) * sin(eta);
                        // yypii = cos(th) * cos(eta) - sin(th) * sin(vp) * sin(eta);
                        // zypii = cos(vp) * sin(eta);
                        // xzpii = sin(eta) * sin(th) - cos(eta) * cos(th) * sin(vp);
                        // yzpii = -sin(eta) * cos(th) - cos(eta) * sin(th) * sin(vp);
                        // zzpii = cos(eta) * cos(vp);

                        for (n2 = 1; n2 <= (*phiNum)[i][j][k]; n2++)
                        {
                            jj = (*phiIdx)[n2][i][j][k];
                            sum1 = 0.0;
                            for (n3 = 1; n3 <= (*phiNum)[i][j][k]; n3++)
                            {
                                kk = (*phiIdx)[n3][i][j][k];

                                // calculate the interface normal and deirivatives of the phase field

                                phidxx = ((*phi)[kk][ip][j][k] + (*phi)[kk][im][j][k] - 2.0 * (*phi)[kk][i][j][k]);
                                phidxy = ((*phi)[kk][ip][jp][k] + (*phi)[kk][im][jm][k] - (*phi)[kk][im][jp][k] - (*phi)[kk][ip][jm][k]) / 4.0;
                                phidxz = ((*phi)[kk][ip][j][kp] + (*phi)[kk][im][j][km] - (*phi)[kk][im][j][kp] - (*phi)[kk][ip][j][km]) / 4.0;

                                phidyx = phidxy;
                                phidyz = ((*phi)[kk][i][jp][kp] + (*phi)[kk][i][jm][km] - (*phi)[kk][i][jm][kp] - (*phi)[kk][i][jp][km]) / 4.0;
                                phidyy = ((*phi)[kk][i][jp][k] + (*phi)[kk][i][jm][k] - 2.0 * (*phi)[kk][i][j][k]);

                                phidzx = phidxz;
                                phidzy = phidyz;
                                phidzz = ((*phi)[kk][i][j][kp] + (*phi)[kk][i][j][km] - 2.0 * (*phi)[kk][i][j][k]);

                                phidx = ((*phi)[kk][ip][j][k] - (*phi)[kk][im][j][k]) / 2.0;
                                phidy = ((*phi)[kk][i][jp][k] - (*phi)[kk][i][jm][k]) / 2.0;
                                phidz = ((*phi)[kk][i][j][kp] - (*phi)[kk][i][j][km]) / 2.0;

                                phiabs2 = phidx * phidx + phidy * phidy + phidz * phidz;
                                phiabs = sqrt(phiabs2);

                                dphiabs2dx = 2.0 * (phidx * phidxx + phidy * phidxy + phidz * phidxz);
                                dphiabs2dy = 2.0 * (phidx * phidxy + phidy * phidyy + phidz * phidyz);
                                dphiabs2dz = 2.0 * (phidx * phidxz + phidy * phidyz + phidz * phidzz);

                                del = 5.0;

                                am = (1.0 + del * sqrt(1.0 + 2.0 * rp0 * rp0) * (1.0 + tan(al0) * tan(al0)));

                                if ((ii + kk == 3))
                                {
                                    th = thij[ii][kk];
                                    vp = vpij[ii][kk];
                                    eta = etaij[ii][kk];

                                    xxp = cos(th) * cos(vp);
                                    yxp = sin(th) * cos(vp);
                                    zxp = sin(vp);
                                    xyp = -sin(th) * cos(eta) - cos(th) * sin(vp) * sin(eta);
                                    yyp = cos(th) * cos(eta) - sin(th) * sin(vp) * sin(eta);
                                    zyp = cos(vp) * sin(eta);
                                    xzp = sin(eta) * sin(th) - cos(eta) * cos(th) * sin(vp);
                                    yzp = -sin(eta) * cos(th) - cos(eta) * sin(th) * sin(vp);
                                    zzp = cos(eta) * cos(vp);

                                    nx = phidx / phiabs;
                                    ny = phidy / phiabs;
                                    nz = phidz / phiabs;

                                    for (l = 0; l <= 3; l++)
                                    {
                                        ux = face[l][0];
                                        uy = face[l][1];
                                        uz = face[l][2];
                                        uu = ux * ux + uy * uy + uz * uz;
                                        al = acos(fabs(nx * ux + ny * uy + nz * uz) / sqrt(uu));
                                        if (l == 0)
                                        {
                                            min_val = al;
                                            min_idx = l;
                                        }
                                        else
                                        {
                                            if (min_val > al)
                                            {
                                                min_val = al;
                                                min_idx = l;
                                            }
                                        }
                                    }

                                    ux = face[min_idx][0];
                                    uy = face[min_idx][1];
                                    uz = face[min_idx][2];
                                    uu = ux * ux + uy * uy + uz * uz;
                                    al = acos((nx * ux + ny * uy + nz * uz) / sqrt(uu));

                                    alm = acos(sqrt((1.0 + rp0 * rp0 * (1.0 - tan(al0) * tan(al0))) / (1.0 + tan(al0) * tan(al0))));

                                    if (fabs(cos(al)) >= cos(alm))
                                    {
                                        epsilon0 = sqrt(aij[ii][kk]);

                                        PP = nx * ny * (ux * uy) + ny * nz * (uy * uz) + nx * nz * (ux * uz);
                                        QQ = pow(nx * ux, 2.0) + pow(ny * uy, 2.0) + pow(nz * uz, 2.0);

                                        CC = sqrt((2.0 * PP + QQ) / uu + rp0 * rp0);
                                        SS = sqrt((uu - QQ - 2.0 * PP) / uu + rp0 * rp0);

                                        ep = epsilon0 * (1.0 + del * (CC + tan(al0) * SS)) / am;

                                        nxx = phidxx / phiabs - phidx / phiabs / phiabs2 * dphiabs2dx / 2.0;
                                        nyx = phidyx / phiabs - phidy / phiabs / phiabs2 * dphiabs2dx / 2.0;
                                        nzx = phidzx / phiabs - phidz / phiabs / phiabs2 * dphiabs2dx / 2.0;

                                        nxy = phidxy / phiabs - phidx / phiabs / phiabs2 * dphiabs2dy / 2.0;
                                        nyy = phidyy / phiabs - phidy / phiabs / phiabs2 * dphiabs2dy / 2.0;
                                        nzy = phidzy / phiabs - phidz / phiabs / phiabs2 * dphiabs2dy / 2.0;

                                        nxz = phidxz / phiabs - phidx / phiabs / phiabs2 * dphiabs2dz / 2.0;
                                        nyz = phidyz / phiabs - phidy / phiabs / phiabs2 * dphiabs2dz / 2.0;
                                        nzz = phidzz / phiabs - phidz / phiabs / phiabs2 * dphiabs2dz / 2.0;

                                        nxphix = -phidx * phidx / phiabs / phiabs2 + 1.0 / phiabs;
                                        nyphix = -phidx * phidy / phiabs / phiabs2;
                                        nzphix = -phidx * phidz / phiabs / phiabs2;

                                        nxphiy = -phidy * phidx / phiabs / phiabs2;
                                        nyphiy = -phidy * phidy / phiabs / phiabs2 + 1.0 / phiabs;
                                        nzphiy = -phidy * phidz / phiabs / phiabs2;

                                        nxphiz = -phidz * phidx / phiabs / phiabs2;
                                        nyphiz = -phidz * phidy / phiabs / phiabs2;
                                        nzphiz = -phidz * phidz / phiabs / phiabs2 + 1.0 / phiabs;

                                        nxphixdx = -1.0 / pow(phiabs2, 3.0) * ((phidxx * phidx + phidx * phidxx) * pow(phiabs, 3.0) - 1.5 * phidx * phidx * phiabs * dphiabs2dx) - 0.5 / phiabs / phiabs2 * dphiabs2dx;
                                        nyphixdx = -1.0 / pow(phiabs2, 3.0) * ((phidxx * phidy + phidx * phidyx) * pow(phiabs, 3.0) - 1.5 * phidx * phidy * phiabs * dphiabs2dx);
                                        nzphixdx = -1.0 / pow(phiabs2, 3.0) * ((phidxx * phidz + phidx * phidzx) * pow(phiabs, 3.0) - 1.5 * phidx * phidz * phiabs * dphiabs2dx);

                                        nxphiydy = -1.0 / pow(phiabs2, 3.0) * ((phidyy * phidx + phidy * phidxy) * pow(phiabs, 3.0) - 1.5 * phidy * phidx * phiabs * dphiabs2dy);
                                        nyphiydy = -1.0 / pow(phiabs2, 3.0) * ((phidyy * phidy + phidy * phidyy) * pow(phiabs, 3.0) - 1.5 * phidy * phidy * phiabs * dphiabs2dy) - 0.5 / phiabs / phiabs2 * dphiabs2dy;
                                        nzphiydy = -1.0 / pow(phiabs2, 3.0) * ((phidyy * phidz + phidy * phidzy) * pow(phiabs, 3.0) - 1.5 * phidy * phidz * phiabs * dphiabs2dy);

                                        nxphizdz = -1.0 / pow(phiabs2, 3.0) * ((phidzz * phidx + phidz * phidxz) * pow(phiabs, 3.0) - 1.5 * phidz * phidx * phiabs * dphiabs2dz);
                                        nyphizdz = -1.0 / pow(phiabs2, 3.0) * ((phidzz * phidy + phidz * phidyz) * pow(phiabs, 3.0) - 1.5 * phidz * phidy * phiabs * dphiabs2dz);
                                        nzphizdz = -1.0 / pow(phiabs2, 3.0) * ((phidzz * phidz + phidz * phidzz) * pow(phiabs, 3.0) - 1.5 * phidz * phidz * phiabs * dphiabs2dz) - 0.5 / phiabs / phiabs2 * dphiabs2dz;

                                        dPdx = nxx * ny * (ux * uy) + nx * nyx * (ux * uy) + nyx * nz * (uy * uz) + ny * nzx * (uy * uz) + nxx * nz * (ux * uz) + nx * nzx * (ux * uz);
                                        dPdy = nxy * ny * (ux * uy) + nx * nyy * (ux * uy) + nyy * nz * (uy * uz) + ny * nzy * (uy * uz) + nxy * nz * (ux * uz) + nx * nzy * (ux * uz);
                                        dPdz = nxz * ny * (ux * uy) + nx * nyz * (ux * uy) + nyz * nz * (uy * uz) + ny * nzz * (uy * uz) + nxz * nz * (ux * uz) + nx * nzz * (ux * uz);

                                        dPdphix = nxphix * ny * (ux * uy) + nx * nyphix * (ux * uy) + nyphix * nz * (uy * uz) + ny * nzphix * (uy * uz) + nxphix * nz * (ux * uz) + nx * nzphix * (ux * uz);
                                        dPdphiy = nxphiy * ny * (ux * uy) + nx * nyphiy * (ux * uy) + nyphiy * nz * (uy * uz) + ny * nzphiy * (uy * uz) + nxphiy * nz * (ux * uz) + nx * nzphiy * (ux * uz);
                                        dPdphiz = nxphiz * ny * (ux * uy) + nx * nyphiz * (ux * uy) + nyphiz * nz * (uy * uz) + ny * nzphiz * (uy * uz) + nxphiz * nz * (ux * uz) + nx * nzphiz * (ux * uz);

                                        dPdphixdx = nxphixdx * ny * (ux * uy) + nxphix * nyx * (ux * uy) + nxx * nyphix * (ux * uy) + nx * nyphixdx * (ux * uy) + nyphixdx * nz * (uy * uz) + nyphix * nzx * (uy * uz) + nyx * nzphix * (uy * uz) + ny * nzphixdx * (uy * uz) + nxphixdx * nz * (ux * uz) + nxphix * nzx * (ux * uz) + nxx * nzphix * (ux * uz) + nx * nzphixdx * (ux * uz);
                                        dPdphiydy = nxphiydy * ny * (ux * uy) + nxphiy * nyy * (ux * uy) + nxy * nyphiy * (ux * uy) + nx * nyphiydy * (ux * uy) + nyphiydy * nz * (uy * uz) + nyphiy * nzy * (uy * uz) + nyy * nzphiy * (uy * uz) + ny * nzphiydy * (uy * uz) + nxphiydy * nz * (ux * uz) + nxphiy * nzy * (ux * uz) + nxy * nzphiy * (ux * uz) + nx * nzphiydy * (ux * uz);
                                        dPdphizdz = nxphizdz * ny * (ux * uy) + nxphiz * nyz * (ux * uy) + nxz * nyphiz * (ux * uy) + nx * nyphizdz * (ux * uy) + nyphizdz * nz * (uy * uz) + nyphiz * nzz * (uy * uz) + nyz * nzphiz * (uy * uz) + ny * nzphizdz * (uy * uz) + nxphizdz * nz * (ux * uz) + nxphiz * nzz * (ux * uz) + nxz * nzphiz * (ux * uz) + nx * nzphizdz * (ux * uz);

                                        dQdx = 2.0 * (ux * ux * nx * nxx + uy * uy * ny * nyx + uz * uz * nz * nzx);
                                        dQdy = 2.0 * (ux * ux * nx * nxy + uy * uy * ny * nyy + uz * uz * nz * nzy);
                                        dQdz = 2.0 * (ux * ux * nx * nxz + uy * uy * ny * nyz + uz * uz * nz * nzz);

                                        dQdphix = 2.0 * (ux * ux * nx * nxphix + uy * uy * ny * nyphix + uz * uz * nz * nzphix);
                                        dQdphiy = 2.0 * (ux * ux * nx * nxphiy + uy * uy * ny * nyphiy + uz * uz * nz * nzphiy);
                                        dQdphiz = 2.0 * (ux * ux * nx * nxphiz + uy * uy * ny * nyphiz + uz * uz * nz * nzphiz);

                                        dQdphixdx = 2.0 * (ux * ux * (nxx * nxphix + nx * nxphixdx) + uy * uy * (nyx * nyphix + ny * nyphixdx) + uz * uz * (nzx * nzphix + nz * nzphixdx));
                                        dQdphiydy = 2.0 * (ux * ux * (nxy * nxphiy + nx * nxphiydy) + uy * uy * (nyy * nyphiy + ny * nyphiydy) + uz * uz * (nzy * nzphiy + nz * nzphiydy));
                                        dQdphizdz = 2.0 * (ux * ux * (nxz * nxphiz + nx * nxphizdz) + uy * uy * (nyz * nyphiz + ny * nyphizdz) + uz * uz * (nzz * nzphiz + nz * nzphizdz));

                                        epdx = epsilon0 * del * (1.0 / (2.0 * uu * CC) - tan(al0) / (2.0 * uu * SS)) * (2.0 * dPdx + dQdx) / am;
                                        epdy = epsilon0 * del * (1.0 / (2.0 * uu * CC) - tan(al0) / (2.0 * uu * SS)) * (2.0 * dPdy + dQdy) / am;
                                        epdz = epsilon0 * del * (1.0 / (2.0 * uu * CC) - tan(al0) / (2.0 * uu * SS)) * (2.0 * dPdz + dQdz) / am;

                                        epdphix = epsilon0 * del * (1.0 / (2.0 * uu * CC) - tan(al0) / (2.0 * uu * SS)) * (2.0 * dPdphix + dQdphix) / am;
                                        epdphiy = epsilon0 * del * (1.0 / (2.0 * uu * CC) - tan(al0) / (2.0 * uu * SS)) * (2.0 * dPdphiy + dQdphiy) / am;
                                        epdphiz = epsilon0 * del * (1.0 / (2.0 * uu * CC) - tan(al0) / (2.0 * uu * SS)) * (2.0 * dPdphiz + dQdphiz) / am;

                                        epdphixdx = 1.0 / 2.0 / uu * epsilon0 * del * ((-1.0 / 2.0 / uu / pow(CC, 3.0) - tan(al0) / (2.0 * uu * pow(SS, 3.0))) * (2.0 * dPdx + dQdx) * (2.0 * dPdphix + dQdphix) + (2.0 * dPdphixdx + dQdphixdx) * (1.0 / CC - tan(al0) / SS)) / am;
                                        epdphiydy = 1.0 / 2.0 / uu * epsilon0 * del * ((-1.0 / 2.0 / uu / pow(CC, 3.0) - tan(al0) / (2.0 * uu * pow(SS, 3.0))) * (2.0 * dPdy + dQdy) * (2.0 * dPdphiy + dQdphiy) + (2.0 * dPdphiydy + dQdphiydy) * (1.0 / CC - tan(al0) / SS)) / am;
                                        epdphizdz = 1.0 / 2.0 / uu * epsilon0 * del * ((-1.0 / 2.0 / uu / pow(CC, 3.0) - tan(al0) / (2.0 * uu * pow(SS, 3.0))) * (2.0 * dPdz + dQdz) * (2.0 * dPdphiz + dQdphiz) + (2.0 * dPdphizdz + dQdphizdz) * (1.0 / CC - tan(al0) / SS)) / am;

                                        termx = ep * ep * phidxx + 2.0 * ep * epdx * phidx + ep * epdx * dphiabs2dx + ep * epdphixdx * phiabs2 + epdx * epdphix * phiabs2;
                                        termy = ep * ep * phidyy + 2.0 * ep * epdy * phidy + ep * epdy * dphiabs2dy + ep * epdphiydy * phiabs2 + epdy * epdphiy * phiabs2;
                                        termz = ep * ep * phidzz + 2.0 * ep * epdz * phidz + ep * epdz * dphiabs2dz + ep * epdphizdz * phiabs2 + epdz * epdphiz * phiabs2;

                                        termiikk = termx + termy + termz;
                                    }
                                    else
                                    {
                                        epsilon0 = sqrt(aij[ii][kk]);
                                        ep = epsilon0 * (1.0 + del * (sqrt(cos(alm) * cos(alm) + rp0 * rp0) + tan(al0) * sqrt(sin(alm) * sin(alm) + rp0 * rp0))) / am;
                                        termiikk = ep * ep * (phidxx + phidyy + phidzz);
                                    }
                                }
                                else
                                {
                                    epsilon0 = sqrt(aij[ii][kk]);
                                    ep = epsilon0 * (1.0 + del * (sqrt(cos(alm) * cos(alm) + rp0 * rp0) + tan(al0) * sqrt(sin(alm) * sin(alm) + rp0 * rp0))) / am;
                                    termiikk = ep * ep * (phidxx + phidyy + phidzz);
                                }

                                if ((jj + kk == 3))
                                {
                                    th = thij[ii][kk];
                                    vp = vpij[ii][kk];
                                    eta = etaij[ii][kk];

                                    xxp = cos(th) * cos(vp);
                                    yxp = sin(th) * cos(vp);
                                    zxp = sin(vp);
                                    xyp = -sin(th) * cos(eta) - cos(th) * sin(vp) * sin(eta);
                                    yyp = cos(th) * cos(eta) - sin(th) * sin(vp) * sin(eta);
                                    zyp = cos(vp) * sin(eta);
                                    xzp = sin(eta) * sin(th) - cos(eta) * cos(th) * sin(vp);
                                    yzp = -sin(eta) * cos(th) - cos(eta) * sin(th) * sin(vp);
                                    zzp = cos(eta) * cos(vp);

                                    nx = phidx / phiabs;
                                    ny = phidy / phiabs;
                                    nz = phidz / phiabs;

                                    for (l = 0; l <= 3; l++)
                                    {
                                        ux = face[l][0];
                                        uy = face[l][1];
                                        uz = face[l][2];
                                        uu = ux * ux + uy * uy + uz * uz;
                                        al = acos(fabs(nx * ux + ny * uy + nz * uz) / sqrt(uu));
                                        if (l == 0)
                                        {
                                            min_val = al;
                                            min_idx = l;
                                        }
                                        else
                                        {
                                            if (min_val > al)
                                            {
                                                min_val = al;
                                                min_idx = l;
                                            }
                                        }
                                    }

                                    ux = face[min_idx][0];
                                    uy = face[min_idx][1];
                                    uz = face[min_idx][2];
                                    uu = ux * ux + uy * uy + uz * uz;
                                    al = acos((nx * ux + ny * uy + nz * uz) / sqrt(uu));

                                    alm = acos(sqrt((1.0 + rp0 * rp0 * (1.0 - tan(al0) * tan(al0))) / (1.0 + tan(al0) * tan(al0))));

                                    if (fabs(cos(al)) >= cos(alm))
                                    {
                                        epsilon0 = sqrt(aij[jj][kk]);

                                        PP = nx * ny * (ux * uy) + ny * nz * (uy * uz) + nx * nz * (ux * uz);
                                        QQ = pow(nx * ux, 2.0) + pow(ny * uy, 2.0) + pow(nz * uz, 2.0);

                                        CC = sqrt((2.0 * PP + QQ) / uu + rp0 * rp0);
                                        SS = sqrt((uu - QQ - 2.0 * PP) / uu + rp0 * rp0);

                                        ep = epsilon0 * (1.0 + del * (CC + tan(al0) * SS)) / am;

                                        nxx = phidxx / phiabs - phidx / phiabs / phiabs2 * dphiabs2dx / 2.0;
                                        nyx = phidyx / phiabs - phidy / phiabs / phiabs2 * dphiabs2dx / 2.0;
                                        nzx = phidzx / phiabs - phidz / phiabs / phiabs2 * dphiabs2dx / 2.0;

                                        nxy = phidxy / phiabs - phidx / phiabs / phiabs2 * dphiabs2dy / 2.0;
                                        nyy = phidyy / phiabs - phidy / phiabs / phiabs2 * dphiabs2dy / 2.0;
                                        nzy = phidzy / phiabs - phidz / phiabs / phiabs2 * dphiabs2dy / 2.0;

                                        nxz = phidxz / phiabs - phidx / phiabs / phiabs2 * dphiabs2dz / 2.0;
                                        nyz = phidyz / phiabs - phidy / phiabs / phiabs2 * dphiabs2dz / 2.0;
                                        nzz = phidzz / phiabs - phidz / phiabs / phiabs2 * dphiabs2dz / 2.0;

                                        nxphix = -phidx * phidx / phiabs / phiabs2 + 1.0 / phiabs;
                                        nyphix = -phidx * phidy / phiabs / phiabs2;
                                        nzphix = -phidx * phidz / phiabs / phiabs2;

                                        nxphiy = -phidy * phidx / phiabs / phiabs2;
                                        nyphiy = -phidy * phidy / phiabs / phiabs2 + 1.0 / phiabs;
                                        nzphiy = -phidy * phidz / phiabs / phiabs2;

                                        nxphiz = -phidz * phidx / phiabs / phiabs2;
                                        nyphiz = -phidz * phidy / phiabs / phiabs2;
                                        nzphiz = -phidz * phidz / phiabs / phiabs2 + 1.0 / phiabs;

                                        nxphixdx = -1.0 / pow(phiabs2, 3.0) * ((phidxx * phidx + phidx * phidxx) * pow(phiabs, 3.0) - 1.5 * phidx * phidx * phiabs * dphiabs2dx) - 0.5 / phiabs / phiabs2 * dphiabs2dx;
                                        nyphixdx = -1.0 / pow(phiabs2, 3.0) * ((phidxx * phidy + phidx * phidyx) * pow(phiabs, 3.0) - 1.5 * phidx * phidy * phiabs * dphiabs2dx);
                                        nzphixdx = -1.0 / pow(phiabs2, 3.0) * ((phidxx * phidz + phidx * phidzx) * pow(phiabs, 3.0) - 1.5 * phidx * phidz * phiabs * dphiabs2dx);

                                        nxphiydy = -1.0 / pow(phiabs2, 3.0) * ((phidyy * phidx + phidy * phidxy) * pow(phiabs, 3.0) - 1.5 * phidy * phidx * phiabs * dphiabs2dy);
                                        nyphiydy = -1.0 / pow(phiabs2, 3.0) * ((phidyy * phidy + phidy * phidyy) * pow(phiabs, 3.0) - 1.5 * phidy * phidy * phiabs * dphiabs2dy) - 0.5 / phiabs / phiabs2 * dphiabs2dy;
                                        nzphiydy = -1.0 / pow(phiabs2, 3.0) * ((phidyy * phidz + phidy * phidzy) * pow(phiabs, 3.0) - 1.5 * phidy * phidz * phiabs * dphiabs2dy);

                                        nxphizdz = -1.0 / pow(phiabs2, 3.0) * ((phidzz * phidx + phidz * phidxz) * pow(phiabs, 3.0) - 1.5 * phidz * phidx * phiabs * dphiabs2dz);
                                        nyphizdz = -1.0 / pow(phiabs2, 3.0) * ((phidzz * phidy + phidz * phidyz) * pow(phiabs, 3.0) - 1.5 * phidz * phidy * phiabs * dphiabs2dz);
                                        nzphizdz = -1.0 / pow(phiabs2, 3.0) * ((phidzz * phidz + phidz * phidzz) * pow(phiabs, 3.0) - 1.5 * phidz * phidz * phiabs * dphiabs2dz) - 0.5 / phiabs / phiabs2 * dphiabs2dz;

                                        dPdx = nxx * ny * (ux * uy) + nx * nyx * (ux * uy) + nyx * nz * (uy * uz) + ny * nzx * (uy * uz) + nxx * nz * (ux * uz) + nx * nzx * (ux * uz);
                                        dPdy = nxy * ny * (ux * uy) + nx * nyy * (ux * uy) + nyy * nz * (uy * uz) + ny * nzy * (uy * uz) + nxy * nz * (ux * uz) + nx * nzy * (ux * uz);
                                        dPdz = nxz * ny * (ux * uy) + nx * nyz * (ux * uy) + nyz * nz * (uy * uz) + ny * nzz * (uy * uz) + nxz * nz * (ux * uz) + nx * nzz * (ux * uz);

                                        dPdphix = nxphix * ny * (ux * uy) + nx * nyphix * (ux * uy) + nyphix * nz * (uy * uz) + ny * nzphix * (uy * uz) + nxphix * nz * (ux * uz) + nx * nzphix * (ux * uz);
                                        dPdphiy = nxphiy * ny * (ux * uy) + nx * nyphiy * (ux * uy) + nyphiy * nz * (uy * uz) + ny * nzphiy * (uy * uz) + nxphiy * nz * (ux * uz) + nx * nzphiy * (ux * uz);
                                        dPdphiz = nxphiz * ny * (ux * uy) + nx * nyphiz * (ux * uy) + nyphiz * nz * (uy * uz) + ny * nzphiz * (uy * uz) + nxphiz * nz * (ux * uz) + nx * nzphiz * (ux * uz);

                                        dPdphixdx = nxphixdx * ny * (ux * uy) + nxphix * nyx * (ux * uy) + nxx * nyphix * (ux * uy) + nx * nyphixdx * (ux * uy) + nyphixdx * nz * (uy * uz) + nyphix * nzx * (uy * uz) + nyx * nzphix * (uy * uz) + ny * nzphixdx * (uy * uz) + nxphixdx * nz * (ux * uz) + nxphix * nzx * (ux * uz) + nxx * nzphix * (ux * uz) + nx * nzphixdx * (ux * uz);
                                        dPdphiydy = nxphiydy * ny * (ux * uy) + nxphiy * nyy * (ux * uy) + nxy * nyphiy * (ux * uy) + nx * nyphiydy * (ux * uy) + nyphiydy * nz * (uy * uz) + nyphiy * nzy * (uy * uz) + nyy * nzphiy * (uy * uz) + ny * nzphiydy * (uy * uz) + nxphiydy * nz * (ux * uz) + nxphiy * nzy * (ux * uz) + nxy * nzphiy * (ux * uz) + nx * nzphiydy * (ux * uz);
                                        dPdphizdz = nxphizdz * ny * (ux * uy) + nxphiz * nyz * (ux * uy) + nxz * nyphiz * (ux * uy) + nx * nyphizdz * (ux * uy) + nyphizdz * nz * (uy * uz) + nyphiz * nzz * (uy * uz) + nyz * nzphiz * (uy * uz) + ny * nzphizdz * (uy * uz) + nxphizdz * nz * (ux * uz) + nxphiz * nzz * (ux * uz) + nxz * nzphiz * (ux * uz) + nx * nzphizdz * (ux * uz);

                                        dQdx = 2.0 * (ux * ux * nx * nxx + uy * uy * ny * nyx + uz * uz * nz * nzx);
                                        dQdy = 2.0 * (ux * ux * nx * nxy + uy * uy * ny * nyy + uz * uz * nz * nzy);
                                        dQdz = 2.0 * (ux * ux * nx * nxz + uy * uy * ny * nyz + uz * uz * nz * nzz);

                                        dQdphix = 2.0 * (ux * ux * nx * nxphix + uy * uy * ny * nyphix + uz * uz * nz * nzphix);
                                        dQdphiy = 2.0 * (ux * ux * nx * nxphiy + uy * uy * ny * nyphiy + uz * uz * nz * nzphiy);
                                        dQdphiz = 2.0 * (ux * ux * nx * nxphiz + uy * uy * ny * nyphiz + uz * uz * nz * nzphiz);

                                        dQdphixdx = 2.0 * (ux * ux * (nxx * nxphix + nx * nxphixdx) + uy * uy * (nyx * nyphix + ny * nyphixdx) + uz * uz * (nzx * nzphix + nz * nzphixdx));
                                        dQdphiydy = 2.0 * (ux * ux * (nxy * nxphiy + nx * nxphiydy) + uy * uy * (nyy * nyphiy + ny * nyphiydy) + uz * uz * (nzy * nzphiy + nz * nzphiydy));
                                        dQdphizdz = 2.0 * (ux * ux * (nxz * nxphiz + nx * nxphizdz) + uy * uy * (nyz * nyphiz + ny * nyphizdz) + uz * uz * (nzz * nzphiz + nz * nzphizdz));

                                        epdx = epsilon0 * del * (1.0 / (2.0 * uu * CC) - tan(al0) / (2.0 * uu * SS)) * (2.0 * dPdx + dQdx) / am;
                                        epdy = epsilon0 * del * (1.0 / (2.0 * uu * CC) - tan(al0) / (2.0 * uu * SS)) * (2.0 * dPdy + dQdy) / am;
                                        epdz = epsilon0 * del * (1.0 / (2.0 * uu * CC) - tan(al0) / (2.0 * uu * SS)) * (2.0 * dPdz + dQdz) / am;

                                        epdphix = epsilon0 * del * (1.0 / (2.0 * uu * CC) - tan(al0) / (2.0 * uu * SS)) * (2.0 * dPdphix + dQdphix) / am;
                                        epdphiy = epsilon0 * del * (1.0 / (2.0 * uu * CC) - tan(al0) / (2.0 * uu * SS)) * (2.0 * dPdphiy + dQdphiy) / am;
                                        epdphiz = epsilon0 * del * (1.0 / (2.0 * uu * CC) - tan(al0) / (2.0 * uu * SS)) * (2.0 * dPdphiz + dQdphiz) / am;

                                        epdphixdx = 1.0 / 2.0 / uu * epsilon0 * del * ((-1.0 / 2.0 / uu / pow(CC, 3.0) - tan(al0) / (2.0 * uu * pow(SS, 3.0))) * (2.0 * dPdx + dQdx) * (2.0 * dPdphix + dQdphix) + (2.0 * dPdphixdx + dQdphixdx) * (1.0 / CC - tan(al0) / SS)) / am;
                                        epdphiydy = 1.0 / 2.0 / uu * epsilon0 * del * ((-1.0 / 2.0 / uu / pow(CC, 3.0) - tan(al0) / (2.0 * uu * pow(SS, 3.0))) * (2.0 * dPdy + dQdy) * (2.0 * dPdphiy + dQdphiy) + (2.0 * dPdphiydy + dQdphiydy) * (1.0 / CC - tan(al0) / SS)) / am;
                                        epdphizdz = 1.0 / 2.0 / uu * epsilon0 * del * ((-1.0 / 2.0 / uu / pow(CC, 3.0) - tan(al0) / (2.0 * uu * pow(SS, 3.0))) * (2.0 * dPdz + dQdz) * (2.0 * dPdphiz + dQdphiz) + (2.0 * dPdphizdz + dQdphizdz) * (1.0 / CC - tan(al0) / SS)) / am;

                                        termx = ep * ep * phidxx + 2.0 * ep * epdx * phidx + ep * epdx * dphiabs2dx + ep * epdphixdx * phiabs2 + epdx * epdphix * phiabs2;
                                        termy = ep * ep * phidyy + 2.0 * ep * epdy * phidy + ep * epdy * dphiabs2dy + ep * epdphiydy * phiabs2 + epdy * epdphiy * phiabs2;
                                        termz = ep * ep * phidzz + 2.0 * ep * epdz * phidz + ep * epdz * dphiabs2dz + ep * epdphizdz * phiabs2 + epdz * epdphiz * phiabs2;

                                        termjjkk = termx + termy + termz;
                                    }
                                    else
                                    {
                                        epsilon0 = sqrt(aij[jj][kk]);
                                        ep = epsilon0 * (1.0 + del * (sqrt(cos(alm) * cos(alm) + rp0 * rp0) + tan(al0) * sqrt(sin(alm) * sin(alm) + rp0 * rp0))) / am;
                                        termjjkk = ep * ep * (phidxx + phidyy + phidzz);
                                    }
                                }
                                else
                                {
                                    epsilon0 = sqrt(aij[jj][kk]);
                                    ep = epsilon0 * (1.0 + del * (sqrt(cos(alm) * cos(alm) + rp0 * rp0) + tan(al0) * sqrt(sin(alm) * sin(alm) + rp0 * rp0))) / am;
                                    termjjkk = ep * ep * (phidxx + phidyy + phidzz);
                                }

                                sum1 += 0.5 * (termiikk - termjjkk) + (wij[ii][kk] - wij[jj][kk]) * (*phi)[kk][i][j][k];
                                // sum1 += 0.5 * (aij[ii][kk] - aij[jj][kk]) * ((*phi)[kk][ip][j][k] + (*phi)[kk][im][j][k] + (*phi)[kk][i][jp][k] + (*phi)[kk][i][jm][k] + (*phi)[kk][i][j][kp] + (*phi)[kk][i][j][km] - 6.0 * (*phi)[kk][i][j][k]) + (wij[ii][kk] - wij[jj][kk]) * (*phi)[kk][i][j][k]; //[式(4.31)の一部]
                            }
                            // if ((ii + jj) == 3)
                            // {
                            //     nxii = phidxii / sqrt(phiabsii);
                            //     nyii = phidyii / sqrt(phiabsii);
                            //     nzii = phidzii / sqrt(phiabsii);

                            //     al111 = acos(fabs(nxii + nyii + nzii) / sqrt(3.0));
                            //     alm111 = acos(fabs(-nxii + nyii + nzii) / sqrt(3.0));
                            //     al1m11 = acos(fabs(nxii - nyii + nzii) / sqrt(3.0));
                            //     al11m1 = acos(fabs(nxii + nyii - nzii) / sqrt(3.0));

                            //     alm511 = acos(fabs(-5.0 * nxii + nyii + nzii) / sqrt(27.0));
                            //     al1m51 = acos(fabs(nxii - 5.0 * nyii + nzii) / sqrt(27.0));
                            //     al11m5 = acos(fabs(nxii + nyii - 5.0 * nzii) / sqrt(27.0));

                            //     al511 = acos(fabs(5.0 * nxii + nyii + nzii) / sqrt(27.0));
                            //     al15m1 = acos(fabs(nxii + 5.0 * nyii - nzii) / sqrt(27.0));
                            //     al1m15 = acos(fabs(nxii - nyii + 5.0 * nzii) / sqrt(27.0));
                            //     al51m1 = acos(fabs(5.0 * nxii + nyii - nzii) / sqrt(27.0));
                            //     al151 = acos(fabs(nxii + 5.0 * nyii + nzii) / sqrt(27.0));
                            //     alm115 = acos(fabs(-nxii + nyii + 5.0 * nzii) / sqrt(27.0));
                            //     al5m11 = acos(fabs(5.0 * nxii - nyii + nzii) / sqrt(27.0));
                            //     alm151 = acos(fabs(-nxii + 5.0 * nyii + nzii) / sqrt(27.0));
                            //     al115 = acos(fabs(nxii + nyii + 5.0 * nzii) / sqrt(27.0));

                            //     double arr[16];
                            //     arr[0] = al111;
                            //     arr[1] = alm111;
                            //     arr[2] = al1m11;
                            //     arr[3] = al11m1;
                            //     arr[4] = alm511;
                            //     arr[5] = al1m51;
                            //     arr[6] = al11m5;
                            //     arr[7] = al511;
                            //     arr[8] = al15m1;
                            //     arr[9] = al1m15;
                            //     arr[10] = al51m1;
                            //     arr[11] = al151;
                            //     arr[12] = alm115;
                            //     arr[13] = al5m11;
                            //     arr[14] = alm151;
                            //     arr[15] = al115;

                            //     min_val = arr[0];
                            //     for (l = 1; l <= 15; l++)
                            //     {
                            //         if (min_val > arr[l])
                            //         {
                            //             min_val = arr[l];
                            //         }
                            //     }

                            //     if (min_val == al111)
                            //     {
                            //         miijj = mij[ii][jj] * (zeta1 + (1 - zeta1) * sqrt(pow(tan(min_val), 2.0) + rp1 * rp1) * tanh(1.0 / sqrt(pow(tan(min_val), 2.0) + rp1 * rp1)));
                            //     }
                            //     if ((min_val == alm111) || (min_val == al1m11) || (min_val == al11m1) || (min_val == alm511) || (min_val == al1m51) || (min_val == al11m5))
                            //     {
                            //         miijj = mij[ii][jj] * (zeta2 + (1 - zeta2) * sqrt(pow(tan(min_val), 2.0) + rp1 * rp1) * tanh(1.0 / sqrt(pow(tan(min_val), 2.0) + rp1 * rp1)));
                            //     }
                            //     if ((min_val == al511) || (min_val == al15m1) || (min_val == al1m15) || (min_val == al51m1) || (min_val == al151) || (min_val == alm115) || (min_val == al5m11) || (min_val == alm151) || (min_val == al115))
                            //     {
                            //         miijj = mij[ii][jj] * (zeta3 + (1 - zeta3) * sqrt(pow(tan(min_val), 2.0) + rp1 * rp1) * tanh(1.0 / sqrt(pow(tan(min_val), 2.0) + rp1 * rp1)));
                            //     }
                            // }
                            // else
                            // {
                            //     miijj = mij[ii][jj];
                            // }
                            miijj = mij[ii][jj];
                            pddtt += -2.0 * miijj / (double)((*phiNum)[i][j][k]) * (sum1 - 8.0 / PI * fij[ii][jj] * sqrt((*phi)[ii][i][j][k] * (*phi)[jj][i][j][k]));
                            //フェーズフィールドの発展方程式[式(4.31)]
                        }
                        (*phi2)[ii][i][j][k] = (*phi)[ii][i][j][k] + pddtt * dtime; //フェーズフィールドの時間発展（陽解法）
                        if ((*phi2)[ii][i][j][k] >= 1.0)
                        {
                            (*phi2)[ii][i][j][k] = 1.0;
                        } //フェーズフィールドの変域補正
                        if ((*phi2)[ii][i][j][k] <= 0.0)
                        {
                            (*phi2)[ii][i][j][k] = 0.0;
                        }
                    }
                } // k
            }     // j
        }         // i

        for (ii = 1; ii <= nm; ii++)
        {
            for (i = start; i <= end; i++)
            {
                for (j = 0; j <= ndmy; j++)
                {
                    for (k = 0; k <= ndmz; k++)
                    {
                        (*phi)[ii][i][j][k] = (*phi2)[ii][i][j][k];
                    }
                }
            }
        }

        //
        for (i = start; i <= end; i++)
        {
            for (j = 0; j <= ndmy; j++)
            {
                for (k = 0; k <= ndmz; k++)
                {
                    sum1 = 0.0;
                    for (ii = 1; ii <= nm; ii++)
                    {
                        sum1 += (*phi)[ii][i][j][k];
                    }
                    for (ii = 1; ii <= nm; ii++)
                    {
                        (*phi)[ii][i][j][k] = (*phi)[ii][i][j][k] / sum1;
                    }
                }
            }
        }

        istep = istep + 1.;
        if (istep < nstep)
        {
            goto start;
        }
    end:;

        for (i = start; i <= end; i++)
        {
            for (j = 0; j <= ndmy; j++)
            {
                for (k = 0; k <= ndmz; k++)
                {
                    for (ii = 1; ii <= nm; ii++)
                    {
                        (*intphi)[i][j][k] += (*phi)[ii][i][j][k] * (*phi)[ii][i][j][k];
                    }
                }
            }
        }
        // Send final result to master
        MPI_Send(&offset, 1, MPI_INT, MASTER, DONE, MPI_COMM_WORLD);
        MPI_Send(&rows, 1, MPI_INT, MASTER, DONE, MPI_COMM_WORLD);
        //// send phase fields
        MPI_Send(&(*intphi)[1], rows * NDY * NDZ, MPI_DOUBLE, MASTER, DONE,
                 MPI_COMM_WORLD);
        MPI_Finalize();
    }
    return 0;
}
